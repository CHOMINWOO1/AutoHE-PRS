from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
import tempfile
import unittest

from autohe_prs.benchmark.dataset import build_dataset
from autohe_prs.benchmark.storage import read_records, write_records
from autohe_prs.benchmark.synthetic import apply_missing_policy, make_synthetic_workload
from autohe_prs.cli import main
from autohe_prs.config import Constraints, FHEPlan
from autohe_prs.fhe.emitter import emit
from autohe_prs.fhe.runner import run_openfhe_plan
from autohe_prs.models.train import CostModelBundle, train_cost_models
from autohe_prs.models.evaluate import evaluate_generalization_splits
from autohe_prs.optimize.bayesian import suggest
from autohe_prs.optimize.baselines import run_baseline_search
from autohe_prs.optimize.ablation import apply_plan_ablation, settings
from autohe_prs.optimize.candidate import CandidateSpace, generate_candidates
from autohe_prs.optimize.pareto import pareto_frontier
from autohe_prs.optimize.rank import rank_candidates
from autohe_prs.optimize.tuner import tune
from autohe_prs.reporting.tables import write_prediction_comparison_table
from genome_he_pilot.openfhe_controlled_backend import (
    _reduce_slots,
    _serialize_key_material,
)


def measured_rows() -> list[dict]:
    rows = []
    for group in range(4):
        for index, scheme in enumerate(("CKKS", "BFV", "BGV")):
            variants = 512 * (group + 1)
            success = not (group == 3 and scheme == "BGV")
            rows.append({
                "workload_group": f"g{group}",
                "status": "ok" if success else "failed",
                "scheme": scheme,
                "samples": 1 + group,
                "variants": variants,
                "evaluation_seconds": 0.01 * (group + 1) * (index + 1) if success else None,
                "total_seconds": 0.02 * (group + 1) * (index + 1) if success else None,
                "peak_memory_mb": 100.0 * (group + 1) if success else None,
                "ciphertext_bytes": 1000 * (group + 1) * (index + 1) if success else None,
                "key_bytes": 2000 * (index + 1) if success else None,
                "max_absolute_error": 1e-8 * (index + 1) if success else None,
                "plan": {
                    "scheme": scheme,
                    "ring_dimension": 8192 << min(group, 1),
                    "window_size": 512,
                    "variants_per_ciphertext": 512,
                    "chunk_count": group + 1,
                    "feature_packing": True,
                    "reduction_tree_strategy": "binary_tree",
                    "chunk_aggregation_strategy": "reduce_then_add",
                    "scale_bits": 40 if scheme == "CKKS" else None,
                    "fixed_point_scaling_factor": None if scheme == "CKKS" else 100000,
                },
            })
    return build_dataset(rows)


class AutoHEPipelineTests(unittest.TestCase):
    def test_storage_roundtrips_single_json_record(self) -> None:
        record = {
            "status": "ok",
            "scheme": "CKKS",
            "plan": {"scheme": "CKKS", "ring_dimension": 8192},
        }
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "measurement.json"
            write_records([record], path)
            restored = read_records(path)
        self.assertEqual(restored, [record])

    def test_parquet_roundtrip_when_optional_dependency_is_available(self) -> None:
        try:
            import pandas  # noqa: F401
            import pyarrow  # noqa: F401
        except ImportError:
            self.skipTest("optional Parquet dependencies are unavailable")
        rows = [
            {"status": "ok", "repeat": 1, "plan": {"scheme": "CKKS"}},
            {"status": "failed", "repeat": "warmup", "plan": {"scheme": "BFV"}},
        ]
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "measurements.parquet"
            write_records(rows, path)
            restored = read_records(path)
        self.assertEqual(len(restored), 2)
        self.assertEqual(restored[0]["status"], "ok")
        self.assertEqual(restored[1]["repeat"], "warmup")

    def test_dataset_storage_model_and_group_split(self) -> None:
        rows = measured_rows()
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            dataset = root / "dataset.jsonl"
            write_records(rows, dataset)
            loaded = read_records(dataset)
            model = train_cost_models(loaded, group_key="workload_group", test_fraction=0.25)
            model.save(root / "model")
            restored = CostModelBundle.load(root / "model")
        self.assertIn("evaluation_seconds", restored.regressors)
        self.assertEqual(restored.split["train_rows"] + restored.split["test_rows"], 11)
        prediction = restored.predict(rows[0])
        self.assertGreaterEqual(prediction["evaluation_seconds"], 0.0)
        self.assertIn("execution_success_probability", prediction)
        self.assertIn("spearman", restored.metrics["evaluation_seconds"])
        self.assertEqual(restored.target_scopes["key_bytes"], "source_defined")

    def test_key_size_model_rejects_mixed_measurement_scopes(self) -> None:
        rows = measured_rows()
        rows[0]["key_bytes_scope"] = "legacy_context+public_key_proxy"
        with self.assertRaisesRegex(ValueError, "incomparable scopes"):
            train_cost_models(
                rows,
                group_key="single_training_group",
                test_fraction=0.25,
            )

    def test_generalization_suite_retains_unavailable_splits(self) -> None:
        rows = measured_rows()
        for index, row in enumerate(rows):
            row["group_variant_count"] = f"variants:{row['matched_variant_count']}"
            row["group_sample_count"] = f"samples:{row['sample_count']}"
            row["group_window_size"] = "window:512"
            row["group_parameter_combination"] = f"p:{index % 3}"
            row["group_scheme_configuration"] = row["scheme"]
            row["group_pgs"] = f"PGS{index % 2}"
            row["group_source_domain"] = "synthetic" if index < 6 else "real"
        evaluations = evaluate_generalization_splits(
            rows, backend="ridge", test_fraction=0.25, seed=9
        )
        self.assertEqual(set(evaluations), {
            "unseen_variant_count",
            "unseen_sample_count",
            "unseen_window_size",
            "unseen_parameter_combination",
            "unseen_scheme_configuration",
            "one_pgs_to_another",
            "synthetic_to_real",
        })
        self.assertEqual(evaluations["synthetic_to_real"].status, "ok")
        self.assertEqual(evaluations["synthetic_to_real"].test_groups, ("real",))
        self.assertEqual(evaluations["unseen_window_size"].status, "unavailable")

    def test_ranker_filters_deterministically(self) -> None:
        rows = measured_rows()
        model = train_cost_models(rows, group_key="workload_group", test_fraction=0.25)
        space = CandidateSpace(
            schemes=("CKKS",),
            ring_dimensions=(8192,),
            window_sizes=(512,),
            scaling_modulus_sizes=(40,),
            aggregation_strategies=("reduce_then_add",),
        )
        candidates = generate_candidates(512, space)
        ranked = rank_candidates(
            candidates,
            model,
            {
                "sample_count": 1,
                "matched_variant_count": 512,
                "genotype_min": 0.0,
                "genotype_max": 2.0,
                "sum_absolute_weights": 1.0,
                "imputed_dosage": False,
            },
            Constraints(max_memory_mb=10_000, max_absolute_error=1.0),
        )
        self.assertEqual(len(ranked), 1)
        self.assertTrue(ranked[0].feasible)

    def test_tuner_pairs_predictions_with_measured_residuals(self) -> None:
        rows = measured_rows()
        model = train_cost_models(rows, group_key="workload_group", test_fraction=0.25)
        candidate = generate_candidates(
            512,
            CandidateSpace(
                schemes=("CKKS",),
                ring_dimensions=(8192,),
                window_sizes=(512,),
                scaling_modulus_sizes=(40,),
                aggregation_strategies=("reduce_then_add",),
            ),
        )[0]

        def evaluator(plan: FHEPlan) -> dict:
            return {
                "status": "ok",
                "measurement_kind": "measured_summary",
                "scheme": plan.scheme,
                "plan": asdict(plan),
                "evaluation_seconds": 0.02,
                "total_seconds": 0.04,
                "peak_memory_mb": 100.0,
                "ciphertext_bytes": 1000,
                "key_bytes": 2000,
                "max_absolute_error": 1e-8,
            }

        result = tune(
            [candidate],
            model,
            {
                "sample_count": 1,
                "matched_variant_count": 512,
                "genotype_min": 0.0,
                "genotype_max": 2.0,
                "sum_absolute_weights": 1.0,
                "imputed_dosage": False,
            },
            Constraints(max_memory_mb=10_000, max_absolute_error=1.0),
            [[0.0] * 512],
            [0.0] * 512,
            top_k=1,
            evaluator=evaluator,
        )
        self.assertEqual(result.selection_status, "measured_and_validated")
        record = result.measured_records[0]
        self.assertEqual(record["search_stage"], "top_k_actual_validation")
        self.assertIn("evaluation_seconds", record["predicted_targets"])
        self.assertIn("evaluation_seconds", record["prediction_residuals"])

    def test_failed_measurement_is_not_paired_with_a_prediction(self) -> None:
        rows = measured_rows()
        model = train_cost_models(rows, group_key="workload_group", test_fraction=0.25)
        candidate = FHEPlan(
            scheme="CKKS",
            ring_dimension=8192,
            slot_count=4096,
            batch_size=512,
            window_size=512,
            variants_per_ciphertext=512,
            rotation_indices=tuple(1 << index for index in range(9)),
            scaling_modulus_size=40,
            first_modulus_size=60,
            scale_bits=40,
        )
        result = tune(
            [candidate],
            model,
            {
                "sample_count": 1,
                "matched_variant_count": 512,
                "genotype_min": 0.0,
                "genotype_max": 2.0,
                "sum_absolute_weights": 1.0,
                "imputed_dosage": False,
            },
            Constraints(max_memory_mb=10_000, max_absolute_error=1.0),
            [[0.0] * 512],
            [0.0] * 512,
            top_k=1,
            evaluator=lambda plan: {
                "status": "failed",
                "scheme": plan.scheme,
                "plan": asdict(plan),
                "total_seconds": 0.1,
                "failure_type": "backend_failure",
            },
        )
        record = result.measured_records[0]
        self.assertEqual(record["prediction_residuals"], {})
        self.assertEqual(
            record["prediction_comparison_status"],
            "unavailable_failed_measurement",
        )
        with tempfile.TemporaryDirectory() as raw:
            table = Path(raw) / "comparison.md"
            write_prediction_comparison_table([record], table)
            content = table.read_text(encoding="utf-8")
        self.assertNotIn("| CKKS |", content)

    def test_openfhe_helper_reduction_and_key_size_accounting(self) -> None:
        class FakeContext:
            def EvalRotate(self, values, amount):
                return values[amount:] + values[:amount]

            def EvalAdd(self, left, right):
                return [a + b for a, b in zip(left, right)]

            def EvalSum(self, values, active):
                return [sum(values[:active]), *([0] * (len(values) - 1))]

            def SerializeEvalMultKey(self, path, _binary):
                Path(path).write_bytes(b"m" * 7)
                return True

            def SerializeEvalAutomorphismKey(self, path, _binary):
                Path(path).write_bytes(b"r" * 11)
                return True

        context = FakeContext()
        self.assertEqual(_reduce_slots(context, [1, 2, 3, 0], 3, "binary_tree")[0], 6)
        self.assertEqual(_reduce_slots(context, [1, 2, 3, 0], 3, "eval_sum")[0], 6)

        class KeyPair:
            publicKey = "public"
            secretKey = "secret"

        lengths = {"public": 3, "secret": 4}
        sizes = _serialize_key_material(
            context,
            KeyPair(),
            lambda value, _binary: b"x" * (
                5 if value is context else lengths[value]
            ),
            object(),
        )
        self.assertEqual(sizes, (5, 3, 4, 7, 11, 25))

    def test_openfhe_absence_is_an_explicit_unavailable_record(self) -> None:
        plan = FHEPlan(
            scheme="CKKS", ring_dimension=8192, slot_count=4096,
            batch_size=4, window_size=4, variants_per_ciphertext=4,
            rotation_indices=(1, 2), scaling_modulus_size=40,
            first_modulus_size=60, scale_bits=40,
        )
        result = run_openfhe_plan(
            plan, [[0.0, 1.0, 2.0, 1.0]], [0.1, 0.2, -0.1, 0.3],
            constraints=Constraints(max_absolute_error=1e-4, max_relative_error=1e-3),
        )
        self.assertEqual(result.measurement_kind, "measured")
        self.assertIn(result.status, {"ok", "unavailable"})
        if result.status == "unavailable":
            self.assertFalse(result.decryption_success)
            self.assertIsNotNone(result.failure_reason)

    def test_pareto_and_bayesian_suggestion(self) -> None:
        rows = [
            {"evaluation_seconds": 1, "peak_memory_mb": 4, "ciphertext_bytes": 3,
             "key_bytes": 2, "max_absolute_error": 0.1},
            {"evaluation_seconds": 2, "peak_memory_mb": 5, "ciphertext_bytes": 4,
             "key_bytes": 3, "max_absolute_error": 0.2},
            {"evaluation_seconds": 3, "peak_memory_mb": 1, "ciphertext_bytes": 2,
             "key_bytes": 1, "max_absolute_error": 0.05},
        ]
        self.assertEqual(len(pareto_frontier(rows)), 2)
        candidates = generate_candidates(
            512,
            CandidateSpace(
                schemes=("CKKS",), ring_dimensions=(8192, 16384),
                window_sizes=(512,), scaling_modulus_sizes=(40,),
                aggregation_strategies=("reduce_then_add",),
            ),
        )
        suggestions = suggest(candidates, [(candidates[0], 1.0)])
        self.assertEqual(len(suggestions), 1)
        self.assertGreaterEqual(suggestions[0].expected_improvement, 0.0)

    def test_search_baselines_ablation_and_synthetic_are_reproducible(self) -> None:
        space = CandidateSpace(
            schemes=("CKKS", "BFV"),
            ring_dimensions=(8192,),
            window_sizes=(512,),
            scaling_modulus_sizes=(40,),
            fixed_point_scales=(100000,),
            aggregation_strategies=("reduce_then_add",),
        )
        candidates = generate_candidates(512, space)
        first = run_baseline_search(
            candidates,
            lambda plan: float(plan.ring_dimension + (0 if plan.scheme == "CKKS" else 1)),
            method="random",
            budget=2,
            seed=7,
        )
        second = run_baseline_search(
            candidates,
            lambda plan: float(plan.ring_dimension + (0 if plan.scheme == "CKKS" else 1)),
            method="random",
            budget=2,
            seed=7,
        )
        self.assertEqual(first.observations, second.observations)
        ablated = apply_plan_ablation(candidates, settings("no_scheme_selection"))
        self.assertTrue(ablated)
        self.assertEqual({plan.scheme for plan in ablated}, {"CKKS"})
        workload_a = make_synthetic_workload(
            samples=3, variants=8, genotype_encoding="imputed",
            missing_rate=0.2, seed=5,
        )
        workload_b = make_synthetic_workload(
            samples=3, variants=8, genotype_encoding="imputed",
            missing_rate=0.2, seed=5,
        )
        self.assertEqual(workload_a, workload_b)
        dosages, weights = apply_missing_policy(workload_a, "zero")
        self.assertEqual(len(dosages), 3)
        self.assertEqual(len(weights), 8)

    def test_emitter_outputs_compilable_script_and_plan(self) -> None:
        plan = FHEPlan(
            scheme="CKKS", ring_dimension=8192, slot_count=4096,
            window_size=512, variants_per_ciphertext=512,
            rotation_indices=tuple(1 << i for i in range(9)),
            scaling_modulus_size=40, first_modulus_size=60, scale_bits=40,
        )
        with tempfile.TemporaryDirectory() as raw:
            config, script = emit(
                plan,
                Path(raw) / "run.py",
                dosage_path=r"data\harmonized\dosage.csv",
                weight_path=r"data\harmonized\weights.csv",
            )
            script_text = script.read_text(encoding="utf-8")
            compile(script_text, str(script), "exec")
            saved = json.loads(config.read_text(encoding="utf-8"))
        self.assertEqual(saved["scheme"], "CKKS")
        self.assertIn("data/harmonized/dosage.csv", script_text)

    def test_cli_help_is_available(self) -> None:
        with self.assertRaises(SystemExit) as context:
            main(["--help"])
        self.assertEqual(context.exception.code, 0)


if __name__ == "__main__":
    unittest.main()
