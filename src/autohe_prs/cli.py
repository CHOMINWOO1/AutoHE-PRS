"""`autohe-prs` command-line interface."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import math
from pathlib import Path
import sys
from time import perf_counter
from typing import Any

from .benchmark.dataset import build_dataset
from .benchmark.runner import run_isolated_repeated
from .benchmark.storage import read_records, write_records
from .config import Constraints, FHEPlan, dump_json, load_experiment_spec, load_mapping
from .fhe.emitter import emit
from .fhe.validator import validate_plan
from .genomics.harmonize import harmonize, write_harmonized
from .genomics.dosage import read_dosage_matrix
from .genomics.io import read_harmonized
from .genomics.pgs_parser import read_pgs
from .genomics.plaintext import prs
from .genomics.vcf_parser import read_vcf
from .models.auto import load_cost_model, train_auto_cost_models
from .models.evaluate import evaluate_generalization_splits
from .optimize.candidate import CandidateSpace, generate_candidates
from .optimize.baselines import run_baseline_search
from .optimize.tuner import tune
from .reporting.summary import write_experiment_summary
from .reporting.figures import plot_latency_memory
from .reporting.tables import write_prediction_comparison_table, write_scheme_table


def _json_print(value: Any) -> None:
    print(json.dumps(value, indent=2, sort_keys=True))


def _workload_features(dosages: list[list[float]], weights: list[float]) -> dict[str, Any]:
    values = [value for row in dosages for value in row]
    weight_mean = math.fsum(weights) / len(weights) if weights else 0.0
    weight_variance = (
        math.fsum((weight - weight_mean) ** 2 for weight in weights) / len(weights)
        if weights else 0.0
    )
    decimal_precision = max(
        (
            len(format(abs(weight), ".15f").rstrip("0").partition(".")[2])
            for weight in weights
        ),
        default=0,
    )
    return {
        "sample_count": len(dosages),
        "matched_variant_count": len(weights),
        "score_count": 1,
        "genotype_min": min(values, default=0.0),
        "genotype_max": max(values, default=0.0),
        "imputed_dosage": any(not math.isclose(value, round(value)) for value in values),
        "genotype_encoding_type": (
            "imputed" if any(not math.isclose(value, round(value)) for value in values)
            else "hard_call"
        ),
        "missing_rate": 0.0,
        "max_absolute_weight": max((abs(weight) for weight in weights), default=0.0),
        "sum_absolute_weights": math.fsum(abs(weight) for weight in weights),
        "weight_mean": weight_mean,
        "weight_standard_deviation": math.sqrt(weight_variance),
        "weight_decimal_precision": decimal_precision,
        "score_sparsity": (
            sum(math.isclose(weight, 0.0, abs_tol=0.0) for weight in weights) / len(weights)
            if weights else 0.0
        ),
        "score_overlap_ratio": 1.0,
        "expected_score_bound": max((abs(value) for value in values), default=0.0)
        * math.fsum(abs(weight) for weight in weights),
    }


def _space(path: str | None, allowed: tuple[str, ...]) -> CandidateSpace:
    if not path:
        return CandidateSpace(schemes=allowed)
    raw = load_mapping(path)
    raw = raw.get("search_space", raw)
    names = CandidateSpace.__dataclass_fields__
    values: dict[str, Any] = {}
    for key, value in raw.items():
        if key in names:
            values[key] = tuple(value) if isinstance(value, list) else value
    values["schemes"] = tuple(scheme for scheme in values.get("schemes", allowed)
                              if scheme.upper() in {item.upper() for item in allowed})
    if not values["schemes"]:
        raise ValueError("search space and allowed_schemes have no scheme in common")
    return CandidateSpace(**values)


def _load_plan_envelope(path: str | Path) -> tuple[FHEPlan, dict[str, Any]]:
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    value = raw.get("selected_plan", raw.get("plan", raw))
    if value is None:
        raise ValueError("plan file contains no selected plan")
    return FHEPlan(**value), raw


def _objective_value(record: dict[str, Any], objective: str) -> float | None:
    keys = {
        "evaluation_latency": "evaluation_seconds",
        "total_latency": "total_seconds",
        "memory": "peak_memory_mb",
        "ciphertext_size": "ciphertext_bytes",
        "key_size": "key_bytes",
    }
    if objective == "balanced_multi_objective":
        targets = (
            "evaluation_seconds",
            "peak_memory_mb",
            "ciphertext_bytes",
            "key_bytes",
            "max_absolute_error",
        )
        if any(record.get(key) is None for key in targets):
            return None
        return math.fsum(
            math.log1p(max(0.0, float(record[key]))) for key in targets
        )
    key = keys.get(objective)
    if key is None:
        raise ValueError(f"unsupported objective: {objective}")
    return float(record[key]) if record.get(key) is not None else None


def command_harmonize(args: argparse.Namespace) -> int:
    total_started = perf_counter()
    pgs_started = perf_counter()
    pgs = read_pgs(args.pgs, chrom=args.chrom)
    pgs_parse_seconds = perf_counter() - pgs_started
    genotype_started = perf_counter()
    if args.vcf:
        vcf = read_vcf(
            args.vcf,
            prefer_ds=not args.prefer_gt,
            genome_build=args.genome_build,
            sample_limit=args.max_samples,
            variant_keys={variant.key for variant in pgs.variants},
        )
    else:
        vcf = read_dosage_matrix(args.dosage_matrix, genome_build=args.genome_build)
    genotype_parse_seconds = perf_counter() - genotype_started
    harmonization_started = perf_counter()
    result = harmonize(
        vcf, pgs, missing_policy=args.missing_policy,
        allow_strand_ambiguous=args.allow_strand_ambiguous,
    )
    harmonization_seconds = perf_counter() - harmonization_started
    write_started = perf_counter()
    write_harmonized(result, args.output, output_format=args.format)
    output_write_seconds = perf_counter() - write_started
    metrics = {
        "preprocessing_measurement_kind": "measured",
        "pgs_parse_seconds": pgs_parse_seconds,
        "genotype_parse_seconds": genotype_parse_seconds,
        "harmonization_seconds": harmonization_seconds,
        "output_write_seconds": output_write_seconds,
        "total_preprocessing_seconds": perf_counter() - total_started,
    }
    dump_json(metrics, Path(args.output) / "preprocessing_metrics.json")
    _json_print({"report": asdict(result.report), "timing": metrics})
    return 0


def command_plaintext(args: argparse.Namespace) -> int:
    sample_ids, dosages, weights = read_harmonized(args.dosage, args.weights)
    scores = prs(dosages, weights)
    result = {
        "sample_count": len(sample_ids),
        "variant_count": len(weights),
        "scores": dict(zip(sample_ids, scores)),
        "measurement_kind": "plaintext_reference",
    }
    if args.output:
        dump_json(result, args.output)
    _json_print(result)
    return 0


def command_benchmark(args: argparse.Namespace) -> int:
    spec = load_experiment_spec(args.spec)
    preprocessing_started = perf_counter()
    _, dosages, weights = read_harmonized(
        spec.workload.genotype_path, spec.workload.weight_path
    )
    preprocessing_seconds = perf_counter() - preprocessing_started
    preprocessing_record: dict[str, Any] = {}
    preprocessing_path = (
        Path(spec.workload.genotype_path).resolve().parent
        / "preprocessing_metrics.json"
    )
    if preprocessing_path.exists():
        preprocessing_record = json.loads(
            preprocessing_path.read_text(encoding="utf-8")
        )
    features = _workload_features(dosages, weights)
    candidates = generate_candidates(len(weights), _space(args.search_space, spec.allowed_schemes))
    if args.max_candidates:
        candidates = candidates[: args.max_candidates]
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    for candidate_id, plan in enumerate(candidates):
        validation = validate_plan(
            plan, samples=len(dosages), variants=len(weights),
            genotype_min=features["genotype_min"], genotype_max=features["genotype_max"],
            sum_abs_weights=features["sum_absolute_weights"], constraints=spec.constraints,
            imputed_dosage=features["imputed_dosage"],
        )
        if not validation.valid:
            records.append({
                "candidate_id": candidate_id,
                "scheme": plan.scheme,
                "backend": "openfhe-python",
                "measurement_kind": "validation_only",
                "status": "rejected",
                "failure_type": "constraint_validation",
                "failure_reason": "; ".join(issue.failure_type for issue in validation.issues),
                "samples": len(dosages),
                "variants": len(weights),
                "plan": asdict(plan),
                "validation": asdict(validation),
                "harmonized_input_loading_seconds": preprocessing_seconds,
                **preprocessing_record,
                **features,
            })
            continue
        summary = run_isolated_repeated(
            plan, dosages, weights, constraints=spec.constraints,
            repeats=args.repeats, warmup=not args.no_warmup,
            timeout_seconds=args.timeout,
        )
        for record in summary.records:
            record.update({
                "candidate_id": candidate_id,
                "harmonized_input_loading_seconds": preprocessing_seconds,
                **preprocessing_record,
                **features,
            })
            records.append(record)
        summaries.append({
            "candidate_id": candidate_id,
            "scheme": plan.scheme,
            "plan": asdict(plan),
            "warmup_status": summary.warmup_status,
            "warmup_record": summary.warmup_record,
            "statistics": {
                name: asdict(values) for name, values in summary.statistics.items()
            },
        })
    write_records(records, output / "benchmark.jsonl")
    write_records(records, output / "benchmark.csv")
    storage_manifest: dict[str, Any] = {
        "jsonl": {"status": "ok", "path": str(output / "benchmark.jsonl")},
        "csv": {"status": "ok", "path": str(output / "benchmark.csv")},
    }
    try:
        write_records(records, output / "benchmark.parquet")
    except RuntimeError as exc:
        storage_manifest["parquet"] = {"status": "unavailable", "reason": str(exc)}
    else:
        storage_manifest["parquet"] = {
            "status": "ok",
            "path": str(output / "benchmark.parquet"),
        }
    dump_json(storage_manifest, output / "storage_manifest.json")
    dump_json(summaries, output / "benchmark_summary.json")
    write_experiment_summary(records, output)
    _json_print({
        "candidates": len(candidates),
        "records": len(records),
        "output": str(output),
        "storage": storage_manifest,
    })
    return 0 if any(row.get("status") == "ok" for row in records) else 2


def command_build_dataset(args: argparse.Namespace) -> int:
    sources = args.results if isinstance(args.results, list) else [args.results]
    paths: list[Path] = []
    for value in sources:
        source = Path(value)
        paths.extend([source] if source.is_file() else sorted(source.rglob("*.jsonl")))
    records = [row for path in paths for row in read_records(path)]
    dataset = build_dataset(records)
    write_records(dataset, args.output)
    _json_print({"input_records": len(records), "dataset_rows": len(dataset), "output": args.output})
    return 0


def command_train(args: argparse.Namespace) -> int:
    rows = read_records(args.dataset)
    model = train_auto_cost_models(
        rows, backend=args.backend, group_key=args.group_key,
        test_fraction=args.test_fraction, seed=args.seed
    )
    model.save(args.output)
    result = {
        "models": sorted(model.regressors),
        "metrics": model.metrics,
        "split": model.split,
        "output": args.output,
    }
    if args.evaluate_splits:
        split_results = evaluate_generalization_splits(
            rows,
            backend=args.backend,
            test_fraction=args.test_fraction,
            seed=args.seed,
        )
        split_path = Path(args.output) / "generalization_evaluation.json"
        dump_json(
            {name: asdict(value) for name, value in split_results.items()},
            split_path,
        )
        result["generalization_evaluation"] = str(split_path)
    _json_print(result)
    return 0


def command_tune(args: argparse.Namespace) -> int:
    spec = load_experiment_spec(args.spec)
    _, dosages, weights = read_harmonized(
        spec.workload.genotype_path, spec.workload.weight_path
    )
    features = _workload_features(dosages, weights)
    model = load_cost_model(args.model)
    candidates = generate_candidates(len(weights), _space(args.search_space, spec.allowed_schemes))
    raw_validation_records: list[dict[str, Any]] = []
    validation_summaries: list[dict[str, Any]] = []

    def isolated_evaluator(plan: FHEPlan) -> dict[str, Any]:
        summary = run_isolated_repeated(
            plan,
            dosages,
            weights,
            constraints=spec.constraints,
            repeats=args.validation_repeats,
            warmup=not args.no_validation_warmup,
            timeout_seconds=args.validation_timeout,
        )
        if summary.warmup_record is not None:
            raw_validation_records.append(dict(summary.warmup_record))
        raw_validation_records.extend(dict(record) for record in summary.records)
        successful = [
            record for record in summary.records
            if record.get("status") == "ok"
        ]
        validation_summaries.append({
            "plan": asdict(plan),
            "warmup_status": summary.warmup_status,
            "repeat_count": len(summary.records),
            "successful_repeats": len(successful),
            "statistics": {
                name: asdict(statistic)
                for name, statistic in summary.statistics.items()
            },
        })
        if len(successful) != args.validation_repeats:
            failed = next(
                (record for record in summary.records if record.get("status") != "ok"),
                summary.records[0] if summary.records else {},
            )
            return {
                **failed,
                "scheme": plan.scheme,
                "plan": asdict(plan),
                "measurement_kind": "measured_summary",
                "status": "failed",
                "failure_type": "repeated_validation_failure",
                "failure_reason": (
                    f"{len(successful)}/{args.validation_repeats} measured repeats succeeded; "
                    f"first failure: {failed.get('failure_type')}: "
                    f"{failed.get('failure_reason')}"
                ),
            }
        aggregate = dict(successful[0])
        for name, statistic in summary.statistics.items():
            value: float | int = statistic.median
            if name.endswith("_bytes"):
                value = int(round(value))
            aggregate[name] = value
        aggregate["measurement_kind"] = "measured_summary"
        aggregate["repeat_count"] = args.validation_repeats
        aggregate["successful_repeats"] = len(successful)
        return aggregate

    result = tune(
        candidates, model, features, spec.constraints, dosages, weights,
        objective=args.objective or spec.objective,
        top_k=args.top_k or spec.top_k_validation,
        use_bayesian_refinement=(
            spec.use_bayesian_refinement and not args.no_bayesian_refinement
        ),
        bayesian_iterations=args.bayesian_iterations,
        evaluator=isolated_evaluator,
    )
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    envelope = {
        "selection_status": result.selection_status,
        "selected_plan": asdict(result.selected_plan) if result.selected_plan else None,
        "workload_features": features,
        "constraints": asdict(spec.constraints),
        "objective": args.objective or spec.objective,
        "selected_measurement": result.selected_record,
    }
    dump_json(envelope, output / "best_plan.json")
    write_records([asdict(item) for item in result.ranking], output / "ranking.jsonl")
    write_records(result.measured_records, output / "top_k_measured.jsonl")
    write_records(raw_validation_records, output / "top_k_raw_measured.jsonl")
    dump_json(validation_summaries, output / "top_k_benchmark_summary.json")
    write_records(result.pareto_records, output / "pareto.jsonl")
    if args.run_baselines:
        objective = args.objective or spec.objective
        deterministic_candidates = []
        for plan in candidates:
            validation = validate_plan(
                plan,
                samples=len(dosages),
                variants=len(weights),
                genotype_min=features["genotype_min"],
                genotype_max=features["genotype_max"],
                sum_abs_weights=features["sum_absolute_weights"],
                constraints=spec.constraints,
                imputed_dosage=features["imputed_dosage"],
            )
            if validation.valid:
                deterministic_candidates.append(plan)
        baseline_records: list[dict[str, Any]] = []
        baseline_raw_records: list[dict[str, Any]] = []
        baseline_results: list[dict[str, Any]] = []

        def run_method(
            method: str,
            *,
            fixed_plan: FHEPlan | None = None,
            result_name: str | None = None,
        ) -> None:
            label = result_name or method

            def evaluate(plan: FHEPlan) -> float | None:
                before = len(raw_validation_records)
                record = isolated_evaluator(plan)
                for raw_record in raw_validation_records[before:]:
                    raw_record["search_stage"] = f"baseline:{label}"
                    baseline_raw_records.append(dict(raw_record))
                record["search_stage"] = f"baseline:{label}"
                baseline_records.append(record)
                return (
                    _objective_value(record, objective)
                    if record.get("status") == "ok"
                    else None
                )

            result_value = run_baseline_search(
                deterministic_candidates,
                evaluate,
                method=method,
                budget=1 if method in {"fixed", "minimal"} else args.baseline_budget,
                seed=args.seed,
                fixed_plan=fixed_plan,
            )
            payload = asdict(result_value)
            payload["method"] = label
            baseline_results.append(payload)

        for scheme in spec.allowed_schemes:
            manual = next(
                (
                    plan for plan in deterministic_candidates
                    if plan.scheme.upper() == scheme.upper()
                    and plan.ring_dimension == 16_384
                    and plan.window_size == 4_096
                    and plan.chunk_aggregation_strategy == "reduce_then_add"
                    and plan.reduction_tree_strategy == "binary_tree"
                    and (
                        plan.scheme.upper() != "CKKS"
                        or plan.scaling_modulus_size == 50
                    )
                    and (
                        plan.scheme.upper() == "CKKS"
                        or plan.fixed_point_scaling_factor == 10_000_000
                    )
                ),
                None,
            )
            if manual is not None:
                run_method("fixed", fixed_plan=manual, result_name=f"manual_{scheme.upper()}")
        for method in ("minimal", "grid", "random", "bayesian"):
            if deterministic_candidates:
                run_method(method)
        successful_values = [
            float(item["best_value"])
            for item in baseline_results
            if item.get("best_value") is not None
        ]
        proposed_best = (
            _objective_value(result.selected_record, objective)
            if result.selected_record is not None
            else None
        )
        comparison_values = successful_values + (
            [proposed_best] if proposed_best is not None else []
        )
        global_best = min(comparison_values) if comparison_values else None
        for item in baseline_results:
            item["regret"] = (
                float(item["best_value"]) - global_best
                if item.get("best_value") is not None and global_best is not None
                else None
            )
        write_records(baseline_records, output / "baseline_measured.jsonl")
        write_records(
            baseline_raw_records,
            output / "baseline_raw_measured.jsonl",
        )
        dump_json(
            {
                "objective": objective,
                "budget_per_search": args.baseline_budget,
                "best_measured_baseline": (
                    min(successful_values) if successful_values else None
                ),
                "best_measured_overall": global_best,
                "proposed": {
                    "best_value": proposed_best,
                    "measured_candidates": len(result.measured_records),
                    "regret": (
                        proposed_best - global_best
                        if proposed_best is not None and global_best is not None
                        else None
                    ),
                },
                "results": baseline_results,
            },
            output / "baseline_comparison.json",
        )
    _json_print({**envelope, "output": str(output)})
    return 0 if result.selection_status == "measured_and_validated" else 2


def command_validate(args: argparse.Namespace) -> int:
    plan, envelope = _load_plan_envelope(args.plan)
    features = envelope.get("workload_features", plan.metadata.get("workload_features", {}))
    constraints = Constraints(**envelope.get("constraints", {}))
    required = {"sample_count", "matched_variant_count", "sum_absolute_weights"}
    if not required <= features.keys():
        raise ValueError(f"plan lacks validation workload features: {sorted(required - features.keys())}")
    result = validate_plan(
        plan, samples=int(features["sample_count"]),
        variants=int(features["matched_variant_count"]),
        genotype_min=float(features.get("genotype_min", 0.0)),
        genotype_max=float(features.get("genotype_max", 2.0)),
        sum_abs_weights=float(features["sum_absolute_weights"]),
        constraints=constraints, imputed_dosage=bool(features.get("imputed_dosage", False)),
    )
    _json_print(asdict(result))
    return 0 if result.valid else 2


def command_emit(args: argparse.Namespace) -> int:
    plan, envelope = _load_plan_envelope(args.plan)
    config_path, script_path = emit(
        plan, args.output, constraints=envelope.get("constraints", {}),
        dosage_path=args.dosage, weight_path=args.weights,
    )
    _json_print({"configuration": str(config_path), "script": str(script_path)})
    return 0


def command_report(args: argparse.Namespace) -> int:
    source = Path(args.experiment)
    paths = [source] if source.is_file() else sorted(source.rglob("*.jsonl"))
    records = [row for path in paths for row in read_records(path)]
    json_path, markdown_path = write_experiment_summary(records, args.output)
    output = Path(args.output)
    scheme_table = write_scheme_table(records, output / "scheme_table.md")
    comparison_table = write_prediction_comparison_table(
        records, output / "prediction_vs_measurement.md"
    )
    figure_status: dict[str, Any]
    try:
        figure = plot_latency_memory(records, output / "latency_memory.png")
    except (ValueError, RuntimeError) as exc:
        figure_status = {"status": "unavailable", "reason": str(exc)}
    else:
        figure_status = {"status": "ok", "path": str(figure)}
    manifest = {
        "records": len(records),
        "json": str(json_path),
        "markdown": str(markdown_path),
        "scheme_table": str(scheme_table),
        "prediction_comparison_table": str(comparison_table),
        "latency_memory_figure": figure_status,
    }
    dump_json(manifest, output / "report_manifest.json")
    _json_print(manifest)
    return 0


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(prog="autohe-prs")
    root.add_argument("--version", action="version", version="autohe-prs 0.1.0")
    commands = root.add_subparsers(dest="command", required=True)

    item = commands.add_parser("harmonize", help="harmonize a VCF with a PGS scoring file")
    source = item.add_mutually_exclusive_group(required=True)
    source.add_argument("--vcf")
    source.add_argument("--dosage-matrix")
    item.add_argument("--pgs", required=True)
    item.add_argument("--output", required=True)
    item.add_argument("--genome-build")
    item.add_argument("--chrom", help="stream only one chromosome from large PGS files")
    item.add_argument("--max-samples", type=int)
    item.add_argument("--missing-policy", choices=("error", "exclude_variant", "zero", "mean_dosage"),
                      default="error")
    item.add_argument("--prefer-gt", action="store_true")
    item.add_argument(
        "--allow-strand-ambiguous",
        action="store_true",
        help=(
            "retain A/T and C/G loci only after the caller has validated them "
            "against a trusted build/frequency reference"
        ),
    )
    item.add_argument("--format", choices=("csv", "parquet", "both"), default="csv")
    item.set_defaults(handler=command_harmonize)

    item = commands.add_parser("plaintext", help="calculate plaintext reference PRS")
    item.add_argument("--dosage", required=True)
    item.add_argument("--weights", required=True)
    item.add_argument("--output")
    item.set_defaults(handler=command_plaintext)

    item = commands.add_parser("benchmark", help="run candidate plans in isolated processes")
    item.add_argument("--spec", required=True)
    item.add_argument("--search-space")
    item.add_argument("--output", required=True)
    item.add_argument("--repeats", type=int, default=3)
    item.add_argument("--timeout", type=float, default=300.0)
    item.add_argument("--no-warmup", action="store_true")
    item.add_argument("--max-candidates", type=int, default=0)
    item.set_defaults(handler=command_benchmark)

    item = commands.add_parser("build-dataset", help="build canonical benchmark dataset")
    item.add_argument(
        "--results",
        required=True,
        action="append",
        help="input file/directory; repeat to merge measured and failure records",
    )
    item.add_argument("--output", required=True)
    item.set_defaults(handler=command_build_dataset)

    item = commands.add_parser("train", help="train grouped cost models")
    item.add_argument("--dataset", required=True)
    item.add_argument("--output", required=True)
    item.add_argument("--group-key", default="workload_group")
    item.add_argument("--test-fraction", type=float, default=0.2)
    item.add_argument("--seed", type=int, default=20260728)
    item.add_argument(
        "--backend",
        choices=("auto", "xgboost", "lightgbm", "hist_gradient_boosting",
                 "random_forest", "mlp", "ridge"),
        default="auto",
    )
    item.add_argument(
        "--evaluate-splits",
        action="store_true",
        help="evaluate all mandated group-held-out and synthetic-to-real splits",
    )
    item.set_defaults(handler=command_train)

    item = commands.add_parser("tune", help="rank, validate, and select a plan")
    item.add_argument("--spec", required=True)
    item.add_argument("--model", required=True)
    item.add_argument("--search-space")
    item.add_argument("--objective")
    item.add_argument("--top-k", type=int)
    item.add_argument("--output", default="outputs")
    item.add_argument("--no-bayesian-refinement", action="store_true")
    item.add_argument("--bayesian-iterations", type=int, default=5)
    item.add_argument("--validation-repeats", type=int, default=3)
    item.add_argument("--validation-timeout", type=float, default=300.0)
    item.add_argument("--no-validation-warmup", action="store_true")
    item.add_argument(
        "--run-baselines",
        action="store_true",
        help="measure manual, minimal, grid, random, and uninformed Bayesian baselines",
    )
    item.add_argument("--baseline-budget", type=int, default=10)
    item.add_argument("--seed", type=int, default=20260728)
    item.set_defaults(handler=command_tune)

    item = commands.add_parser("validate", help="deterministically validate a plan")
    item.add_argument("--plan", required=True)
    item.set_defaults(handler=command_validate)

    item = commands.add_parser("emit", help="emit executable OpenFHE-Python code")
    item.add_argument("--plan", required=True)
    item.add_argument("--output", required=True)
    item.add_argument("--dosage", default="data/harmonized/dosage.csv")
    item.add_argument("--weights", default="data/harmonized/weights.csv")
    item.set_defaults(handler=command_emit)

    item = commands.add_parser("report", help="summarize an experiment directory")
    item.add_argument("--experiment", required=True)
    item.add_argument("--output", required=True)
    item.set_defaults(handler=command_report)
    return root


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        return int(args.handler(args))
    except (ValueError, RuntimeError, OSError) as exc:
        print(f"autohe-prs: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
