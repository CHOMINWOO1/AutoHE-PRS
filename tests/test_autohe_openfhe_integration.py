from __future__ import annotations

from dataclasses import asdict
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from autohe_prs.config import Constraints, FHEPlan
from autohe_prs.fhe.emitter import emit
from autohe_prs.fhe.runner import run_openfhe_plan


OPENFHE_AVAILABLE = importlib.util.find_spec("openfhe") is not None


@unittest.skipUnless(
    OPENFHE_AVAILABLE,
    "OpenFHE integration tests run in the pinned OpenFHE container",
)
class OpenFHEIntegrationTests(unittest.TestCase):
    dosages = [
        [0.0, 1.0, 2.0, 1.0, 0.0, 2.0],
        [2.0, 0.0, 1.0, 2.0, 1.0, 0.0],
    ]
    weights = [0.1, -0.2, 0.05, 0.4, -0.3, 0.2]
    constraints = Constraints(
        max_absolute_error=2e-5,
        max_relative_error=2e-4,
        max_memory_mb=None,
    )

    def _plan(
        self,
        scheme: str,
        aggregation: str,
        reduction: str,
    ) -> FHEPlan:
        integer = scheme in {"BFV", "BGV"}
        ring_dimension = 16_384 if scheme == "BGV" else 8_192
        return FHEPlan(
            scheme=scheme,
            ring_dimension=ring_dimension,
            slot_count=ring_dimension if integer else ring_dimension // 2,
            batch_size=4,
            window_size=4,
            variants_per_ciphertext=4,
            chunk_count=2,
            rotation_indices=(1, 2),
            chunk_aggregation_strategy=aggregation,
            reduction_tree_strategy=reduction,
            scaling_modulus_size=None if integer else 40,
            first_modulus_size=None if integer else 60,
            scale_bits=None if integer else 40,
            scaling_technique=None if integer else "FLEXIBLEAUTO",
            plaintext_modulus=2_147_352_577 if integer else None,
            fixed_point_scaling_factor=1_000_000 if integer else None,
        )

    def test_ckks_bfv_bgv_end_to_end_with_multiple_chunks(self) -> None:
        cases = (
            ("CKKS", "reduce_then_add", "binary_tree"),
            ("BFV", "add_then_reduce", "eval_sum"),
            ("BGV", "reduce_then_add", "binary_tree"),
        )
        for scheme, aggregation, reduction in cases:
            with self.subTest(scheme=scheme):
                result = run_openfhe_plan(
                    self._plan(scheme, aggregation, reduction),
                    self.dosages,
                    self.weights,
                    constraints=self.constraints,
                )
                self.assertEqual(result.status, "ok", result.failure_reason)
                self.assertTrue(result.decryption_success)
                self.assertTrue(result.error_constraint_feasible)
                self.assertEqual(result.key_bytes_scope, "public+secret+multiplication+rotation")
                self.assertGreater(result.key_bytes or 0, 0)
                self.assertGreater(result.rotation_key_bytes or 0, 0)
                self.assertLessEqual(
                    result.max_absolute_error or float("inf"),
                    self.constraints.max_absolute_error or 0.0,
                )

    def test_generated_script_executes_without_manual_edits(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            dosage = root / "dosage.csv"
            weight = root / "weights.csv"
            dosage.write_text(
                "sample_id,v1,v2,v3,v4,v5,v6\n"
                "s1,0,1,2,1,0,2\n"
                "s2,2,0,1,2,1,0\n",
                encoding="utf-8",
            )
            weight.write_text(
                "effect_weight\n0.1\n-0.2\n0.05\n0.4\n-0.3\n0.2\n",
                encoding="utf-8",
            )
            _, script = emit(
                self._plan("CKKS", "reduce_then_add", "binary_tree"),
                root / "run.py",
                constraints=asdict(self.constraints),
                dosage_path=str(dosage),
                weight_path=str(weight),
            )
            environment = dict(os.environ)
            environment["PYTHONPATH"] = str(
                Path(__file__).resolve().parents[1] / "src"
            )
            process = subprocess.run(
                [sys.executable, str(script)],
                capture_output=True,
                text=True,
                timeout=180,
                env=environment,
                check=False,
            )
            self.assertEqual(process.returncode, 0, process.stderr or process.stdout)
            payload = json.loads(process.stdout)
            self.assertEqual(payload["status"], "ok")
            self.assertEqual(payload["measurement_kind"], "measured")
