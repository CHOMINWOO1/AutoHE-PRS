from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "configs" / "regression" / "pgs004941_chr22.json"


@unittest.skipUnless(
    (ROOT / "data" / "autohe_smoke" / "pgs004941_chr22" / "dosage.csv").exists(),
    "real-data regression artifacts are not present",
)
class PGS004941RegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))

    def test_harmonization_counts_and_plaintext_prs_are_stable(self) -> None:
        report = json.loads(
            (
                ROOT
                / "data"
                / "autohe_smoke"
                / "pgs004941_chr22"
                / "harmonization_report.json"
            ).read_text(encoding="utf-8")
        )
        for name, expected in self.manifest["harmonization"].items():
            self.assertEqual(report[name], expected, name)
        plaintext = json.loads(
            (ROOT / "results" / "autohe_smoke" / "plaintext.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(plaintext["variant_count"], 51_225)
        for sample, expected in self.manifest["plaintext_scores"].items():
            self.assertTrue(
                math.isclose(
                    plaintext["scores"][sample],
                    expected,
                    rel_tol=0.0,
                    abs_tol=1e-15,
                ),
                sample,
            )

    def test_harmonized_artifact_hashes_and_chunk_count_are_stable(self) -> None:
        for name in ("dosage", "weights"):
            relative = self.manifest["outputs"][f"{name}_path"]
            expected = self.manifest["outputs"][f"{name}_sha256"]
            digest = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
            self.assertEqual(digest, expected)
        variants = self.manifest["harmonization"]["final_scored_variants"]
        self.assertEqual(
            math.ceil(variants / 4096),
            self.manifest["window_4096_chunk_count"],
        )
