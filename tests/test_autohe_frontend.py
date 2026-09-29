from __future__ import annotations

import gzip
from pathlib import Path
import tempfile
import unittest

from autohe_prs.config import Constraints, FHEPlan, load_experiment_spec
from autohe_prs.fhe.packing import (
    chunk_ranges,
    padding_mask,
    reduction_rotation_indices,
)
from autohe_prs.fhe.reduction import binary_tree_reduce
from autohe_prs.fhe.validator import validate_plan
from autohe_prs.genomics.harmonize import harmonize
from autohe_prs.genomics.dosage import read_dosage_matrix
from autohe_prs.genomics.pgs_parser import read_pgs
from autohe_prs.genomics.plaintext import compare, prs
from autohe_prs.genomics.vcf_parser import gt_to_alt_dosage, read_vcf
from autohe_prs.ir.graph import build_prs_graph
from autohe_prs.ir.nodes import NodeKind
from autohe_prs.optimize.candidate import CandidateSpace, generate_candidates


class AutoHEFrontendTests(unittest.TestCase):
    def _inputs(self, directory: Path) -> tuple[Path, Path]:
        vcf = directory / "tiny.vcf.gz"
        pgs = directory / "PGS999999.txt"
        with gzip.open(vcf, "wt", encoding="utf-8") as handle:
            handle.write(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
                "22\t100\trs1\tA\tG\t.\tPASS\t.\tGT:DS\t0/0:0.0\t0/1:0.8\n"
                "22\t101\trs2\tC\tT\t.\tPASS\t.\tGT\t1/1\t0/1\n"
                "22\t102\trs3\tG\tA\t.\tPASS\t.\tGT\t./.\t0/0\n"
            )
        pgs.write_text(
            "#pgs_id=PGS999999\n"
            "#HmPOS_build=GRCh37\n"
            "chr_name\tchr_position\trsID\teffect_allele\tother_allele\teffect_weight\n"
            "22\t100\trs1\tG\tA\t0.1\n"
            "22\t101\trs2\tC\tT\t-0.2\n"
            "22\t102\trs3\tA\tG\t0.3\n"
            "22\t999\trs_missing\tA\tC\t0.4\n",
            encoding="utf-8",
        )
        return vcf, pgs

    def test_gt_and_ds_parsing_and_harmonization(self) -> None:
        self.assertEqual(gt_to_alt_dosage("0|1"), 1.0)
        self.assertIsNone(gt_to_alt_dosage("./."))
        with tempfile.TemporaryDirectory() as raw:
            vcf_path, pgs_path = self._inputs(Path(raw))
            vcf = read_vcf(vcf_path, genome_build="GRCh37")
            scoring = read_pgs(pgs_path)
            result = harmonize(vcf, scoring, missing_policy="zero")
        self.assertEqual(result.sample_ids, ("S1", "S2"))
        self.assertEqual(result.report.final_scored_variants, 3)
        self.assertEqual(result.report.alt_oriented_variants, 2)
        self.assertEqual(result.report.ref_oriented_variants, 1)
        self.assertEqual(result.report.missing_variants, 1)
        self.assertEqual(result.dosage_matrix[0], [0.0, 0.0, 0.0])
        self.assertEqual(result.dosage_matrix[1], [0.8, 1.0, 0.0])
        self.assertEqual(prs(result.dosage_matrix, result.weights), [0.0, -0.12])

    def test_missing_policy_is_not_silent(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            vcf_path, pgs_path = self._inputs(Path(raw))
            vcf = read_vcf(vcf_path, genome_build="GRCh37")
            scoring = read_pgs(pgs_path)
            with self.assertRaisesRegex(ValueError, "missing genotype"):
                harmonize(vcf, scoring, missing_policy="error")

    def test_preprocessed_dosage_matrix_preserves_allele_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            matrix = root / "dosage.csv"
            pgs_path = root / "score.txt"
            matrix.write_text(
                "sample_id,22:100:A:G:rs1,22:101:C:T:rs2\n"
                "S1,0,2\nS2,1.5,1\n",
                encoding="utf-8",
            )
            pgs_path.write_text(
                "#pgs_id=PGS1\n#genome_build=GRCh37\n"
                "chr_name\tchr_position\trsID\teffect_allele\tother_allele\teffect_weight\n"
                "22\t100\trs1\tG\tA\t0.2\n"
                "22\t101\trs2\tC\tT\t-0.1\n",
                encoding="utf-8",
            )
            dosage = read_dosage_matrix(matrix, genome_build="GRCh37")
            result = harmonize(dosage, read_pgs(pgs_path), missing_policy="error")
        self.assertEqual(result.report.final_scored_variants, 2)
        self.assertEqual(result.dosage_matrix, [[0.0, 0.0], [1.5, 1.0]])

    def test_unambiguous_strand_flip_and_ambiguous_exclusion(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            matrix = root / "dosage.csv"
            pgs_path = root / "score.txt"
            matrix.write_text(
                "sample_id,22:100:A:C,22:101:A:T\nS1,1,2\n",
                encoding="utf-8",
            )
            pgs_path.write_text(
                "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n"
                "22\t100\tG\tT\t0.2\n"
                "22\t101\tT\tA\t0.1\n",
                encoding="utf-8",
            )
            result = harmonize(
                read_dosage_matrix(matrix), read_pgs(pgs_path),
                missing_policy="error",
            )
        self.assertEqual(result.report.final_scored_variants, 1)
        self.assertEqual(result.report.strand_flipped_variants, 1)
        self.assertEqual(result.report.strand_ambiguous_variants, 1)
        self.assertTrue(result.variants[0].strand_flipped)

    def test_packing_reduction_and_ir(self) -> None:
        self.assertEqual(chunk_ranges(9, 4), ((0, 4), (4, 8), (8, 9)))
        self.assertEqual(padding_mask(2, 4), (1, 1, 0, 0))
        self.assertEqual(reduction_rotation_indices(5), (1, 2, 4))
        self.assertEqual(binary_tree_reduce([1, 2, 3, 4, 5], 5)[0], 15)
        graph = build_prs_graph(2, 9, 4)
        graph.validate()
        self.assertEqual(graph.nodes[-1].shape, (2, 1))
        self.assertIn(NodeKind.WEIGHT_INPUT, {node.kind for node in graph.nodes})
        self.assertIn(NodeKind.ENCODE_PLAIN, {node.kind for node in graph.nodes})

    def test_validator_rejects_overflow_and_imputed_fixed_point(self) -> None:
        plan = FHEPlan(
            scheme="BFV",
            ring_dimension=8192,
            slot_count=8192,
            batch_size=4,
            window_size=4,
            variants_per_ciphertext=4,
            rotation_indices=(1, 2),
            plaintext_modulus=257,
            fixed_point_scaling_factor=100,
        )
        result = validate_plan(
            plan,
            samples=1,
            variants=4,
            genotype_min=0.0,
            genotype_max=2.0,
            sum_abs_weights=2.0,
            imputed_dosage=True,
        )
        self.assertFalse(result.valid)
        self.assertIn("fixed_point_wraparound", {issue.failure_type for issue in result.issues})
        self.assertIn("unsupported_imputed_fixed_point", {issue.failure_type for issue in result.issues})

        bgv = FHEPlan(
            scheme="BGV",
            ring_dimension=8192,
            slot_count=8192,
            batch_size=512,
            window_size=512,
            variants_per_ciphertext=512,
            rotation_indices=tuple(1 << index for index in range(9)),
            plaintext_modulus=2_147_352_577,
            fixed_point_scaling_factor=1_000_000,
        )
        bgv_result = validate_plan(
            bgv,
            samples=1,
            variants=512,
            genotype_min=0.0,
            genotype_max=2.0,
            sum_abs_weights=1.0,
            constraints=Constraints(),
        )
        self.assertIn(
            "bgv_ring_dimension_below_openfhe_minimum",
            {issue.failure_type for issue in bgv_result.issues},
        )

    def test_candidate_generation_is_deterministic(self) -> None:
        space = CandidateSpace(
            schemes=("CKKS",),
            ring_dimensions=(8192,),
            window_sizes=(512,),
            scaling_modulus_sizes=(40, 50),
            aggregation_strategies=("reduce_then_add",),
        )
        candidates = generate_candidates(1000, space)
        self.assertEqual(len(candidates), 2)
        self.assertEqual(candidates[0].rotation_indices, tuple(1 << i for i in range(9)))

    def test_plaintext_error_metrics(self) -> None:
        metrics = compare([1.0, 0.0], [1.1, 0.0])
        self.assertAlmostEqual(metrics.max_absolute_error, 0.1)
        self.assertTrue(metrics.finite)

    def test_yaml_configuration_without_required_third_party_parser(self) -> None:
        content = """\
workload:
  operation: prs_weighted_sum
  genotype_path: data/dosage.csv
  weight_path: data/weights.csv
  genotype_encoding: hard_call
constraints:
  security_bits: 128
  max_absolute_error: 1.0e-6
optimization:
  objective: evaluation_latency
  allowed_schemes:
    - CKKS
    - BFV
  top_k_validation: 3
hardware:
  profile: local_cpu
"""
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "spec.yaml"
            path.write_text(content, encoding="utf-8")
            spec = load_experiment_spec(path)
        self.assertEqual(spec.allowed_schemes, ("CKKS", "BFV"))
        self.assertEqual(spec.top_k_validation, 3)


if __name__ == "__main__":
    unittest.main()
