from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path
import platform
import statistics
import sys
from time import perf_counter


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from genome_he_pilot import load_vcf_dosage_subset, make_synthetic_dataset, prs_scores  # noqa: E402
from genome_he_pilot.openfhe_controlled_backend import (  # noqa: E402
    DEFAULT_OPENFHE_CKKS_SCALING_MOD_SIZE,
    DEFAULT_OPENFHE_FIXED_POINT_SCALE,
    DEFAULT_OPENFHE_PLAIN_MODULUS,
    DEFAULT_OPENFHE_RING_DIM,
    DEFAULT_OPENFHE_WINDOW_SNPS,
    openfhe_ckks_prs_benchmark,
    openfhe_fixed_point_prs_benchmark,
)


DEFAULT_VCF = (
    ROOT
    / "data"
    / "1000genomes"
    / "phase3"
    / "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
)
DEFAULT_RAW_OUTPUT = ROOT / "results" / "validation_inputs_openfhe_controlled_raw.csv"
DEFAULT_SUMMARY_OUTPUT = ROOT / "results" / "validation_inputs_openfhe_controlled_summary.csv"
DEFAULT_REPORT = ROOT / "submission" / "validation_inputs_openfhe_controlled_audit.md"

SCHEMES = ("CKKS", "BFV", "BGV")
SCHEME_ORDER = {scheme: index for index, scheme in enumerate(SCHEMES)}
SYNTHETIC_CONFIGS = (
    ("HE1", 8, 128, 20260701),
    ("HE2", 16, 256, 20260702),
    ("HE3", 32, 512, 20260703),
    ("HE4", 64, 1024, 20260704),
    ("HE5", 96, 2048, 20260705),
)
PUBLIC_CONFIGS = (
    ("G3", 16, 256, 20260733),
    ("G4", 32, 512, 20260734),
)
SUMMARY_METRICS = (
    "load_seconds",
    "plain_prs_seconds",
    "context_seconds",
    "encryption_seconds",
    "evaluation_seconds",
    "decryption_seconds",
    "serialization_seconds",
    "total_he_compute_seconds",
    "max_abs_error",
    "he_samples_per_second",
    "he_snp_products_per_second",
    "input_ciphertext_total_bytes",
    "input_ciphertext_mean_bytes_per_sample",
    "result_ciphertext_mean_bytes_per_sample",
    "public_context_bytes",
    "fixed_point_scale",
    "plain_modulus",
    "integer_dot_abs_bound",
    "safe_modulus_margin",
)


@dataclass(frozen=True)
class Check:
    name: str
    passed: bool
    detail: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, default=DEFAULT_VCF)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--dataset", choices=["all", "synthetic", "public"], default="all")
    parser.add_argument("--window-snps", type=int, default=DEFAULT_OPENFHE_WINDOW_SNPS)
    parser.add_argument("--ring-dim", type=int, default=DEFAULT_OPENFHE_RING_DIM)
    parser.add_argument("--plain-modulus", type=int, default=DEFAULT_OPENFHE_PLAIN_MODULUS)
    parser.add_argument("--fixed-point-scale", type=int, default=DEFAULT_OPENFHE_FIXED_POINT_SCALE)
    parser.add_argument("--ckks-scaling-mod-size", type=int, default=DEFAULT_OPENFHE_CKKS_SCALING_MOD_SIZE)
    parser.add_argument("--raw-output", type=Path, default=DEFAULT_RAW_OUTPUT)
    parser.add_argument("--summary-output", type=Path, default=DEFAULT_SUMMARY_OUTPUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def format_float(value: float | int | None) -> str:
    if value is None:
        return ""
    return f"{float(value):.12g}"


def sample_sd(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    return statistics.stdev(values)


def float_or_none(value: str) -> float | None:
    if value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def passfail(value: bool) -> str:
    return "PASS" if value else "FAIL"


def run_scheme(args: argparse.Namespace, dataset, scheme: str):
    if scheme == "CKKS":
        return openfhe_ckks_prs_benchmark(
            dataset.genotypes,
            dataset.weights,
            window_snps=args.window_snps,
            ring_dim=args.ring_dim,
            scaling_mod_size=args.ckks_scaling_mod_size,
        )
    return openfhe_fixed_point_prs_benchmark(
        dataset.genotypes,
        dataset.weights,
        scheme=scheme,
        window_snps=args.window_snps,
        ring_dim=args.ring_dim,
        plain_modulus=args.plain_modulus,
        fixed_point_scale=args.fixed_point_scale,
    )


def dataset_rows(args: argparse.Namespace):
    if args.dataset in {"all", "synthetic"}:
        for run, samples, snps, seed in SYNTHETIC_CONFIGS:
            start = perf_counter()
            dataset = make_synthetic_dataset(n_samples=samples, n_snps=snps, seed=seed)
            yield "synthetic", run, samples, snps, seed, "synthetic_binomial", dataset, perf_counter() - start
    if args.dataset in {"all", "public"}:
        for run, samples, snps, seed in PUBLIC_CONFIGS:
            start = perf_counter()
            dataset = load_vcf_dosage_subset(
                args.vcf,
                n_samples=samples,
                n_snps=snps,
                seed=seed,
                min_minor_allele_count=1,
            )
            yield "public_1000g_chr22", run, samples, snps, seed, "1000g_chr22_vcf", dataset, perf_counter() - start


def row_for_result(
    args: argparse.Namespace,
    dataset_group: str,
    run: str,
    samples: int,
    snps: int,
    seed: int,
    source: str,
    dataset,
    load_seconds: float,
    repeat: int,
    scheme: str,
) -> dict[str, str]:
    start = perf_counter()
    plain_scores = prs_scores(dataset.genotypes, dataset.weights)
    plain_seconds = perf_counter() - start

    result = run_scheme(args, dataset, scheme)
    max_error = max(abs(actual - expected) for actual, expected in zip(result.scores, plain_scores))
    total_he_seconds = result.encryption_seconds + result.evaluation_seconds + result.decryption_seconds

    return {
        "dataset_group": dataset_group,
        "run": run,
        "repeat": str(repeat),
        "samples": str(samples),
        "snps": str(snps),
        "seed": str(seed),
        "source": source,
        "scheme": scheme,
        "backend": f"OpenFHE {scheme}RNS" if scheme != "CKKS" else "OpenFHE CKKSRNS",
        "status": "ok",
        "load_seconds": format_float(load_seconds),
        "plain_prs_seconds": format_float(plain_seconds),
        "context_seconds": format_float(result.context_seconds),
        "encryption_seconds": format_float(result.encryption_seconds),
        "evaluation_seconds": format_float(result.evaluation_seconds),
        "decryption_seconds": format_float(result.decryption_seconds),
        "serialization_seconds": format_float(result.serialization_seconds),
        "total_he_compute_seconds": format_float(total_he_seconds),
        "max_abs_error": format_float(max_error),
        "he_samples_per_second": format_float(samples / result.evaluation_seconds if result.evaluation_seconds > 0 else 0.0),
        "he_snp_products_per_second": format_float(samples * snps / result.evaluation_seconds if result.evaluation_seconds > 0 else 0.0),
        "input_ciphertext_total_bytes": str(result.input_ciphertext_bytes),
        "input_ciphertext_mean_bytes_per_sample": format_float(result.input_ciphertext_bytes / samples if samples else 0.0),
        "result_ciphertext_total_bytes": str(result.result_ciphertext_bytes),
        "result_ciphertext_mean_bytes_per_sample": format_float(result.result_ciphertext_bytes / samples if samples else 0.0),
        "public_context_bytes": str(result.public_context_bytes),
        "window_snps": str(result.window_snps),
        "windows_per_sample": str(result.windows_per_sample),
        "input_ciphertexts": str(result.input_ciphertexts),
        "result_ciphertexts": str(result.result_ciphertexts),
        "ring_dim": str(result.ring_dim),
        "batch_size": str(result.batch_size),
        "plain_modulus": "" if result.plain_modulus is None else str(result.plain_modulus),
        "fixed_point_scale": "" if result.fixed_point_scale is None else str(result.fixed_point_scale),
        "ckks_scaling_mod_size": "" if result.ckks_scaling_mod_size is None else str(result.ckks_scaling_mod_size),
        "integer_weight_abs_sum": "" if result.integer_weight_abs_sum is None else str(result.integer_weight_abs_sum),
        "integer_dot_abs_bound": "" if result.integer_dot_abs_bound is None else str(result.integer_dot_abs_bound),
        "safe_modulus_margin": "" if result.safe_modulus_margin is None else str(result.safe_modulus_margin),
        "openfhe_version": version("openfhe"),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "note": "Same-backend OpenFHE validation-input PRS comparison",
    }


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("cannot write an empty CSV")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    groups: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[(row["dataset_group"], row["run"], row["scheme"])].append(row)

    summary_rows: list[dict[str, str]] = []
    for (dataset_group, run, scheme), group in sorted(
        groups.items(),
        key=lambda item: (item[0][0], item[0][1], SCHEME_ORDER[item[0][2]]),
    ):
        first = group[0]
        summary = {
            "dataset_group": dataset_group,
            "run": run,
            "samples": first["samples"],
            "snps": first["snps"],
            "source": first["source"],
            "scheme": scheme,
            "backend": first["backend"],
            "ring_dim": first["ring_dim"],
            "batch_size": first["batch_size"],
            "window_snps": first["window_snps"],
            "openfhe_version": first["openfhe_version"],
            "repeats": str(len(group)),
            "ok_repeats": str(sum(1 for row in group if row["status"] == "ok")),
            "status": "ok" if all(row["status"] == "ok" for row in group) else "mixed",
            "note": first["note"],
        }
        for metric in SUMMARY_METRICS:
            values = [value for row in group if (value := float_or_none(row.get(metric, ""))) is not None]
            if values:
                summary[f"{metric}_mean"] = format_float(statistics.mean(values))
                summary[f"{metric}_sd"] = format_float(sample_sd(values))
                summary[f"{metric}_min"] = format_float(min(values))
                summary[f"{metric}_max"] = format_float(max(values))
            else:
                summary[f"{metric}_mean"] = ""
                summary[f"{metric}_sd"] = ""
                summary[f"{metric}_min"] = ""
                summary[f"{metric}_max"] = ""
        summary_rows.append(summary)
    return summary_rows


def scheme_sets(summary_rows: list[dict[str, str]]) -> dict[tuple[str, str], set[str]]:
    found: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in summary_rows:
        found[(row["dataset_group"], row["run"])].add(row["scheme"])
    return found


def collect_checks(args: argparse.Namespace, rows: list[dict[str, str]], summary_rows: list[dict[str, str]]) -> list[Check]:
    expected_runs = []
    if args.dataset in {"all", "synthetic"}:
        expected_runs.extend(("synthetic", run) for run, *_ in SYNTHETIC_CONFIGS)
    if args.dataset in {"all", "public"}:
        expected_runs.extend(("public_1000g_chr22", run) for run, *_ in PUBLIC_CONFIGS)
    sets = scheme_sets(summary_rows)
    errors = [float(row["max_abs_error_mean"]) for row in summary_rows if row["status"] == "ok"]
    fixed_rows = [row for row in rows if row["scheme"] in {"BFV", "BGV"}]
    return [
        Check("Raw rows completed", len(rows) == len(expected_runs) * len(SCHEMES) * args.repeats, f"rows={len(rows)}"),
        Check("Summary rows completed", len(summary_rows) == len(expected_runs) * len(SCHEMES), f"rows={len(summary_rows)}"),
        Check("Each validation run has CKKS/BFV/BGV", all(sets.get(run) == set(SCHEMES) for run in expected_runs), str(dict(sets))),
        Check("All rows are ok", all(row["status"] == "ok" for row in rows), "all ok"),
        Check("All rows use OpenFHE 1.5.1 package", all(row["openfhe_version"].startswith("1.5.1") for row in rows), rows[0]["openfhe_version"] if rows else "missing"),
        Check("Maximum validation-input error below 1e-3", max(errors) < 1e-3 if errors else False, f"max={max(errors):.3g}" if errors else "missing"),
        Check("BFV/BGV modulus margins are positive", all(float(row["safe_modulus_margin"]) > 0 for row in fixed_rows), "positive"),
    ]


def write_report(args: argparse.Namespace, summary_rows: list[dict[str, str]], checks: list[Check]) -> None:
    failures = [check for check in checks if not check.passed]
    lines = [
        "# OpenFHE Controlled Validation-Input Scheme Comparison Audit",
        "",
        "## Verdict",
        "",
        "PASS: Synthetic and public validation inputs were benchmarked under CKKS, BFV, and BGV in the same OpenFHE runtime."
        if not failures
        else "FAIL: validation-input scheme comparison has unresolved issues.",
        "",
        "## Files",
        "",
        f"- Raw CSV: `{args.raw_output.relative_to(ROOT).as_posix()}`",
        f"- Summary CSV: `{args.summary_output.relative_to(ROOT).as_posix()}`",
        "",
        "## Checks",
        "",
        "| Check | Status | Detail |",
        "| --- | --- | --- |",
    ]
    lines.extend(f"| {check.name} | {passfail(check.passed)} | `{check.detail}` |" for check in checks)
    lines.extend(
        [
            "",
            "## Summary",
            "",
            "| Dataset | Run | Samples | SNPs | Scheme | Eval sec mean | Max error mean | Input KB/sample mean |",
            "| --- | --- | ---: | ---: | --- | ---: | ---: | ---: |",
        ]
    )
    for row in summary_rows:
        input_kb = float(row["input_ciphertext_mean_bytes_per_sample_mean"]) / 1024
        lines.append(
            f"| {row['dataset_group']} | {row['run']} | {row['samples']} | {row['snps']} | {row['scheme']} | "
            f"{float(row['evaluation_seconds_mean']):.4f} | {float(row['max_abs_error_mean']):.2e} | {input_kb:.1f} |"
        )
    if failures:
        lines.extend(["", "## Failures", ""])
        lines.extend(f"- {check.name}: `{check.detail}`" for check in failures)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    rows: list[dict[str, str]] = []
    for dataset_group, run, samples, snps, seed, source, dataset, load_seconds in dataset_rows(args):
        for repeat in range(1, args.repeats + 1):
            for scheme in SCHEMES:
                print(f"OpenFHE validation {dataset_group} {run} {scheme}: repeat={repeat}/{args.repeats}", flush=True)
                rows.append(
                    row_for_result(
                        args,
                        dataset_group,
                        run,
                        samples,
                        snps,
                        seed,
                        source,
                        dataset,
                        load_seconds,
                        repeat,
                        scheme,
                    )
                )

    summary_rows = summarize_rows(rows)
    write_csv(args.raw_output, rows)
    write_csv(args.summary_output, summary_rows)
    checks = collect_checks(args, rows, summary_rows)
    write_report(args, summary_rows, checks)
    failures = [check for check in checks if not check.passed]
    if failures:
        print("OpenFHE validation-input audit failures:")
        for check in failures:
            print(f"- {check.name}: {check.detail}")
        return 1
    print("OpenFHE validation-input scheme comparison audit passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
