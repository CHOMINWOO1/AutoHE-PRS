"""SQLite cache for repeated VCF dosage lookups."""

from __future__ import annotations

from dataclasses import dataclass
import gzip
from pathlib import Path
import sqlite3

from .pgs import PGSScoringFile
from .pgs_vcf import PGSVCFMatchResult
from .vcf import (
    VCFGenomicDataset,
    VariantRecord,
    _is_biallelic_snp,
    _parse_gt_dosage,
    _select_vcf_samples,
)


@dataclass(frozen=True)
class VCFCacheBuildResult:
    cache_path: Path
    samples: int
    variants: int


def build_vcf_dosage_cache(
    vcf_path: Path,
    cache_path: Path,
    n_samples: int = 32,
    sample_ids: list[str] | None = None,
    chrom: str = "22",
) -> VCFCacheBuildResult:
    """Build a SQLite cache of biallelic SNP alternate allele dosages."""

    target_chrom = chrom.removeprefix("chr")
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    if cache_path.exists():
        cache_path.unlink()

    connection = sqlite3.connect(cache_path)
    try:
        connection.execute("PRAGMA journal_mode=OFF")
        connection.execute("PRAGMA synchronous=OFF")
        connection.execute("PRAGMA temp_store=MEMORY")
        connection.execute(
            "CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)"
        )
        connection.execute(
            """
            CREATE TABLE variants (
                chrom TEXT NOT NULL,
                position INTEGER NOT NULL,
                variant_id TEXT NOT NULL,
                ref TEXT NOT NULL,
                alt TEXT NOT NULL,
                allele_frequency REAL NOT NULL,
                dosages TEXT NOT NULL,
                PRIMARY KEY (chrom, position, ref, alt)
            )
            """
        )

        selected_sample_ids: list[str] | None = None
        sample_field_indices: list[int] | None = None
        rows: list[tuple[str, int, str, str, str, float, str]] = []
        inserted = 0

        with gzip.open(vcf_path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    header = line.rstrip("\n").split("\t")
                    selected_sample_ids, sample_field_indices = _select_vcf_samples(
                        vcf_sample_ids=header[9:],
                        n_samples=n_samples,
                        requested_sample_ids=sample_ids,
                    )
                    continue
                if line.startswith("#"):
                    continue
                if selected_sample_ids is None or sample_field_indices is None:
                    raise ValueError("VCF header with sample IDs was not found")

                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10:
                    continue
                vcf_chrom = fields[0].removeprefix("chr")
                if vcf_chrom != target_chrom:
                    continue

                position = int(fields[1])
                variant_id, ref, alt = fields[2], fields[3].upper(), fields[4].upper()
                if not _is_biallelic_snp(ref, alt):
                    continue

                format_fields = fields[8].split(":")
                try:
                    gt_index = format_fields.index("GT")
                except ValueError:
                    continue

                dosages: list[int] = []
                missing = False
                sample_fields = fields[9:]
                for sample_index in sample_field_indices:
                    dosage = _parse_gt_dosage(sample_fields[sample_index], gt_index)
                    if dosage is None:
                        missing = True
                        break
                    dosages.append(dosage)
                if missing:
                    continue

                allele_frequency = sum(dosages) / (2 * len(dosages))
                rows.append(
                    (
                        vcf_chrom,
                        position,
                        variant_id,
                        ref,
                        alt,
                        allele_frequency,
                        ",".join(str(value) for value in dosages),
                    )
                )
                if len(rows) >= 5000:
                    connection.executemany(
                        "INSERT OR IGNORE INTO variants VALUES (?, ?, ?, ?, ?, ?, ?)",
                        rows,
                    )
                    inserted += len(rows)
                    rows.clear()

        if selected_sample_ids is None:
            raise ValueError("VCF header with sample IDs was not found")
        if rows:
            connection.executemany(
                "INSERT OR IGNORE INTO variants VALUES (?, ?, ?, ?, ?, ?, ?)",
                rows,
            )
            inserted += len(rows)

        connection.executemany(
            "INSERT INTO metadata VALUES (?, ?)",
            [
                ("source_path", str(vcf_path)),
                ("chrom", target_chrom),
                ("sample_count", str(len(selected_sample_ids))),
                ("sample_ids", ",".join(selected_sample_ids)),
            ],
        )
        connection.commit()
        variant_count = connection.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
    finally:
        connection.close()

    return VCFCacheBuildResult(
        cache_path=cache_path,
        samples=n_samples,
        variants=int(variant_count if "variant_count" in locals() else inserted),
    )


def load_pgs_matched_cached_dataset(
    cache_path: Path,
    scoring: PGSScoringFile,
    chrom: str = "22",
    max_snps: int | None = None,
) -> PGSVCFMatchResult:
    """Load effect-allele dosage for PGS variants from a SQLite VCF cache."""

    target_chrom = chrom.removeprefix("chr")
    weights_by_position = {
        weight.position: weight
        for weight in scoring.weights
        if weight.chrom == target_chrom
    }
    if not weights_by_position:
        raise ValueError(f"no PGS weights found on chromosome {target_chrom}")

    connection = sqlite3.connect(cache_path)
    try:
        metadata = dict(connection.execute("SELECT key, value FROM metadata").fetchall())
        sample_ids = metadata["sample_ids"].split(",")
        source_path = Path(metadata.get("source_path", ""))

        variant_columns: list[list[int]] = []
        weights: list[float] = []
        variants: list[VariantRecord] = []
        skipped_allele_mismatch = 0

        for position, weight in weights_by_position.items():
            if max_snps is not None and len(weights) >= max_snps:
                break
            rows = connection.execute(
                """
                SELECT chrom, position, variant_id, ref, alt, allele_frequency, dosages
                FROM variants
                WHERE chrom = ? AND position = ?
                """,
                (target_chrom, position),
            ).fetchall()
            for row in rows:
                if max_snps is not None and len(weights) >= max_snps:
                    break
                row_chrom, row_position, variant_id, ref, alt, allele_frequency, dosage_text = row
                alt_dosages = [int(value) for value in dosage_text.split(",")]
                if weight.effect_allele == alt:
                    effect_dosages = alt_dosages
                    effect_frequency = float(allele_frequency)
                elif weight.effect_allele == ref:
                    effect_dosages = [2 - dosage for dosage in alt_dosages]
                    effect_frequency = sum(effect_dosages) / (2 * len(effect_dosages))
                else:
                    skipped_allele_mismatch += 1
                    continue

                variant_columns.append(effect_dosages)
                weights.append(weight.effect_weight)
                variants.append(
                    VariantRecord(
                        chrom=row_chrom,
                        position=int(row_position),
                        variant_id=variant_id if variant_id != "." else weight.rsid,
                        ref=ref,
                        alt=alt,
                        allele_frequency=effect_frequency,
                    )
                )
    finally:
        connection.close()

    if not variant_columns:
        raise ValueError("no PGS variants matched the VCF cache with compatible alleles")

    genotypes = [
        [variant_columns[j][i] for j in range(len(variant_columns))]
        for i in range(len(sample_ids))
    ]
    raw_scores = [sum(g * w for g, w in zip(row, weights)) for row in genotypes]
    threshold = sorted(raw_scores)[len(raw_scores) // 2]
    phenotypes = [int(score >= threshold) for score in raw_scores]

    return PGSVCFMatchResult(
        dataset=VCFGenomicDataset(
            genotypes=genotypes,
            weights=weights,
            phenotypes=phenotypes,
            allele_frequencies=[variant.allele_frequency for variant in variants],
            sample_ids=sample_ids,
            variants=variants,
            source_path=source_path,
        ),
        matched_weights=len(weights),
        skipped_allele_mismatch=skipped_allele_mismatch,
        requested_weights_on_chrom=len(weights_by_position),
    )
