# PGS Catalog Integration Plan

## Why This Is Separate

The current KCI pilot uses synthetic sparse weights. This is honest and useful
for validating encrypted computation, but a stronger public-data validation
needs a real PGS Catalog scoring file.

## Required Matching Logic

The 1000 Genomes parser currently encodes genotype as alternate allele dosage.
PGS Catalog scoring files encode the dosage of the `effect_allele`.

For each matched variant:

- If `effect_allele == ALT`, use `ALT dosage`.
- If `effect_allele == REF`, use `2 - ALT dosage`.
- Otherwise, skip or handle strand/allele harmonization explicitly.

The first real-PRS implementation should avoid ambiguous strand cases and report
how many variants were matched, skipped, or flipped.

## Added Utilities

Inspect a downloaded PGS scoring file:

```powershell
.\.venv\Scripts\python.exe scripts\inspect_pgs_file.py data\pgs_catalog\PGSXXXXXX\PGSXXXXXX_hmPOS_GRCh37.txt.gz --chrom 22
```

Download a PGS score by ID:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\download_pgs_score.ps1 -PgsId PGSXXXXXX -Build GRCh37
```

Run the current smoke test:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\download_pgs_score.ps1 -PgsId PGS000348 -Build formatted
.\.venv\Scripts\python.exe scripts\run_real_pgs_chr22.py --backend tenseal
```

Current smoke-test result:

- PGS: `PGS000348`
- Trait: prostate cancer
- Genome build: GRCh37
- chr22 matched variants: 1
- CKKS max absolute error: about `3.90e-8`
- Limitation: not enough chr22 variants for a meaningful PRS validation table

## Current Larger CAD Candidate

`PGS004941` was selected as a stronger real-PRS validation candidate:

- Trait: coronary artery disease
- Genome build: GRCh37
- Catalog metadata variants: 3,711,629
- chr22 scoring-file weights: 51,391
- chr22 VCF/cache matched SNPs: 51,227

Results:

```text
results/real_pgs_chr22_pgs004941_plain.csv
results/real_pgs_chr22_pgs004941_ckks_matrix.csv
```

The CKKS pilot matrix used `--max-snps 512 1024 2048 4096`:

- samples: 32
- matched SNPs: 512-4,096
- encrypted evaluation: about 0.340-0.499 seconds
- max absolute error: about `2.49e-6`

## Recommended Next Trait

Use coronary artery disease or type 2 diabetes for the first real PRS validation.
Both are clinically recognizable and have established PRS literature. Choose a
PGS with:

- GRCh37 harmonized scoring file
- SNP-level positions
- additive effect weights
- enough chromosome 22 variants for a small validation
- license terms suitable for research use

## SQLite VCF Cache

The first implementation scanned the gzipped chr22 VCF sequentially. A SQLite
dosage cache has now been added:

```powershell
.\.venv\Scripts\python.exe scripts\build_vcf_cache.py --samples 32
```

Current cache:

- `data/cache/1000g_chr22_samples32.sqlite`
- biallelic chr22 SNP rows: 1,055,454
- size: about 118 MB
- build time: about 71 seconds

This cache is good enough for KCI experiments. A future SCI implementation
should still consider tabix/htslib or a more compact columnar cache.

## Provenance Check

The public-data provenance audit ties `PGS004941`, the 1000 Genomes chr22 files,
the SQLite cache, and the current CAD result CSVs together:

```powershell
.\.venv\Scripts\python.exe -m scripts.check_public_data_provenance
```
