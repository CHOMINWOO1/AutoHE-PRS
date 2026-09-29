# Data Sources for Public Validation

## Do Not Download Everything Yet

The encrypted experiment matrix can be completed with synthetic data. Public data
is needed after the KCI method and benchmark table are stable.

## First Public Dataset to Download

Use a small 1000 Genomes Phase 3 chromosome subset first:

- 1000 Genomes Phase 3 chr22 VCF
- chr22 tabix index
- sample panel file

Recommended files:

```text
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
```

Approximate size:

- chr22 VCF: 196 MB
- tabix index: 35 KB
- sample panel: 54 KB

Why chr22 first:

- It is the smallest autosome in the Phase 3 VCF set.
- It is enough to validate VCF parsing, genotype dosage extraction, and allele
  count computation.
- It avoids downloading multi-GB whole-genome files before the method is stable.

## PRS Weight Data

For real PRS weights, use the PGS Catalog after choosing a trait. PGS scoring
files provide columns such as:

- `rsID`
- `chr_name`
- `chr_position`
- `effect_allele`
- `other_allele`
- `effect_weight`

Start with harmonized GRCh37 scoring files when matching the 1000 Genomes Phase
3 VCF release. If a selected PGS has too few chr22 variants, use it later for a
larger multi-chromosome run and keep chr22 validation with synthetic weights.

## Provenance Audit

The local public-data provenance audit records downloaded file sizes, SHA-256
values, SQLite cache table counts, and result-lineage checks without
redistributing raw genotype or scoring files:

```powershell
.\.venv\Scripts\python.exe -m scripts.check_public_data_provenance
```

Raw 1000 Genomes files, downloaded PGS scoring files, and local SQLite caches
are kept outside the anonymous review ZIP. The submission bundle includes
scripts, result CSVs, checksums, and provenance reports so reviewers can rebuild
or inspect the public-data path without redistributing large public inputs.

## Suggested Project Sequence

1. Finish synthetic encrypted benchmark table.
2. Download 1000 Genomes chr22 files.
3. Parse a small subset: 32 to 128 samples, 512 to 4,096 SNPs.
4. Re-run allele count and PRS using public genotype dosage.
5. Add a selected PGS Catalog scoring file after choosing the disease/trait.
