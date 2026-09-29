# SCI Extension Plan

This note defines the next paper after the KCI pilot. The goal is to move from
"HE works for pilot-scale genomic risk scoring" to a publishable framework with
scale, deployment assumptions, and a clear technical contribution.

## Working SCI Claim

> A federated and packed homomorphic-encryption framework can compute genomic
> risk scores and selected association statistics over public or institutional
> genotype partitions while preserving raw genotype privacy and maintaining
> reproducible numerical accuracy.

## What Must Be New Beyond the KCI Paper

The KCI paper should remain a compact pilot. The SCI version needs at least one
of the following as the central contribution:

1. Packed CKKS implementation that reduces ciphertext count and amortizes PRS
   evaluation across samples or SNP blocks.
2. Federated protocol where multiple sites keep genotype data local and submit
   only encrypted genotype encodings or encrypted partial statistics.
3. A second exact-arithmetic backend for allele counts or burden statistics,
   most likely BFV/BGV through SEAL/OpenFHE if the Windows toolchain is stable.
4. Multi-chromosome or whole-genome benchmark using public genotypes and PGS
   Catalog scores, with a documented storage/runtime tradeoff.
5. Stronger reproducibility package with Docker/Conda, pinned data manifests,
   and independent verification scripts.

## Initial Federated Experiment Scaffold

The first runnable scaffold is:

```powershell
.\.venv\Scripts\python.exe -m scripts.run_federated_prs_simulation
```

It writes:

```text
results/sci_federated_prs_simulation.csv
```

The current scaffold simulates multiple institutions by partitioning synthetic
samples across sites. Each site is evaluated with a per-site CKKS context, and
the script reports total context, encryption, evaluation, decryption,
serialization time, ciphertext bytes, and maximum PRS error against plaintext.

This is not yet the final SCI contribution. It is a control experiment that
will become the baseline for packed CKKS and true multi-site protocol variants.

## Initial Packed CKKS Prototype

The first runnable packed prototype is:

```powershell
.\.venv\Scripts\python.exe -m scripts.run_packed_ckks_matrix
```

It writes:

```text
results/packed_ckks_prs_matrix.csv
```

The current prototype now has two packed modes:

1. Individual PRS recovery concatenates multiple samples into one CKKS vector
   and uses a plaintext block-diagonal matrix to recover one PRS score per
   packed sample. This reduces ciphertext bytes but is slower than the
   one-vector-per-sample baseline because TenSEAL plain matrix multiplication
   has substantial rotation overhead in this layout.
2. Aggregate PRS computes a block-level PRS sum with one encrypted dot product
   per packed block. This is not a substitute for per-sample PRS recovery, but
   it is directly useful for site/cohort-level statistics and federated
   aggregate reporting.

Current observation:

| Experiment | Samples | SNPs | Block size | Baseline eval sec | Packed individual eval sec | Packed aggregate eval sec | Aggregate speedup |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| PACK1 | 8 | 128 | 2 | 0.048539 | 0.395510 | 0.028190 | 1.72x |
| PACK2 | 16 | 256 | 4 | 0.113064 | 1.797892 | 0.033918 | 3.33x |
| PACK3 | 32 | 512 | 4 | 0.310357 | 8.091463 | 0.095137 | 3.26x |

Interpretation: individual packed PRS recovery is not yet faster, but packed
aggregate PRS already gives a speed and size signal for the SCI extension. The
next optimization problem should therefore be split into two tracks: diagonal
packing for individual PRS and federated aggregate statistics for site-level
reporting.

## Federated Packed Aggregate Prototype

The first site-level packed aggregate scaffold is:

```powershell
.\.venv\Scripts\python.exe -m scripts.run_federated_packed_aggregate
```

It writes:

```text
results/sci_federated_packed_aggregate.csv
```

Current observation:

| Experiment | Sites | Samples/Site | SNPs | Block size | Baseline eval sec | Packed aggregate eval sec | Aggregate speedup | Site-sum max error |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| SCI_FED_PACK_AGG_1 | 2 | 8 | 256 | 4 | 0.108465 | 0.034230 | 3.17x | 3.03e-6 |
| SCI_FED_PACK_AGG_2 | 4 | 8 | 512 | 4 | 0.256646 | 0.096363 | 2.66x | 5.46e-6 |
| SCI_FED_PACK_AGG_3 | 4 | 16 | 1,024 | 4 | 0.658983 | 0.200336 | 3.29x | 1.42e-5 |

Interpretation: this is the first result where packing and federation meet in
one runnable protocol. It targets site-level PRS sum/mean reporting rather than
individual PRS disclosure. At block size 4, it reduces input and result
ciphertext bytes by about 4x and gives a 2.66-3.29x evaluation speedup over the
per-site individual-score CKKS baseline.

## Real PGS Federated Packed Aggregate Prototype

The first public-data version is:

```powershell
.\.venv\Scripts\python.exe -m scripts.run_real_pgs_federated_packed_aggregate
```

It writes:

```text
results/real_pgs_federated_packed_aggregate.csv
```

Current observation on 1000 Genomes chr22 matched to CAD PGS004941:

| Matched SNPs | Sites | Samples | Block size | Slots used | Baseline eval sec | Packed aggregate eval sec | Aggregate speedup | Site-sum max error |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 512 | 4 | 32 | 4 | 2,048 | 0.307100 | 0.090971 | 3.38x | 3.02e-6 |
| 1,024 | 4 | 32 | 4 | 4,096 | 0.317537 | 0.099505 | 3.19x | 3.20e-6 |
| 2,048 | 4 | 32 | 2 | 4,096 | 0.345695 | 0.198642 | 1.74x | 9.69e-6 |
| 4,096 | 4 | 32 | 1 | 4,096 | 0.389816 | 0.393364 | 0.99x | 1.39e-5 |

Interpretation: the public-data run confirms the packed aggregate advantage at
512-1,024 SNPs and exposes the CKKS slot-capacity limit. Once `block_size *
matched_snps` reaches 4,096 slots, the method must reduce block size; the speed
and ciphertext-size advantage shrinks accordingly. This makes SNP-window
partitioning or larger parameter sets a concrete SCI design question rather
than a generic future-work statement.

## Real PGS Windowed Aggregate Prototype

The first slot-aware SNP-window version is:

```powershell
.\.venv\Scripts\python.exe -m scripts.run_real_pgs_windowed_federated_packed_aggregate
```

It writes:

```text
results/real_pgs_windowed_federated_packed_aggregate.csv
```

Current observation on 1000 Genomes chr22 matched to CAD PGS004941:

| Matched SNPs | Direct block | Windowed block | Windows/Site | Baseline eval sec | Direct aggregate eval sec | Windowed aggregate eval sec | Windowed speedup vs baseline | Windowed speedup vs direct |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 512 | 4 | 4 | 1 | 0.298279 | 0.093640 | 0.100391 | 2.97x | 0.93x |
| 1,024 | 4 | 4 | 1 | 0.372227 | 0.110889 | 0.091799 | 4.05x | 1.21x |
| 2,048 | 2 | 4 | 2 | 0.407743 | 0.220941 | 0.221812 | 1.84x | 1.00x |
| 4,096 | 1 | 4 | 4 | 0.445004 | 0.442281 | 0.435843 | 1.02x | 1.01x |

Interpretation: slot-aware windows prove that the protocol can keep block size
4 beyond the full-vector slot limit, but this first implementation evaluates
each site-window independently. It motivated the batched context-reuse variant
below.

## Real PGS Batched Windowed Aggregate Prototype

The context-reused batched SNP-window version is:

```powershell
.\.venv\Scripts\python.exe -m scripts.run_real_pgs_batched_windowed_federated_packed_aggregate
```

It writes:

```text
results/real_pgs_batched_windowed_federated_packed_aggregate.csv
```

Current observation on 1000 Genomes chr22 matched to CAD PGS004941:

| Matched SNPs | Windows/Site | Independent eval sec | Batched eval sec | Compute speedup vs independent | Context reduction | Site-sum max error |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 512 | 1 | 0.093102 | 0.093244 | 2.83x | 3.78x | 5.29e-6 |
| 1,024 | 1 | 0.100839 | 0.099910 | 2.85x | 3.89x | 2.09e-6 |
| 2,048 | 2 | 0.176353 | 0.200479 | 4.02x | 7.91x | 2.93e-6 |
| 4,096 | 4 | 0.388682 | 0.393431 | 5.44x | 15.59x | 2.16e-5 |

Interpretation: reusing one CKKS context across all site-window aggregate
ciphertexts strongly reduces context overhead and public context bytes. It does
not yet reduce the count of encrypted dot products, so evaluation-only speed is
similar to independent windows. The next SCI method step should target either
fused/diagonal aggregate evaluation or rotation-aware individual-score packing.

## Experiment Ladder

| Stage | Dataset | Method | Expected paper value |
| --- | --- | --- | --- |
| KCI pilot | Synthetic, chr22 1000G, CAD PGS chr22 | One sample vector per CKKS ciphertext | Demonstrates feasibility and reproducibility |
| SCI baseline | Synthetic multi-site partitions | Per-site CKKS contexts | Establishes federated baseline overhead |
| SCI packed baseline | Synthetic packed sample blocks | Naive CKKS plaintext matrix multiply plus aggregate dot product | Shows per-sample bottleneck and aggregate speed advantage |
| SCI method | Synthetic and 1000G multi-site partitions | Packed CKKS site-level aggregate statistics with context-reused SNP windows, dot-product fusion, and optimized individual-score packing | Main technical novelty |
| SCI validation | Multi-chromosome PGS Catalog scores | Packed CKKS plus plaintext comparison | Shows practical genomic scale |
| SCI extension | Case/control allele or burden statistic | Exact BFV/BGV or CKKS approximation | Broadens beyond PRS |

## Decision Gates

Proceed to an SCI manuscript only after these are true:

- Packed CKKS gives a measurable speed or size advantage over the current
  one-vector-per-sample backend for either per-sample PRS recovery or a clearly
  scoped aggregate genomic statistic.
- At least one public-data benchmark uses more than chromosome 22 or uses a
  clearly justified multi-chromosome subset.
- The security model states who holds secret keys, what the evaluator sees, and
  what collusion assumptions are excluded.
- The reproducibility package can rebuild all figures/tables from raw public
  data or from documented cache files.
- The paper compares against plaintext and non-packed HE baselines.

## Likely SCI Manuscript Structure

1. Introduction and clinical privacy motivation.
2. Related work: HE PRS, privacy-preserving GWAS, federated genomics, secure
   multiparty computation baselines.
3. Threat model and system architecture.
4. Packed encoding and protocol design.
5. PRS and genomic statistic algorithms.
6. Experimental setup and public datasets.
7. Runtime, memory, ciphertext size, and numerical accuracy results.
8. Discussion: deployment, limits, and controlled-access data.
9. Conclusion.
