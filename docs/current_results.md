# Current Results Snapshot

## Synthetic OpenFHE Validation Matrix

Result file:

```text
results/validation_inputs_openfhe_controlled_summary.csv
```

Key result:

| Experiment | Samples | SNPs | Scheme | Encrypted Eval Seconds | Max Abs Error | Input KB/Sample |
| --- | ---: | ---: | --- | ---: | ---: | ---: |
| HE1 | 8 | 128 | CKKS | 1.2247 | 1.96e-12 | 769.0 |
| HE1 | 8 | 128 | BFV | 1.2873 | 5.00e-8 | 512.8 |
| HE1 | 8 | 128 | BGV | 1.1942 | 5.00e-8 | 769.0 |
| HE5 | 96 | 2,048 | CKKS | 20.2933 | 3.96e-12 | 769.0 |
| HE5 | 96 | 2,048 | BFV | 24.9972 | 2.79e-7 | 512.8 |
| HE5 | 96 | 2,048 | BGV | 20.7453 | 2.79e-7 | 769.0 |

Interpretation:

- Synthetic validation now compares CKKS, BFV, and BGV in the same OpenFHE Docker runtime.
- CKKS error is near `1e-12`; BFV/BGV error is fixed-point quantization error and stays below `2.8e-7`.
- The current implementation encrypts one sample vector per ciphertext with one 4,096-SNP window for these synthetic sizes.

## 1000 Genomes chr22 Public Validation

Result file:

```text
results/validation_inputs_openfhe_controlled_summary.csv
```

VCF filtering:

- Source: 1000 Genomes Phase 3 chr22 VCF
- Variant type: biallelic SNP only
- Genotype encoding: alternate allele dosage `0/1/2`
- Minimum minor allele count in selected samples: 1
- PRS weights: reproducible synthetic sparse weights

Key result:

| Experiment | Samples | SNPs | Scheme | Load Seconds | Eval Seconds | Max Abs Error | Input KB/Sample |
| --- | ---: | ---: | --- | ---: | ---: | ---: | ---: |
| G3 | 16 | 256 | CKKS | 0.188345 | 2.6052 | 3.33e-12 | 769.0 |
| G3 | 16 | 256 | BFV | 0.188345 | 3.0204 | 1.38e-7 | 512.8 |
| G3 | 16 | 256 | BGV | 0.188345 | 2.6347 | 1.38e-7 | 769.0 |
| G4 | 32 | 512 | CKKS | 0.352994 | 5.4655 | 3.10e-12 | 769.0 |
| G4 | 32 | 512 | BFV | 0.352994 | 6.5825 | 2.23e-7 | 512.8 |
| G4 | 32 | 512 | BGV | 0.352994 | 5.6408 | 2.23e-7 | 769.0 |

Interpretation:

- The VCF parser successfully extracts real genotype dosage from chr22.
- The same encrypted PRS workflow works on public genotype dosage under CKKS, BFV, and BGV.
- CKKS error is near `1e-12`; BFV/BGV error stays below `2.3e-7` in the public validation rows.
- The public validation now uses the same OpenFHE backend as the real-PGS comparison.

## Environment Report

Result files:

```text
results/environment_report.csv
results/environment_report.md
```

Key environment:

- OS: Windows 10
- Logical CPU count: 24
- RAM: about 68.1 GB
- Controlled validation/comparison runtime: Docker image `genome-he-openfhe:1.5.1`
- Docker server/context: Docker Desktop 29.5.3, `desktop-linux`
- Docker daemon kernel: Linux `6.18.33.1-microsoft-standard-WSL2`
- Docker cgroup version: 2
- Docker-reported CPU count: 24
- Docker-reported memory visible to containers: 33,367,441,408 bytes, about 31.1 GiB
- Container CPU quota/cpuset: no explicit quota, `cpu.max=max 100000`; effective cpuset `0-23`
- Container memory limit: no explicit per-container limit, `memory.max=max`; `/proc/meminfo` MemTotal 32,585,392 kB, about 31.1 GiB
- Runtime OS/Python: Linux WSL2, Python 3.10.20
- OpenFHE-Python: 1.5.1.0.22.4
- Common OpenFHE ring dimension: `N=16384`
- Common batch/window size: 4,096 SNPs
- CKKS approximate `log2(Q)=129.087`
- BFV approximate `log2(Q)=120.000`, plaintext modulus `p=2147352577`, fixed-point scale `S=10000000`
- BGV approximate `log2(Q)=114.000`, plaintext modulus `p=2147352577`, fixed-point scale `S=10000000`

## Real PGS Smoke Test

Result file:

```text
results/real_pgs_chr22_summary.csv
```

PGS scoring file:

- PGS ID: `PGS000348`
- Reported trait: prostate cancer
- Genome build: GRCh37
- Source file: `data/pgs_catalog/PGS000348/PGS000348.txt.gz`

Key result:

| Samples | Chromosome | Matched SNPs | Backend | VCF Match Seconds | Eval Seconds | Max Abs Error |
| ---: | --- | ---: | --- | ---: | ---: | ---: |
| 32 | 22 | 1 | CKKS | 46.273172 | 0.027215 | 3.90e-8 |

Interpretation:

- The real PGS parser works on a downloaded PGS Catalog scoring file.
- Effect-allele dosage matching with 1000 Genomes VCF works.
- `PGS000348` has only one chr22 variant, so this is a smoke test rather than
  a strong biological PRS validation.
- Sequential VCF matching took about 46 seconds for this smoke test.
- A SQLite dosage cache was later added to avoid repeated full-VCF scans.

## Real CAD PGS Validation

Candidate files:

```text
data/pgs_catalog/PGS004941/PGS004941.txt.gz
data/cache/1000g_chr22_samples32.sqlite
```

PGS scoring file:

- PGS ID: `PGS004941`
- Reported trait: coronary artery disease
- Genome build: GRCh37
- Total variants in metadata: 3,711,629
- chr22 weights in scoring file: 51,391

Key result:

| Run | Samples | Matched SNPs | Backend | Match Mode | Match Seconds | Eval Seconds | Max Abs Error |
| --- | ---: | ---: | --- | --- | ---: | ---: | ---: |
| Full chr22 PRS | 32 | 51,227 | plaintext | SQLite cache | 1.936398 |  |  |
| CKKS subset | 32 | 512 | CKKS | SQLite cache | 0.043093 | 0.340483 | 2.49e-6 |
| CKKS subset | 32 | 1,024 | CKKS | SQLite cache | 0.044361 | 0.446604 | 1.46e-6 |
| CKKS subset | 32 | 2,048 | CKKS | SQLite cache | 0.088446 | 0.499111 | 1.93e-6 |
| CKKS subset | 32 | 4,096 | CKKS | SQLite cache | 0.157311 | 0.424159 | 1.08e-6 |

Interpretation:

- The SQLite chr22 cache reduced repeated PGS-VCF matching from tens of seconds
  to sub-second or low-second lookup.
- A real CAD PGS scoring file can be matched against the 1000 Genomes chr22 VCF.
- The CKKS workflow works on real PGS-weight subsets from 512 to 4,096 matched SNPs, not only synthetic weights.
- Full 51k-SNP encrypted CAD PRS is now evaluated for the allele-compatible chromosome-22 subset with chunked CKKS, BFV, and BGV. Genome-wide PGS004941 remains future work because the full score has millions of variants.

## Real CAD PGS HE Scheme Comparison

Result files:

```text
results/real_pgs_chr22_pgs004941_openfhe_controlled_raw.csv
results/real_pgs_chr22_pgs004941_openfhe_controlled_summary.csv
results/openfhe_controlled_runtime_parameters.csv
submission/openfhe_controlled_scheme_comparison_audit.md
submission/openfhe_controlled_full_run.log
```

Controlled OpenFHE parameters:

| Scheme | Ring dimension N | Batch/window SNPs | Multiplicative depth | log2(Q) approx. | Plain modulus | Fixed-point scale |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| CKKS | 16,384 | 4,096 | 1 | 129.087 |  |  |
| BFV | 16,384 | 4,096 | 1 | 120.000 | 2,147,352,577 | 10,000,000 |
| BGV | 16,384 | 4,096 | 1 | 114.000 | 2,147,352,577 | 10,000,000 |

Key result:

| SNPs | Scheme | Status | Eval Seconds | Max Abs Error | Input KB/Sample |
| ---: | --- | --- | ---: | ---: | ---: |
| 512 | CKKS | ok | 6.6320 +/- 0.7558 | 3.29e-12 +/- 6.24e-13 | 769.0 |
| 512 | BFV | ok | 8.2608 +/- 0.1609 | 1.31e-6 +/- 0.00e+0 | 512.8 |
| 512 | BGV | ok | 6.4171 +/- 0.3058 | 1.31e-6 +/- 0.00e+0 | 769.0 |
| 1,024 | CKKS | ok | 6.7677 +/- 0.2876 | 3.21e-12 +/- 4.14e-13 | 769.0 |
| 1,024 | BFV | ok | 8.0806 +/- 0.6203 | 6.46e-7 +/- 0.00e+0 | 512.8 |
| 1,024 | BGV | ok | 7.0961 +/- 0.0458 | 6.46e-7 +/- 0.00e+0 | 769.0 |
| 2,048 | CKKS | ok | 7.6904 +/- 0.3177 | 3.68e-12 +/- 1.43e-13 | 769.0 |
| 2,048 | BFV | ok | 9.0626 +/- 0.3313 | 2.16e-6 +/- 0.00e+0 | 512.8 |
| 2,048 | BGV | ok | 7.8102 +/- 0.1941 | 2.16e-6 +/- 0.00e+0 | 769.0 |
| 4,096 | CKKS | ok | 8.3866 +/- 0.2379 | 3.22e-12 +/- 8.92e-13 | 769.0 |
| 4,096 | BFV | ok | 9.5287 +/- 0.4469 | 3.14e-6 +/- 0.00e+0 | 512.8 |
| 4,096 | BGV | ok | 8.1034 +/- 0.8457 | 3.14e-6 +/- 0.00e+0 | 769.0 |

Interpretation:

- CKKS, BFV, and BGV were rerun on the same PGS004941 chromosome-22 matched subsets with three repeats in one Docker/OpenFHE runtime.
- BFV and BGV use fixed-point integer weights with scale `10,000,000`; their reported errors include weight quantization.
- In this controlled OpenFHE single-window range, BFV was slower than CKKS and BGV, while CKKS and BGV had similar evaluation times.
- Microsoft SEAL itself supports BFV, BGV, and CKKS; the reason for using OpenFHE here is that the TenSEAL 0.3.16 Python vector API used in the early baseline did not expose the BGV vector interface needed for this benchmark.

## Real CAD PGS Chunked HE Scheme Comparison

Result files:

```text
results/real_pgs_chr22_pgs004941_openfhe_controlled_raw.csv
results/real_pgs_chr22_pgs004941_openfhe_controlled_summary.csv
submission/openfhe_controlled_scheme_comparison_audit.md
```

Key result:

| SNPs | Scheme | Windows/Sample | Eval Seconds | Max Abs Error | Input MB/Sample |
| ---: | --- | ---: | ---: | ---: | ---: |
| 8,192 | CKKS | 2 | 15.5666 +/- 1.1430 | 6.35e-12 +/- 5.71e-13 | 1.502 |
| 8,192 | BFV | 2 | 18.9087 +/- 1.1383 | 4.89e-6 +/- 0.00e+0 | 1.002 |
| 8,192 | BGV | 2 | 15.7052 +/- 0.1316 | 4.89e-6 +/- 0.00e+0 | 1.502 |
| 16,384 | CKKS | 4 | 32.4911 +/- 0.7576 | 5.22e-12 +/- 4.85e-13 | 3.004 |
| 16,384 | BFV | 4 | 40.6016 +/- 2.9875 | 8.55e-6 +/- 0.00e+0 | 2.003 |
| 16,384 | BGV | 4 | 31.7587 +/- 0.8472 | 8.55e-6 +/- 0.00e+0 | 3.004 |
| 32,768 | CKKS | 8 | 70.7985 +/- 4.4521 | 1.00e-11 +/- 8.91e-13 | 6.008 |
| 32,768 | BFV | 8 | 87.9265 +/- 4.2657 | 9.66e-6 +/- 0.00e+0 | 4.007 |
| 32,768 | BGV | 8 | 74.0597 +/- 2.8591 | 9.66e-6 +/- 0.00e+0 | 6.007 |
| 51,227 | CKKS | 13 | 107.4396 +/- 6.6052 | 1.27e-11 +/- 3.60e-13 | 9.762 |
| 51,227 | BFV | 13 | 129.0441 +/- 7.9971 | 1.04e-5 +/- 0.00e+0 | 6.511 |
| 51,227 | BGV | 13 | 103.9833 +/- 9.0158 | 1.04e-5 +/- 0.00e+0 | 9.762 |

Interpretation:

- The full allele-compatible chromosome-22 subset is now evaluated under all three schemes using 4,096-SNP chunks in the same OpenFHE runtime.
- At 51,227 SNPs, mean evaluation time is 107.4396 seconds for CKKS, 129.0441 seconds for BFV, and 103.9833 seconds for BGV.
- BFV and BGV have identical fixed-point quantization errors because they use the same quantized weights and plaintext modulus safety checks.

## SCI Federated PRS Scaffold

Result file:

```text
results/sci_federated_prs_simulation.csv
```

Key result:

| Experiment | Sites | Samples/Site | SNPs | Encrypted Eval Seconds | Max Abs Error |
| --- | ---: | ---: | ---: | ---: | ---: |
| SCI_FED_1 | 2 | 8 | 256 | 0.136992 | 1.58e-6 |
| SCI_FED_2 | 4 | 8 | 512 | 0.289624 | 2.48e-6 |
| SCI_FED_3 | 4 | 16 | 1,024 | 0.851901 | 1.40e-6 |

Interpretation:

- This is an SCI extension scaffold, not a KCI manuscript claim.
- Synthetic samples are partitioned across simulated sites.
- Each site uses a per-site CKKS context, which gives a conservative baseline
  before packed CKKS optimization.
- Accuracy remains in the same approximate range as the single-site pilot.

## SCI Packed CKKS Prototype

Result file:

```text
results/packed_ckks_prs_matrix.csv
```

Key result:

| Experiment | Samples | SNPs | Block Size | Baseline Eval Seconds | Packed Individual Eval Seconds | Packed Aggregate Eval Seconds | Aggregate Speedup vs Baseline | Packed Individual Max Error | Packed Aggregate Max Error |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| PACK1 | 8 | 128 | 2 | 0.048539 | 0.395510 | 0.028190 | 1.72x | 1.04e-8 | 8.22e-7 |
| PACK2 | 16 | 256 | 4 | 0.113064 | 1.797892 | 0.033918 | 3.33x | 2.35e-8 | 1.07e-6 |
| PACK3 | 32 | 512 | 4 | 0.310357 | 8.091463 | 0.095137 | 3.26x | 5.97e-8 | 7.37e-7 |

Size signal:

| Experiment | Baseline Input Bytes | Packed Individual Input Bytes | Packed Aggregate Input Bytes | Baseline Result Bytes | Packed Individual Result Bytes | Packed Aggregate Result Bytes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| PACK1 | 2,674,290 | 1,337,226 | 1,337,089 | 1,883,027 | 941,485 | 941,390 |
| PACK2 | 5,350,545 | 1,336,772 | 1,337,578 | 3,765,819 | 941,268 | 941,444 |
| PACK3 | 10,698,093 | 2,674,725 | 2,674,944 | 7,531,479 | 1,882,949 | 1,882,949 |

Interpretation:

- The first packed prototype is a working SCI baseline.
- It packs multiple samples per CKKS vector and returns one PRS score per sample.
- Per-sample packed PRS recovery is still slower than the simple backend in the
  current matrix layout.
- The new packed aggregate path computes a block-level PRS sum with one
  encrypted dot product per block; it is faster than the simple backend in all
  three current experiments while keeping aggregate error around `1e-6`.
- This separates two SCI directions: optimized diagonal/rotation scheduling for
  per-sample PRS, and federated site-level aggregate statistics where packed
  dot products already show a speed and size advantage.

## SCI Federated Packed Aggregate Prototype

Result file:

```text
results/sci_federated_packed_aggregate.csv
```

Key result:

| Experiment | Sites | Samples/Site | SNPs | Block Size | Baseline Eval Seconds | Packed Aggregate Eval Seconds | Aggregate Speedup | Aggregate Site-Sum Max Error |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| SCI_FED_PACK_AGG_1 | 2 | 8 | 256 | 4 | 0.108465 | 0.034230 | 3.17x | 3.03e-6 |
| SCI_FED_PACK_AGG_2 | 4 | 8 | 512 | 4 | 0.256646 | 0.096363 | 2.66x | 5.46e-6 |
| SCI_FED_PACK_AGG_3 | 4 | 16 | 1,024 | 4 | 0.658983 | 0.200336 | 3.29x | 1.42e-5 |

Size signal:

| Experiment | Baseline Input Bytes | Packed Aggregate Input Bytes | Input Reduction | Baseline Result Bytes | Packed Aggregate Result Bytes | Result Reduction |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| SCI_FED_PACK_AGG_1 | 5,348,686 | 1,337,416 | 4.00x | 3,765,605 | 941,377 | 4.00x |
| SCI_FED_PACK_AGG_2 | 10,698,113 | 2,674,737 | 4.00x | 7,531,490 | 1,882,940 | 4.00x |
| SCI_FED_PACK_AGG_3 | 21,395,219 | 5,347,618 | 4.00x | 15,063,223 | 3,765,707 | 4.00x |

Interpretation:

- This is an SCI extension result, not a KCI manuscript claim.
- Each simulated site keeps its own CKKS context and sends packed encrypted
  genotype blocks. The evaluator computes encrypted block PRS sums, and each
  site decrypts its own aggregate.
- Compared with the per-site individual PRS baseline, packed aggregate
  site-level statistics are about 2.66-3.29x faster in evaluation and reduce
  input/result ciphertext bytes by about 4x at block size 4.
- This is now a plausible SCI technical contribution track for site-level
  privacy-preserving genomic risk summaries.

## Real CAD PGS Federated Packed Aggregate

Result file:

```text
results/real_pgs_federated_packed_aggregate.csv
```

Key result:

| Matched SNPs | Sites | Samples | Block Size | Slots Used | Baseline Eval Seconds | Packed Aggregate Eval Seconds | Aggregate Speedup | Aggregate Site-Sum Max Error |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 512 | 4 | 32 | 4 | 2,048 | 0.307100 | 0.090971 | 3.38x | 3.02e-6 |
| 1,024 | 4 | 32 | 4 | 4,096 | 0.317537 | 0.099505 | 3.19x | 3.20e-6 |
| 2,048 | 4 | 32 | 2 | 4,096 | 0.345695 | 0.198642 | 1.74x | 9.69e-6 |
| 4,096 | 4 | 32 | 1 | 4,096 | 0.389816 | 0.393364 | 0.99x | 1.39e-5 |

Size signal:

| Matched SNPs | Baseline Input Bytes | Packed Aggregate Input Bytes | Input Reduction | Baseline Result Bytes | Packed Aggregate Result Bytes | Result Reduction |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 512 | 10,697,890 | 2,674,460 | 4.00x | 7,531,484 | 1,883,011 | 4.00x |
| 1,024 | 10,698,668 | 2,674,159 | 4.00x | 7,531,900 | 1,883,012 | 4.00x |
| 2,048 | 10,697,715 | 5,349,453 | 2.00x | 7,531,711 | 3,765,939 | 2.00x |
| 4,096 | 10,696,329 | 10,698,732 | 1.00x | 7,531,296 | 7,531,014 | 1.00x |

Interpretation:

- This moves the SCI aggregate protocol from synthetic data to public
  1000 Genomes chr22 dosages matched with the CAD PGS Catalog score PGS004941.
- The speed and ciphertext-size benefit tracks the available CKKS slot budget:
  block size 4 at 512/1,024 SNPs gives about 3.19-3.38x faster aggregate
  evaluation, block size 2 at 2,048 SNPs gives 1.74x, and block size 1 at 4,096
  SNPs removes the packing advantage.
- The result motivates SNP blocking or chromosome/window partitioning as the
  next SCI design step for larger PGS panels.

## Real CAD PGS Windowed Federated Packed Aggregate

Result file:

```text
results/real_pgs_windowed_federated_packed_aggregate.csv
```

Key result:

| Matched SNPs | Direct Block | Windowed Block | Windows/Site | Baseline Eval Seconds | Direct Aggregate Eval Seconds | Windowed Aggregate Eval Seconds | Windowed Speedup vs Baseline | Windowed Speedup vs Direct | Windowed Site-Sum Max Error |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 512 | 4 | 4 | 1 | 0.298279 | 0.093640 | 0.100391 | 2.97x | 0.93x | 5.32e-6 |
| 1,024 | 4 | 4 | 1 | 0.372227 | 0.110889 | 0.091799 | 4.05x | 1.21x | 5.53e-6 |
| 2,048 | 2 | 4 | 2 | 0.407743 | 0.220941 | 0.221812 | 1.84x | 1.00x | 1.03e-5 |
| 4,096 | 1 | 4 | 4 | 0.445004 | 0.442281 | 0.435843 | 1.02x | 1.01x | 1.32e-5 |

Interpretation:

- Windowing keeps the packed aggregate block size at 4 by splitting SNPs into
  1,024-SNP windows, even when the full matched-SNP vector would force direct
  packing down to block size 2 or 1.
- On the current TenSEAL implementation, windowing preserves numerical accuracy
  but gives only modest evaluation gains at 2,048-4,096 SNPs because each
  window is evaluated as a separate packed aggregate call.
- This result sharpens the next SCI engineering target: reuse one CKKS context
  across windows and batch/fuse window-level operations before claiming
  large-scale windowed speedup.

## Real CAD PGS Batched Windowed Federated Packed Aggregate

Result file:

```text
results/real_pgs_batched_windowed_federated_packed_aggregate.csv
```

Key result:

| Matched SNPs | Windows/Site | Independent Windowed Eval Seconds | Batched Windowed Eval Seconds | Batched Eval Speedup vs Independent | Batched Compute Speedup vs Independent | Context Reduction | Site-Sum Max Error |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 512 | 1 | 0.093102 | 0.093244 | 1.00x | 2.83x | 3.78x | 5.29e-6 |
| 1,024 | 1 | 0.100839 | 0.099910 | 1.01x | 2.85x | 3.89x | 2.09e-6 |
| 2,048 | 2 | 0.176353 | 0.200479 | 0.88x | 4.02x | 7.91x | 2.93e-6 |
| 4,096 | 4 | 0.388682 | 0.393431 | 0.99x | 5.44x | 15.59x | 2.16e-5 |

Interpretation:

- The batched implementation encrypts all site-window aggregate blocks under
  one reusable CKKS context instead of creating one context per site-window call.
- Pure evaluation time remains similar because the number of encrypted dot
  products is essentially unchanged, but compute time including context,
  encryption, evaluation, and decryption improves by 2.83-5.44x.
- Public context material falls by about the number of formerly independent
  site-window calls, from about 4x at one window per site to about 16x at four
  windows per site.
- This completes the first context-reuse step. The next SCI optimization should
  reduce the number of encrypted dot products or add rotation-aware individual
  PRS packing rather than only reusing context.

## Submission Bundle Integrity

Generated files:

```text
submission/bundle_integrity.md
submission/bundle_checksums.sha256
submission/public_data_provenance_audit.md
submission/public_data_provenance_audit.json
submission/ethics_data_use_audit.md
submission/text_encoding_readability_audit.md
submission/reviewer_response_readiness_audit.md
```

Key result:

- ZIP entries and planned-file coverage are reported in `submission/bundle_integrity.md`.
- The current bundle SHA-256 is reported only in `submission/bundle_integrity.md`
  to avoid a self-changing hash inside bundled documentation.
- ZIP coverage: all hashed planned files are present in the ZIP.
- Public-data provenance records local SHA-256 hashes for the 1000 Genomes chr22
  files, PGS Catalog scores/metadata, and SQLite cache while keeping raw public
  inputs outside anonymous review ZIPs.
- Ethics/data-use checks verify public-data scope, raw-input non-redistribution,
  local-only author/portal evidence, and institutional-review caveats.
- Text encoding/readability checks verify UTF-8 Markdown, portal Korean text,
  and DOCX-extracted Korean text before submission packaging.
- Reviewer response readiness checks verify that likely reviewer concerns are
  backed by manuscript text, novelty positioning, supplementary documentation,
  and local audits.

## What This Supports in the KCI Paper

The current code supports three paper claims:

1. Genomic dosage vectors can be encrypted and used for PRS computation without
   exposing raw genotypes to the analysis server.
2. CKKS approximate arithmetic preserves PRS with very small numerical error in
   synthetic, public genotype, and real PGS settings.
3. BFV and BGV fixed-point PRS can be run on the same real PGS inputs, including
   the full allele-compatible chromosome-22 matched subset through chunking.

## Paper Assets

Generated figures:

```text
paper_assets/figures/fig1_encrypted_prs_workflow.svg
paper_assets/figures/fig2_vcf_preprocessing.svg
paper_assets/final_figures/fig3_openfhe_controlled_scheme_comparison_300dpi.png
paper_assets/figures/fig4_real_cad_pgs_validation.svg
paper_assets/final_figures/he_ciphertext_space_map_ckks_bfv_bgv_300dpi.png
paper_assets/final_figures/fig6_chunked_he_scheme_comparison_300dpi.png
```

Generated tables:

```text
paper_assets/tables/generated_tables.md
```

Rendered figure previews:

```text
paper_assets/previews/
```
