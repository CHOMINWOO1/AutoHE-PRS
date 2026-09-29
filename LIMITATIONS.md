# Limitations

## Cryptographic scope

- The threat model is a semi-honest evaluator. The implementation does not
  defend against a malicious server, chosen-ciphertext attacks, or collusion.
- The data owner retains the secret key. PGS weights, variant order, matched
  count, sample count, score identifier, execution graph, parameters, and
  ciphertext count are public in the MVP.
- Static modulus bounds are conservative early-rejection checks. OpenFHE
  context generation is the authoritative security and parameter-compatibility
  check.
- Fresh AutoHE-PRS runs serialize the context, public key, secret key,
  multiplication key, and rotation keys separately. The 63 reused historical
  rows expose only a context-plus-public-key proxy. Their key-size target is
  therefore labeled `legacy_context+public_key_proxy` and is not compared with
  fresh total-key measurements.

## Genomic scope

- The MVP supports additive biallelic variants. Dominant, recessive,
  interaction, haplotype, and multiallelic terms are excluded.
- A/T and C/G variants are strand-ambiguous and are excluded by default.
  Frequency-assisted strand resolution is not implemented.
- PGS weights are not private.
- Missing genotypes require an explicit policy. `zero` and `mean_dosage` can
  change the scientific estimand and are reported in the harmonization record.
- The current BFV/BGV path accepts hard-call dosages in `{0,1,2}`. Fixed-point
  imputed-dosage encoding is deferred.
- Sample packing and multi-PGS reuse are represented in the plan schema but are
  not implemented in the deterministic OpenFHE runner.
- The implementation computes PRS. It does not assess clinical validity,
  ancestry calibration, or portability across populations.

## Current real-data discrepancy

The historical pipeline reports 51,227 allele-compatible PGS004941 chr22
variants. The new conservative frontend scores 51,225. It excludes two PGS
positions because the VCF contains duplicate records at those positions. It
also excludes 164 multiallelic VCF records. The difference is recorded in
`data/autohe_smoke/pgs004941_chr22/harmonization_report.json`; it is not
silently forced to match the historical count.

## Benchmark and model scope

- The cost-model smoke dataset contains 63 successful historical OpenFHE
  records plus 18 current local `OpenFHEUnavailable` repeats. It is not
  sufficient evidence for hardware-transfer claims or a final autotuner
  comparison.
- The group-held-out smoke evaluation trains on public 1000 Genomes rows and
  evaluates on synthetic rows. Evaluation-latency MAPE is 13.78% and Spearman
  correlation is 0.980 for this split. These values are an implementation
  check, not a generalization guarantee.
- Relative percentage error is unstable when numerical-error targets are near
  zero. Numerical-error models must be judged with MAE/RMSE and rank metrics in
  addition to MAPE.
- On the same small group-held-out split, the tried XGBoost and
  HistGradientBoosting models had approximately 54% evaluation-latency MAPE,
  so the ridge model was retained. This is a negative result, not evidence
  that linear models are generally preferable.
- The one-PGS-to-another and unseen-window evaluations are explicitly
  unavailable because the measured dataset has only one identified PGS and
  one window group. No metric is synthesized for either split.
- Bayesian refinement uses a small discrete Gaussian-process surrogate. It has
  not yet been benchmarked against grid, random, and Bayesian baselines on the
  full candidate space.

## Runtime verification

OpenFHE-Python is not installed in the local Windows Python environment. Local
benchmark attempts therefore remain explicit `unavailable` records. All final
FHE measurements were made in `genome-he-openfhe:1.5.1`.

The narrow 512-window smoke search initially predicted 13.54 seconds and
measured 86.09 seconds. Expanding the configured search space exposed the
much faster `add_then_reduce` 4,096-window plans. The final selected plan was
predicted at 2.3486 seconds and measured at 0.3586 seconds, a 554.9% absolute
percentage error relative to the measurement. Thus the model helped place the
fast plan family in its top three but did not meet the 20% calibration target.
Claims about latency prediction transfer must remain limited.

## Deployment scope

This repository is a research implementation. It is not a clinical system,
does not implement authentication or access control, and does not protect
sample count, variant count, timing, or ciphertext-size side channels.
