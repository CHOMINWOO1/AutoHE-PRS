# AutoHE-PRS PGS004941 chr22 experiment

Date: 2026-07-28

## Verification scope

The local Python 3.11.9 and Python 3.14 suites each ran 46 tests successfully.
Two OpenFHE-dependent tests were skipped locally. In the pinned OpenFHE image,
both container tests passed: the CKKS/BFV/BGV multi-chunk test and the
generated-script test.

The Docker image digest was
`sha256:08f79bfdef923aeca1f50197aa7901b8a3f33320ca00759d65e010aa0aa9495f`.
All FHE results below are measurements from that image, not simulations.

## Genomic frontend

The streaming frontend read 3,711,629 PGS004941 rows, including 51,391 chr22
rows. It harmonized 51,225 variants to four 1000 Genomes samples: 4,048 were
ALT-oriented and 47,177 were REF-oriented. It excluded 164 multiallelic
records and two duplicate VCF positions.

Preprocessing took 61.6062 seconds. The plaintext reference scores were:

| Sample | PRS |
| --- | ---: |
| HG00096 | -0.017874403152041602 |
| HG00097 | -0.006901380096257772 |
| HG00099 | 0.026931796023706808 |
| HG00100 | 0.0047631991680389145 |

## Model evaluation and transfer result

The retained historical ridge model trained on 18 public-data rows and tested
on 45 synthetic rows. Evaluation-time MAPE was 13.78%, R² was 0.8850, and
Spearman correlation was 0.9798. XGBoost and HistGradientBoosting each had
approximately 54% evaluation-time MAPE on this split.

These held-out results did not translate into calibrated absolute predictions
for the full real workload. For the final top-K group, the model predicted
2.3486 evaluation seconds, while the selected plan measured 0.3586 seconds
(554.9% absolute percentage error relative to the measurement). The model
nevertheless placed all three fastest tested CKKS variants in its top three;
rank 2 was the measured winner. The result supports useful coarse ranking in
this case but rejects the 20% absolute-error target.

The model predicted 19,603,396 ciphertext bytes versus 14,728,728 measured.
Its key-size target is not compared: historical rows contain a
context-plus-public-key proxy, whereas fresh runs serialize all public,
secret, multiplication, and rotation key material.

## Top-K actual validation

All three ranked candidates completed three independent measured runs and
passed decryption, security, memory, and numerical-error constraints.

| Predicted rank | Scaling technique | First modulus | Evaluation median (s) | Total median (s) |
| ---: | --- | ---: | ---: | ---: |
| 1 | FIXEDAUTO | 50 | 0.3720 | 1.0375 |
| 2 | FLEXIBLEAUTO | 50 | 0.3586 | 0.9871 |
| 3 | FIXEDAUTO | 60 | 0.3721 | 1.0392 |

All three plans used CKKS, ring dimension 8,192, a 4,096-variant window,
13 chunks, `add_then_reduce`, a binary reduction tree, and BV key switching.

## Same-input scheme comparison and failures

| Scheme | Valid plan | Evaluation median (s) | Total median (s) | Peak MB | Ciphertext bytes | Total key bytes | Max abs error |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| CKKS | N=8,192; window=4,096 | 0.3586 | 0.9871 | 78.13 | 14,728,728 | 7,217,382 | 3.29e-10 |
| BFV | N=16,384; window=4,096; S=10⁸ | 22.9247 | 24.8219 | 125.77 | 29,408,680 | 21,244,306 | 8.19e-7 |
| BGV | N=16,384; window=4,096; S=10⁸ | 9.6960 | 10.3997 | 161.13 | 43,047,248 | 35,272,941 | 8.19e-7 |

Two failed plan families were retained:

- BFV with ring 8,192, window 512, and scale 10⁶ decrypted successfully but
  violated both numerical-error constraints (`9.16e-5` maximum absolute error
  and `0.0164` maximum relative error).
- BGV with the same small ring was rejected by OpenFHE because ring 8,192 does
  not meet its 128-bit recommendation for the configured large plaintext
  modulus. This observation was added to the deterministic validator.

The scale-10⁸, ring-16,384 recovery plans passed all constraints in 3/3 runs.

## Baseline comparison

| Method | Best measured evaluation time (s) | Regret (s) |
| --- | ---: | ---: |
| Proposed learned top-K | 0.3586 | 0 |
| Manual CKKS baseline | 15.8768 | 15.5182 |
| Uninformed Bayesian, budget 2 | 0.9761 | 0.6175 |
| Random, budget 1 | no valid result | not defined |

The proposed result reduced median evaluation time by 97.74% relative to the
manual configuration. The random candidate was a BFV scale-10⁵ plan and
failed the measured error constraints; it remains in the dataset as a failed
attempt. The small budgets make this an implementation comparison, not a
general claim about search methods.

## Final emitted plan

The final selection is the measured rank-2 CKKS plan. Its authoritative
three-run record is stored in
`results/autohe_smoke/final_selection/authoritative_best_plan.json`. The
cross-platform generated script at
`generated/autohe_smoke/run_best_plan.py` was executed without edits in the
pinned Linux container after generation on Windows. That final smoke execution
completed in 1.0769 seconds end to end, with 0.4065 seconds of homomorphic
evaluation and `1.27e-10` maximum absolute error.

Paper-ready tables, prediction comparisons, and the latency-memory figure are
under `reports/autohe_smoke/final`. Failures, local runtime unavailability,
and scope-incomparable key predictions are retained rather than omitted.
