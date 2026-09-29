# Related work

## Approximate and exact homomorphic arithmetic

CKKS supports approximate arithmetic over encrypted real or complex vectors
[1]. This property matches PRS weights and imputed dosages, but the decoded
result must be evaluated with an explicit error budget. BFV and BGV instead
operate over plaintext rings. AutoHE-PRS represents real-valued effect weights
with a fixed-point integer scale for these schemes and rejects a plan when the
conservative dot-product bound can wrap modulo the plaintext modulus.

The current implementation uses OpenFHE for controlled CKKS, BFV, and BGV
execution. Earlier repository experiments use TenSEAL [4] for CKKS and BFV.
Records identify the backend and measurement provenance, so a TenSEAL result is
not treated as an OpenFHE measurement.

## Privacy-preserving genomic computation

Prior work has evaluated homomorphic encryption for GWAS statistics [5].
AutoHE-PRS addresses a narrower workload, the weighted sum that defines an
additive polygenic score. Its optimization variables include the HE scheme,
ring dimension, modulus and scale settings, SNP window size, chunk count,
packing layout, reduction schedule, and rotation-key set. The method does not
claim to compile arbitrary genomic analysis programs.

## Public genomic and scoring resources

The 1000 Genomes Project provides the public genotype source used in the
repository [6]. PGS Catalog provides scoring-file metadata and the
`effect_allele` and `effect_weight` fields used during harmonization [2,3].
PGS004941 is a coronary artery disease score linked to the China Kadoorie
Biobank study [7,8].

The frontend matches genome build, chromosome, position, REF/ALT alleles, and
effect-allele direction before encryption. It records duplicates,
multiallelic records, missing genotypes, unsupported weight types, and allele
mismatches. Strand-ambiguous A/T and C/G records are excluded by default unless
the user explicitly accepts them. This preprocessing is part of the evaluated
system because a numerically accurate HE computation cannot correct a
misoriented allele.

## Autotuning position

AutoHE-PRS separates three decisions:

1. A deterministic validator rejects structurally invalid or conservatively
   unsafe plans.
2. A learned cost model predicts latency, memory, serialized size, numerical
   error, and execution success for plans that pass static validation.
3. The highest-ranked plans are executed, and final selection uses successful
   measured results when they exist.

The learned model does not decide cryptographic security or arithmetic
validity. OpenFHE context generation remains the authoritative runtime check.
When no top-ranked plan completes in OpenFHE, the output is marked
`predicted_only_no_successful_actual_validation`; it is not reported as a
validated optimum.

## References

1. Cheon, J. H., Kim, A., Kim, M., and Song, Y. (2017). Homomorphic Encryption
   for Arithmetic of Approximate Numbers. *ASIACRYPT 2017*, 409-437.
   https://doi.org/10.1007/978-3-319-70694-8_15
2. Lambert, S. A., Gil, L., Jupp, S., et al. (2021). The Polygenic Score
   Catalog as an open database for reproducibility and systematic evaluation.
   *Nature Genetics*, 53, 420-425.
   https://doi.org/10.1038/s41588-021-00783-5
3. PGS Catalog. Download Information. https://www.pgscatalog.org/downloads/
4. Benaissa, A., Retiat, B., Cebere, B., and Belfedhal, A. E. (2021). TenSEAL:
   A Library for Encrypted Tensor Operations Using Homomorphic Encryption.
   arXiv:2104.03152.
5. Sim, J. J., Chan, F. M., Chen, S., Tan, B. H. M., and Aung, K. M. M. (2019).
   Achieving GWAS with Homomorphic Encryption. arXiv:1902.04303.
6. The 1000 Genomes Project Consortium. (2015). A global reference for human
   genetic variation. *Nature*, 526, 68-74.
   https://doi.org/10.1038/nature15393
7. PGS Catalog. PGS004941: Coronary artery disease.
   https://www.pgscatalog.org/score/PGS004941/
8. The China Kadoorie Biobank Collaborative Group. (2024). Joint impact of
   polygenic risk score and lifestyles on early- and late-onset cardiovascular
   diseases. *Nature Human Behaviour*, 8, 1810-1818.
   https://doi.org/10.1038/s41562-024-01923-7
