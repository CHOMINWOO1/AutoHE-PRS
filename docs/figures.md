# Draft Figures

Rendered SVG files are available in:

```text
paper_assets/figures/
```

Generated files:

- `fig1_encrypted_prs_workflow.svg`
- `fig2_vcf_preprocessing.svg`
- `fig3_synthetic_ckks_scaling.svg`
- `fig4_real_cad_pgs_validation.svg`

## Figure 1. Privacy-Preserving PRS Workflow

```mermaid
flowchart LR
    A["Data owner\nGenotype dosage 0/1/2"] --> B["CKKS encryption"]
    B --> C["Analysis server\nNo secret key"]
    C --> D["Encrypted dot product\nEnc(genotype) · weights"]
    D --> E["Encrypted PRS result"]
    E --> F["Data owner\nDecryption"]
    F --> G["PRS score"]
```

Caption draft:

> Overview of the proposed privacy-preserving PRS computation workflow. The
> analysis server receives encrypted genotype dosage vectors and computes PRS
> without access to the secret key or raw genotype values.

## Figure 2. 1000 Genomes VCF Preprocessing

```mermaid
flowchart LR
    A["1000 Genomes chr22 VCF"] --> B["Biallelic SNP filter"]
    B --> C["GT field parsing"]
    C --> D["Alternate allele dosage\n0/1/2"]
    D --> E["Minor allele count filter"]
    E --> F["Dosage matrix\nsamples x SNPs"]
    F --> G["Plaintext and CKKS PRS experiments"]
```

Caption draft:

> Preprocessing pipeline for public-data validation. Biallelic SNPs are selected
> from the chr22 VCF, genotype fields are converted to alternate allele dosage,
> and monomorphic variants in the selected cohort are filtered before PRS
> computation.

## Figure 3. KCI-to-SCI Extension Path

```mermaid
flowchart TD
    A["KCI pilot\nsample-wise CKKS PRS"] --> B["Real PGS Catalog weights"]
    B --> C["Packed CKKS optimization"]
    C --> D["Federated multi-institution protocol"]
    D --> E["GWAS-lite and covariate-adjusted statistics"]
    E --> F["SCI framework paper"]
```

Caption draft:

> Planned extension path from the KCI pilot implementation to a scalable
> SCI-level privacy-preserving genomic analysis framework.
