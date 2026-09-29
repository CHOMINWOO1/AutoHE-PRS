# Threat model

## Parties

- The data owner holds genotype or effect-allele dosage, generates the
  OpenFHE key pair, and retains the secret key.
- The evaluator receives encrypted dosage, public/evaluation keys, and
  plaintext PGS weights. It evaluates `Enc(D) * B` and returns encrypted PRS.
- The result recipient is the data owner in the MVP. It decrypts and decodes
  the PRS and compares it with an authorized plaintext reference during
  validation.

The evaluator is semi-honest. It follows the emitted plan but may inspect all
messages and metadata available to it.

## Confidential values

- Per-sample genotype and effect-allele dosage are encrypted before the
  evaluator receives them.
- PRS results remain encrypted until they return to the secret-key holder.
- The evaluator does not receive the secret key.

## Public values and leakage

The MVP treats these values as public:

- PGS effect weights and PGS identifier;
- matched variant order and count;
- sample and score counts;
- HE scheme, parameters, execution graph, packing and chunk plan;
- ciphertext and evaluation-key counts and serialized sizes;
- operation timing and success/failure status.

The protocol therefore does not hide participation, workload dimensions,
access patterns, or timing/size side channels.

## Assumptions

- OpenFHE's 128-bit classic security profile is used. Static validation is an
  early filter; successful OpenFHE context generation is required.
- The data owner protects its secret key and plaintext input.
- PGS and VCF harmonization occurs in the data-owner boundary before
  encryption.
- The evaluator does not alter the computation, submit malformed ciphertexts,
  or collude with the secret-key holder.

## Out of scope

- malicious evaluators and active protocol attacks;
- chosen-ciphertext security and decryption oracles;
- key compromise, multi-party threshold keys, and collusion;
- private PGS weights;
- private variant or sample counts;
- authentication, authorization, audit logging, and transport security;
- differential privacy for a PRS released in plaintext;
- clinical inference attacks against the final PRS.

These exclusions are protocol limitations, not properties supplied by the
cost model or optimizer.
