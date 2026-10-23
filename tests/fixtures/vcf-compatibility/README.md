# VCF component regression fixtures

These hand-authored files are synthetic, test-scale and non-clinical. They are
not public patient data, human reference sequences, or peptide truth sets.

- `reference.fa` and its exact FAI contain one short `chrT` sequence. The initial
  `CAAA` tract makes insertion/deletion left alignment independently checkable.
- `short-format.vcf` is valid: VCF permits omitted trailing FORMAT fields. The
  declared `AD` is missing, not zero; the original `0/1` text must be retained.
- `trailing-alt.vcf` is invalid: a trailing comma denotes an empty ALT token. It
  must not be silently dropped by tokenization.

The remaining small positive/negative records are generated deterministically
in `test-vcf-robustness.R`. Tests call internal VCF/phase components directly and
stop before sequence reconstruction, HLA candidate generation or prediction.
See [VCF compatibility validation](../../../docs/user/vcf-compatibility.md) for
public-data provenance, independent oracles and the opt-in larger suite.
