# Expanded indirect class II validation

This is technical interoperability evidence with genuine official IPD-IMGT/HLA
proteins and authentic pretrained NetMHCIIpan 4.3k, **fictional donor/recipient
typing**, not patient observations or clinical performance.

| Case | Group | Donor Origins | Self Origins | Unique queries | Donor pMHCs | Candidate-Origin links | Counted units | Inverted query rows |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| DQ_DP | DP | 351 | 351 | 1188 | 702 | 702 | 31 | 98 |
| DQ_DP | DQ | 351 | 351 | 1188 | 702 | 702 | 6 | 0 |
| DRB1 | DR | 351 | 351 | 594 | 351 | 351 | 44 | 0 |
| DRB3_all13 | DR | 3271 | 3271 | 3198 | 2994 | 3271 | 3 | 0 |
| DRB4 | DR | 603 | 603 | 846 | 603 | 603 | 21 | 0 |
| DRB5 | DR | 603 | 603 | 846 | 603 | 603 | 23 | 0 |

The DRB3_all13 case configures all 13 source loci; DRB4/5 are explicitly absent
there and independently generated from full proteins in their own cases. All
13 source loci were independently verified across cases. Counts/rows repeated
across cases are not independent antigens. DQ and DP totals are separate research
molecular counts, never a pooled clinical score. DP inversion rows include donor
and self queries; each maps to the reported reversed-ligand register.

Final source tests: 12280 assertions PASS, 0 skips;
installed focused tests: 421 PASS; full installed suite:
12274 PASS, 1 expected repository-only skip, 0 failures/warnings.
R CMD check --no-manual: 0 errors / 0 warnings / 0 notes, full tests enabled.
A metadata-only utils import declaration corrected a check NOTE after the full source run; all R source stayed identical, and the final full installed suite checked the final namespace. No TeX is installed here; PDF manual compilation was not tested. New focused
coverage has 92 passing assertions. Synthetic tests are kept separate from the
actual model execution evidence.

Independent Python checks verify official source hashes/full-CDS completeness,
all full precursor windows, exact allele identity, raw BA/EL strings and numeric
projection, primary/inverted core coordinates, self subtraction and all links.
Fresh installed R processes reopen all five cases without inference; files and
mtimes stay unchanged. The previous 6836 DRB1 saved result and the 3.2 complete
baseline/public result/explanation are strictly identical. All prior R modules
remain byte-identical to the starting revision. Protein-number guard extensions
do not change any of the earlier molecular outputs.

Authentic negative probes: exact partial allele A*01:01:04
remains incomplete in an offline plan; DRB1*01:142 is
absent from the executed model list and remains unsupported with no fabricated
scores. Sequence-unavailable and unknown-typing negative cases are independently
authored software fixtures, not fabricated biological observations.

Source release 3.66.0 is pinned at `b2055ea33ef464781b9e8ecab20b275b6704f878`. Scripts:
[real local workflow](../../scripts/development/indirect_class2_ipd_validation.R),
[independent read-only verifier](../../scripts/development/verify_indirect_class2_ipd.py).
[Usage and rights](../user/indirect-class2.md). Machine-readable receipts and raw
stdout SHA-256 hashes: [evidence](indirect-class2-validation.json). Full processed
IPD/raw peptide outputs stay local under upstream non-derivative redistribution
terms; they are not bundled. No predictor/weights or clinical data are published.
The exact source is the commit containing this evidence; runtime source archive
hashes are retained in the local review record. No official PIRCHE parity claimed.
