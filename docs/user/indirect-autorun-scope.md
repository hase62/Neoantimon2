# Indirect autorun scope and future development candidates

This is a limited literature-derived workflow, not official PIRCHE software or
an exhaustive product parity audit. No numerical/clinical equivalence is established.

| Area | Current autorun | Future candidate / evidence boundary |
|---|---|---|
| Binding model | 4.3k BA standard with explicit BA1000 assumption; separately selected 3.2 strict legacy | Official v4 uses proprietary Frost and different rank semantics; any authorized integration needs a separately named/validated profile, not copied weights/formulas |
| Presentation | Explicit DRA*01:01–two-field DRB1 | DQ, DP, DRB3/4/5 and additional explicit alpha/beta combinations require new scope validation; broader presentation has been described in PIRCHE literature |
| Direction | Human SOT donor-to-recipient | HSCT/reverse/bidirectional contexts need separately justified self/presenter roles and profiles |
| MHC class | Class II, 15mer ligand / 9mer primary core; 4.3k retains orientation | Class I/PIRCHE-I would need an independent workflow, not reinterpretation of Class II counts |
| HLA/reference uncertainty | Exact typed alleles and literal supplied precursor sequences; no guessed completion | Ambiguous typing, missing sequence completion and alternative resource policies require explicit uncertainty handling and independent evidence |
| Self universe | Supplied scoped HLA precursor set; 4.3k exact core/orientation/presenter, legacy forward core/presenter | Broader protein universes/minor antigens are separate scope choices; not established as current product parity here |
| Official numerical agreement | Not tested | Authorized matched-input comparison is needed; matching a count alone would not validate clinical performance |
| Clinical use | Not evaluated | Independent clinical evidence and appropriate approvals are required; no clinical score/risk claim |
| Execution recovery | Dry-run and checksum-verified completed-run reopening | Interrupted-batch recovery is a software enhancement; not a claim about official PIRCHE capabilities |

These are candidates, not features implemented in this verification task. Existing
foundation/import workflows can represent states beyond this runner's scope;
that does not establish automatic execution support.

Sources:
- [Official PIRCHE TxPredictor handbook, peptide-HLA binding section](https://s3.eu-central-1.amazonaws.com/content.pirche.com/TxPredictor_The_Handbook_web.pdf)
- [Primary technical review of PIRCHE-II, 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC6971131/)
- [Primary v3/v4.2 methods comparison, 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC11850546/)
- [DTU NetMHCIIpan 3.2](https://services.healthtech.dtu.dk/services/NetMHCIIpan-3.2/)

The external sources establish selected reported capabilities, not this package's
implementation or accuracy. Current local behavior is defined by code/tests and
the verified profile. No new predictor, clinical study or scientific algorithm
was implemented for this checklist.
