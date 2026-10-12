# T002 independent synthetic truth

All eight snapshots are **synthetic**, not mock tool output and not real data.
No HCC1395, legacy biological example, downloaded sequence or model output was
used. `synthetic-*` checksums/locators are declared identifiers in this contract
fixture, not checksums of a real reference or evidence of a predictor run.
Assignments are manually stipulated possible registers, not predictions.

| Fixture | Hand-verifiable expected relationship |
| --- | --- |
| simple_somatic | One event, one origin, one ligand/candidate; unspecified MHC, no prediction |
| multiple_origins | Two origins of ACDEFGHIK; first has DNA=5, TPM=1; second DNA=1, RNA=9, TPM=20. Projected evidence belongs to the first origin only |
| class_ii_shared_core | AAACDEFGHIKLMMM and TTTACDEFGHIKLNNN are distinct ligands, with ACDEFGHIK at positions 3–11 and 4–12; same presenter, distinct candidates, same epitope |
| multiple_registers | AAACDEFGHIKLMMM has ACDEFGHIK at 3–11 and CDEFGHIKL at 4–12; one candidate, two assignments/epitopes, no selected primary |
| subject_roles | Donor source; recipient self and presenter; one synthetic HLA source, no allo grading |
| linked_antigens | Two donor protein origins for one ligand/candidate; B evidence links only the first source antigen |
| zero_vs_missing | Measured RNA zero versus missing normal count; measured expression zero versus NOT_RUN/null |

`multi_context.json` holds the same sample/molecular candidate in two somatic contexts, two independent summaries, and an additional allo context.
`identity_vectors.json` pins revised candidate, molecular candidate and epitope IDs.

`expected.json` fixes human-written cardinalities and provenance.
`canonical_vectors.json` fixes null/empty/logical/integer/binary64/Unicode/object
serialization and SHA-256 vectors. Snapshot IDs are computed by the independent
Python standard-library oracle in `author_fixtures.py`; it never imports or calls
Neoantimon2. The contract is consulted only to fill required null/empty slots.
Sequences, coordinates, subjects, relationships, statuses, identity argument
orders and expected cardinalities are independently specified in that script.
Production R tests read the committed files; they never regenerate expected values.
Python is an optional authoring aid, not a package dependency.

Review procedure: inspect sequences and coordinates above, inspect origin/evidence
associations, compare explicit identity argument order with the contract, and
review any JSON diff before accepting a fixture update. An R failure is never
resolved by exporting the current R output as expected truth. A separate human
review has not been performed by this implementation task.

A future **mock** fixture must label simulated adapter execution as mock; a future
**real** dataset must retain public source/release/checksum and stay outside this
synthetic fixture set. T001 public data are not T002 truth.
