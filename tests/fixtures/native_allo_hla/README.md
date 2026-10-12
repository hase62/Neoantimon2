# Native allo-HLA synthetic fixtures

All observations and assay identities are synthetic. No external patient data
or external assay/predictor execution is represented.

`complete-observations.json` supplies one observation per native dimension over
the frozen Phase C complete graph: source and presenter expression, scoped
allele loss, individual Class II APM gene, MS with both assignments unresolved,
and missing TCR evidence. The frozen complete count is five.

`wrong-scope-observations.json` supplies wrong-subject expression, tumor LOH in
an unrelated APC compartment and a wrong-responder TCR assay. All are retained
with provenance and not applicable; there is no count or negative-penalty effect.

`focused-observations.json` contains 27 explicit measurement/scope examples and
independently declared expected applicability states. Remaining acceptance cases
and negative schema examples are exercised in `test-native-allo.R`.

`assignment-evidence-cases.json` records explicit resolved-assignment acceptance
and rejection cases: missing evidence, arbitrary paths/strings, legacy escape
fields, missing semantic fields, unsupported/failed results, independent source
and presenter evidence, mismatched claims, frozen references, physical fields,
and the separate supplied exhaustive assay contract. Schema structural checks
and runtime cross-record/graph integrity results are declared separately.
Inline assignment records are semantic; `provenance.assignment_evidence` retains
physical artifact metadata outside the native observation ID.
