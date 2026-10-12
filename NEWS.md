# Neoantimon2 0.2.0

Released on GitHub on 2026-10-12: https://github.com/hase62/Neoantimon2/releases/tag/v0.2.0.

* Add four public APIs: `neo_export_bundle()`, `neo_analyze_external()`,
  `neo_external_results()` and `neo_external_explain()`.
* Completed native runs export verified literal peptide/HLA queries and the
  complete origin/annotation graph. External observations match exact sequences
  and full canonical presenters, retaining every origin sharing a peptide.
* Import the strict external-prediction-1.0.0 format; derive explicitly selected
  evidence/ranking and prerequisite-gated compatible MHC percentile somatic
  grading. Raw scores never become invented percentiles or grades.
* Run user-created local R/Python converters, or look up immutable community
  converter IDs with SHA256 verification, explicit first-use execution consent
  and approved offline cache reuse.
* Use the independent Community Converter Registry at
  https://www.hase62.jp/Neoantimon2-Converters/. Community feedback is a manual
  snapshot, not scientific validation, script safety or execution consent.
* Transient external lookup indices preserve output identity and ordering.
  Mutable feedback remains separate from immutable converter definitions.
* Preserve the eight existing native APIs, native predictor workflows, biological
  identities and frozen scientific contracts. Native analysis requires neither
  external predictions nor the Registry. Installed packages retain required NOTICE.
* External analysis requires a completed local neo_project. Independent
  bundle-only downstream execution is unsupported. Vaccine, longitudinal,
  receptor, presentation and allo policies continue to consume native results.
* Validation uses synthetic fixtures; licensed predictors, clinical accuracy and
  other platforms are not certified. Distribution retains PolyForm Noncommercial
  1.0.0 and is not a CRAN submission.

# Neoantimon2 0.1.0

* Initial native package with eight public APIs. The external prediction workflow
  described above is introduced in 0.2.0.
