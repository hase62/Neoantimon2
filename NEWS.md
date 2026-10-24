# Neoantimon2 0.3.0 (unreleased release candidate)

* Make the README easier to navigate and add four complete, offline examples
  covering multiple patients and separately declared sample VCFs.
* Improve VCF input robustness: harden ALT syntax, duplicate INFO/FORMAT keys
  and Number=A/R cardinality checks; accept valid omitted trailing FORMAT values.
  Preserve existing normalization, genotype/phase scope, IDs and missingness.
* Add `neo_download_resources()`, `neo_check_resources()`,
  `neo_check_predictors()` and `neo_setup_predictor()`. Reuse the shared resource
  catalog and retain the existing Python/shell utilities. Licensed predictors
  are acquired separately; Python environment installation requires explicit consent.
* Add an editable scientific workflow overview showing retained origins,
  Native/External boundaries, and explicit specialized-workflow prerequisites.
* Validate source/installed examples, clean installation, public VCF import,
  resource setup and backward compatibility. Public-data checks assess ingestion
  and supported processing, not peptide or predictor accuracy.
* Keep scientific/API contracts intact. Multi-sample-column VCFs remain
  unsupported. Real licensed predictor execution, a real MHCflurry environment
  and other platforms remain unvalidated. Research outputs do not establish
  clinical immunogenicity, treatment efficacy or clinical risk.

This candidate is local pending owner approval; no v0.3.0 tag or GitHub Release
has been created. The immutable v0.2.0 release and assets are unchanged.

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
