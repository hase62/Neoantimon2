# Supplied peptide candidates, reference percentiles and private review

These development features use the ordinary completed project and existing
external/reporting APIs. The immutable v0.3.0 and reproducibility-v1 snapshots do
not include them. There are still 22 exported APIs; `neo_analyze_external()` has
one new optional argument at the end. The current manuscript and historical
scientific results are outside this feature-validation change.

## Full peptide FASTA input

```r
cfg <- neo_config(paths = list(run_root = "/explicit/output/runs"),
  peptides = list(class_i_lengths = 9L, class_ii_lengths = 15L))
p <- neo_run(neo_project(
  samples = data.frame(sample_id = "sample", patient_id = "subject"),
  hla = list(subject = "HLA-A*02:01"), config = cfg,
  inputs = list(list(input_id = "supplied", input_role = "peptide_fasta",
    sample_id = "sample", path = "/explicit/input/peptides.fasta"))))
```

Plain and gzip files, multiline sequences and LF/CRLF are accepted. The first
whitespace-delimited header token is the record ID; duplicate IDs within a file
are errors with the current and first header lines. Sequences are joined and
uppercased, using the existing AA20+U validation. Empty files/records, sequence
before a header, whitespace inside sequence lines, X, stop characters and other
invalid residues are rejected with file, line and reason. Blank sequence lines
are rejected. Nothing is translated, truncated or split into shorter peptides.

Each record has its own `sequence` Origin with `source_domain=provided_peptide`;
its full header, input ID, actual file SHA256, record ID/order/header line are
retained in the saved peptide artifact's private supplied-peptide provenance.
`source_annotation_ref` identifies that artifact record, and private HTML Origin
details include it. Shared HTML withholds full headers. Candidate keys stay the
existing sequence/HLA/context/sample keys, so duplicate records create one
candidate/query with multiple Origins; different samples remain separate.

FASTA-only runs need no fabricated VCF, coordinates, transcript or reference
bundle. Existing reference validation still applies to declared resources and
mixed molecular inputs. Supplied sequence does not establish mutation, phase,
wild type, tumor specificity or DNA support. Its molecular eligibility/Grade
remain unevaluated. A molecular Origin can independently support a candidate
that also has supplied Origins; supplied Origins do not add event evidence.

Configured Class I/II lengths decide candidacy. Unsupported lengths retain their
whole sequence/Origin and an explicit saved status, even if no HLA candidate is
created. `neo_report_bundle()` includes `supplied_peptide_status` so the original
record, sequence and no-candidate reason can be inspected without a candidate ID.
Predictor-specific residue/length limits remain separate adapter
unsupported statuses. Content SHA256 binding/checkpoints prevent resuming stale
FASTA candidates after the input changes. Preserve original input files.

## Optional local reference percentile

```r
a <- neo_analyze_external(p, "/explicit/input/standard.tsv",
  policy = list(mode = "evidence_only"),
  reference_percentiles = "/explicit/reference/manifest.json")
neo_external_explain(a, neo_results(p)$candidate_id[[1]])
report <- neo_report_bundle(p, analyses = list(a))
report$tables$external_reference_percentiles
```

The separate `reference_percentiles.tsv`/analysis payload has **one row per input
measurement**, an evidence ID, original source row and score text/numeric value,
all associated candidate IDs, result or reason, reference ID/SHA256/count and
method. Fan-out to samples/candidates does not inflate measurement counts. Original
provider raw scores and native percentiles stay in `pmhc_predictions`; derived
values appear beside them in the viewer/explanation, never replace them or enter
ranking/Grade automatically. Existing policies retain their current semantics.

`reference-percentiles-1.0.0` JSON contains a `distributions` array. Each record
requires `reference_id`, `provider`, `predictor_version`, `model_version`, `task`,
`metric`, `unit`, `score_direction`, canonical `mhc_key`, integer `peptide_length`,
`species="NCBITaxon:9606"`, `hla_class="I"`, `reference_set`, `source`,
`source_version`, `original_sha256`, `scores_file`, `scores_sha256`, `count`,
and `license`. Score files have a single `score` header and finite numeric values,
one per line. Paths are resolved relative to the local manifest. No network is
used by this operation, and Core requires neither Python nor HDF5.

Provider/software/model/task/metric/unit/direction/HLA/length must all match.
Unknown model identities, other tasks/alleles/lengths, nonraw metrics and other
scopes are not silently substituted; missing, unsupported and not-evaluated input
statuses remain explicit. Hash/count/nonfinite/duplicate-key errors refuse the
analysis before its new directory is published. Existing analyses reload without
the reference source files using their saved, checksummed artifacts. NULL leaves
the pre-feature fingerprint unchanged.

Distributions are sorted once and queried by binary left/right insertion bounds.
The method is `strict-better-plus-half-ties-v1`:

`100 * (number strictly better + 0.5 * number equal) / reference count`.

Better is lower or higher according to the raw channel direction; zero is the
best end of the resulting percentile. Duplicate reference values retain their
frequency. Results are stored without rounding; three-decimal comparisons are
only for checking pVACtools, whose default mixes lengths. Reference hashes,
metadata and method are included in the optional analysis fingerprint.

`scripts/data/prepare_reference_percentiles.py` is an explicit **offline** HDF5
preparation utility requiring optional h5py/numpy. It checks the supplied source
SHA256 and extracts exactly one specified length dataset. All provenance/licensing
and a verified model identity must be supplied; it refuses unknown model versions.
It never downloads a moving main or runs during Core execution. A dataset name
or a predictor software version alone is not evidence of matching model weights.

## Review JSON save, reopen and evidence updates

Source the existing `inst/examples/report-html.R` helper. To enable human review:

```r
write_neo_report_html(p, "/explicit/output/private-review.html",
  analyses = list(a), privacy = "private", confirm_private = TRUE)
```

The ordinary default remains sharing mode. In private mode, the list shows
未レビュー / 採用 / 除外 / 保留 and provides a review filter. Candidate details
have a separate human-review panel, comment editor and UTC history. Recording a
decision never edits Native eligibility, Grade, rank or external results.

Use **レビューJSONを保存** after editing. JSON is the durable authority; there is
no localStorage/cookie/database reliance under `file://`. Open a new browser
session and use **レビューJSONを読み込む**. The versioned root identifies the run,
analysis and authoritative source-manifest SHA256. Each history record retains
the candidate ID, state/comment/time, normalized evidence SHA256, exact canonical
evidence JSON and displayed external analysis IDs. Record/identifier collections
and object keys are normalized before hashing; presentation sorting/filtering
does not change evidence. Projection-dependent prediction row IDs are excluded;
biological/evidence IDs and all scientific values remain. Numeric positional
arrays retain their order.

A newly generated HTML selecting additional analyses can show changed evidence.
Prior decisions/history remain visible with **根拠更新・再確認が必要**; prior and
current evidence snapshots/IDs are available together. An accepted decision with
changed evidence stays out of confirmed-accepted TSV until the reviewer explicitly
records the current evidence again. The TSV includes candidate ID, peptide, HLA,
decision evidence hash, analysis IDs and UTC time. It omits comments/history and
neutralizes spreadsheet formula prefixes and control separators.

Files for a different project/base run, malformed/unknown versions, invalid
history fields or mismatching stored evidence hashes are refused before review
state changes. Partial HTML retains and resaves validated history for candidates
outside its scope without applying it to other candidates. Comments render as
literal text. Review JSON and snapshots are private; sharing HTML and ordinary
data exports never contain human review history/comments or its original-ID map.
The generic embedded-view download is a scientific view, not a review backup.

## Reproducible software verification

`scripts/development/check_pvac_percentiles.py` compiles the two unmodified methods
from explicitly supplied pVACtools source at commit
`871f2b3d979d5d1988af7a2d3ef2d52e179d1524` and compares them with independent counts
against the freshly checksum-verified HDF5 from data commit
`24fdb2aed1c575d1421eae909a473706edda84dc`. The supplied checksum is enforced.
Length-specific checkpoints are MHCflurry BA 500 nM -> 5.068 and BigMHC EL 0.5 ->
0.617; pVACtools' default mixed-length outputs are 2.017 and 0.267. Core's unrounded
calculation is also checked on those same exported values by the R acceptance
driver. This is a computation comparison, not evidence of model equivalence.

The HDF5 root/group/dataset attributes and pinned README do not identify the
original model weights. Historical HCC1395 model identity is therefore not filled
in. `validate_reference_model.py` instead uses an explicit MHCflurry 2.0.6 model
directory, hashes every model file, evaluates a deterministic fixed 256-peptide
reference set and two queries **in the same loaded model**, and preserves native
output, original Converter output and an explicit model-identity receipt. Its
reference-set percentiles verify software integration; they do not estimate
biological prevalence, immunogenicity, clinical benefit or an established
population calibration.

Reproduction needs explicit local model assets, their upstream terms and the
optional Python packages. Scripts do not assume an existing download or fetch
models automatically. The review ZIP provides exact observed hashes, environment,
commands and inputs; do not substitute a newer model silently. See its report for
actual PASS/FAIL/SKIP results, browser sessions and remaining limitations.
