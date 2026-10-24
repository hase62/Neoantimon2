# External prediction interoperability

The [workflow overview](workflow-overview.md) shows the Native and independent
External routes and their evidence boundaries.

Neoantimon2 supports three workflows. External prediction is optional; native
analysis never needs an Internet connection, registry or new executable.

External import and its downstream views currently require the original
completed `neo_project`. `annotations.json` and predictor output X alone are
not a supported standalone import. Preserve original inputs and the rerun
recipe outside disposable runs so that the completed project can be recreated.

| Evidence or downstream view | External evidence support |
| --- | --- |
| Peptide-intrinsic and peptide/HLA predictions | Exact sequence/full-presenter matching; all candidate-origin associations retained |
| `neo_external_results()` / `neo_external_explain()` | Explicit external analysis handle or fingerprint, with the original completed project |
| User-selected external ranking | Declared channel/direction; ranks within sample/class; raw scores receive no fabricated percentile or grade |
| Opt-in external somatic grading | Declared MHC percentiles and original eligible core somatic/technical prerequisites |
| Default `neo_run()`, `neo_results()` / `neo_explain()` | Existing native evidence and semantics |
| Vaccine, longitudinal, receptor/TCR, presentation and allo overlays | Existing native evidence only |
| Sample/tumor-context/TCR-specific external scopes | Unsupported in external-prediction-1.0.0 |

## Native analysis

Keep existing `neo_project()` / `neo_run()` calls and predictor configuration.
Configured NetMHCpan, MHCflurry and NetMHCIIpan adapters work as before.
`neo_results(p)` and `neo_explain(p, candidate_id = ...)` keep native semantics.
Every completed run also creates a verified bundle under `results/export_bundle`.
Runs without predictors or somatic evidence still export their candidates.

```r
p <- neo_run(p)
bundle <- neo_export_bundle(p) # verifies and reopens the automatic export
copy <- neo_export_bundle(p, dir = "/selected/new/export-directory")
```

The bundle has distinct **literal** peptide queries and distinct full
peptide/MHC queries, alongside every candidate, every origin and every
candidate-origin association. Array/object TSV cells are JSON. Complete
persisted records and peptide provenance occurrences are in annotations.json,
including source/long-read observations, events, hypotheses and origin paths.
Extension event TSVs retain their original source_event_id and
long_read_source_event_id namespaces. query_coverage.tsv labels legacy or
unresolved presenters UNSUPPORTED_PRESENTER; those candidates still appear in
the full annotation and peptide query list, without fabricated HLA queries.
Candidate annotation is never deduplicated by sequence. The manifest pins the
completed source run and SHA256 of every file; changes are rejected on reuse.
An empty run has valid header-only tables.

## Direct standardized input

Run any predictor outside Neoantimon2. Normalize its output to the separately
versioned external-prediction-1.0.0 TSV. The algorithm does not need to be in
the native predictor registry, and the output needs no candidate/origin IDs.
For a full sample schema, see
[`inst/extdata/external-standard-example.tsv`](../../inst/extdata/external-standard-example.tsv).

```r
a <- neo_analyze_external(p, "/preserved/input/standard.tsv",
  policy = list(mode = "external_rank", provider = "synthetic-model",
    run_id = "demo-run-1", task = "external:synthetic_score", metric = "score",
    score_direction = "higher_is_better"))
results <- neo_external_results(a)
explanation <- neo_external_explain(a, results$candidate_id[[1]])
reopened <- neo_external_results(a$analysis_id, project = neo_project(path = p$path))
stopifnot(identical(results, reopened))
```

Input columns are exactly:

| Column | Meaning |
| --- | --- |
| schema_version | `1.0.0`, including in every row |
| peptide_sequence | Literal, nonempty uppercase AA20+U; no substitutions, case folding or trimming |
| scope | `peptide` (sequence-intrinsic) or `peptide_hla` (sequence and presenter intrinsic) |
| mhc_key | Empty for peptide; full canonical key for peptide_hla |
| task | Native task name or a new `external:...` namespace; no native-task aliasing |
| provider | Predictor identifier |
| predictor_version, model_version | Declared versions; empty when unknown |
| run_id | Provider run/source identifier, separate from the Neoantimon2 run |
| metric | Named metric within the task/provider/model/run channel |
| metric_kind | `raw_score` or truly declared `percentile_rank` |
| metric_semantics | `raw`, `mhc_binding_percentile`, `mhc_presentation_percentile` or `external_percentile` |
| score | Finite numeric when evaluated; empty otherwise |
| score_direction | `higher_is_better` or `lower_is_better` |
| unit | Explicit units when meaningful; percent ranks require `percent` |
| status | `evaluated`, `unsupported`, `not_evaluated` or `failed` |
| reason | Required for every non-evaluated row |
| source_row | Original source row/reference for audit, independent of output order |

Quote complete TSV cells with double quotes and double internal quotes (R:
`write.table(..., sep = "\t", qmethod = "double")`). Do not put tabs/newlines
inside values. Missing evaluated scores, NA/NaN/Inf, malformed values,
unknown directions and duplicate/conflicting logical keys fail validation.
Percent ranks require lower-is-better, percent units and values in [0,100].
Metric meaning must be consistent within each provider/model/run/task/metric.

Matching is exact. Class II needs the full alpha~beta key, for example
`HLA-DQA1*05:01~HLA-DQB1*02:01`; a beta-only token is rejected. A converter
must explicitly map reversible tool allele spellings into canonical keys;
Neoantimon2 never guesses missing chains. MHC tasks need peptide_hla scope.
Known class-specific tasks attach only to applicable HLA classes. TCR-specific,
sample-specific and tumor-context-specific scopes are unsupported in v1: do
not declare such scores peptide-intrinsic merely to make them importable.
For namespaced tasks the author must truthfully declare the applicable scope.

Queries may be sorted, reordered and deduplicated externally. One prediction
can attach to multiple samples and origins only because the declared scope is
sequence/presenter intrinsic. Fan-out copies share an evidence ID and represent
one measurement. A full candidate-anchored result retains NOT_PREDICTED rows;
absence is never negative evidence. Non-evaluated statuses keep empty scores.
Unmatched queries go to unmatched.tsv and the import report. Original annotation
IDs and all links remain visible.

## Interpretation policies

`evidence_only` retains original native grade/rank and attaches every imported
channel. It requires only `list(mode = "evidence_only")`.

`external_rank` requires explicit provider/run/task/metric/direction. Select
predictor_version and model_version too when multiple versions coexist.
Evaluated scores order within each sample and HLA class, with lexical candidate
ID ties; missing scores remain visible and unranked. This is labelled a
**user-selected external rank**. Native final_grade stays unchanged. There is
no averaging, scale harmonization or clinical probability.

`somatic_grade` is opt-in `external-percentile-grade-v1`. It requires a selected
Class I or II binding/presentation task, peptide_hla scope, truly declared MHC
percent-rank semantics and original core somatic/technical evidence:

```r
a <- neo_analyze_external(p, "/preserved/input/declared-percent-ranks.tsv",
  policy = list(mode = "somatic_grade", provider = "my-provider",
    run_id = "prediction-run-1", task = "mhc1_presentation", metric = "my-percent-rank",
    score_direction = "lower_is_better", hla_class = "I"))
```

This versioned selection policy reuses the existing somatic grading constructor,
original thresholds, technical eligibility, origin selection and RNA rules.
It does not change `somatic-priority-v1` or its native primary constraints.
An arbitrary score cannot enter A-E grading. Derived grades appear in
external_final_grade and external_grade_status, with the full selected-predictor
trace in external_priority. Native final_grade, grading_status, primary fields
and native rank remain intact. Candidates outside the selected core class/prerequisites
are NOT_APPLICABLE with no external grade. For candidates shared
with extensions, grading uses core somatic origin paths explicitly labelled
CORE_SOMATIC_ORIGINS, while explanation retains the complete origin graph.
A-E are evidence tiers, not calibrated clinical outcomes.

Each selected analysis persists separately in
`results/external_predictions/<fingerprint>/`, with normalized input, complete
graph/evidence JSON, results.tsv, matched_evidence.tsv, unmatched.tsv,
import_report.json, report.md and checksummed manifest. Multiple models coexist.
Use the returned handle or exact fingerprint to reopen it; it never becomes
the default native view. Vaccine, longitudinal, receptor, presentation and allo
subsystems continue to use native results and are not scientifically integrated
with this new external analysis view.

## Local and registry converters

Provide an R or Python script obeying:

```text
Rscript converter.R --input X --output standard.tsv --bundle export
python3 converter.py --input X --output standard.tsv --bundle export
```

The script reads X and writes the standard table; it does not run the predictor,
reconstruct origins or rank candidates. The bundle argument may be unused.
No private annotation should be uploaded to any service. Predictors need not be
installed to import their precomputed output. A new local converter works offline:

```r
a <- neo_analyze_external(p, "/preserved/input/X.tsv",
  converter = "/preserved/converters/my_converter.R", timeout = 60,
  policy = list(mode = "evidence_only"))
```

Converter dependencies are user-installed. Execution is arbitrary local code,
not a security sandbox. Neoantimon2 copies input/script into a disposable run
workspace, uses a finite timeout, validates output and captures interpreter
version/hash, script and input SHA256, exit status, stdout/stderr and command
placeholders. Logs stay local and may contain whatever the converter prints.
Failures preserve diagnostics under cache/external_conversions/failed-*.

Community lookup is optional, with an explicitly configured local index or HTTPS
URL. The separate public [community registry](https://www.hase62.jp/Neoantimon2-Converters/)
uses the canonical HTTPS index
`https://www.hase62.jp/Neoantimon2-Converters/index.json`. Use that URL: the
github.io alias was observed redirecting through HTTP, and the R downloader
does not independently restrict redirect hops. The canonical origin was
HTTPS/byte-identity and completed-project workflow checked; arbitrary
URLs are not thereby trusted.
See its [submission/execution guide](https://github.com/hase62/Neoantimon2-Converters/blob/main/README.md)
for proposals, GitHub OK/NG reactions and license boundaries. Votes do not
certify science or authorize execution. The registry stores issue links and vote snapshots in
feedback.json, separately from immutable converter files. Rebuilding index.json
preserves those values. The discovery page labels counts as snapshots with a
recorded date or unknown date; unlinked converters have unavailable feedback.
Reaction refresh remains manual. Prior maintainer tests verified hosted GitHub PR checks,
immutability rejection, Pages and HTTPS import using maintainer synthetic tests;
this is not independent community or licensed-predictor validation. Existing
approved run pins retain their original discovery snapshot and exact script.
The registry guide covers PR submission, exact-ID/script lookup, checksum review,
OK/NG comments and reporting incompatible or broken converters. Submission
requires metadata/documentation and data/safety checks, not scientific peer
review or official algorithm approval. CI never runs contributed scripts.

```r
args <- list(project = p, predictions = "/preserved/input/X.tsv",
  converter_id = "synthetic-example@1.0.0",
  registry = "https://www.hase62.jp/Neoantimon2-Converters/index.json",
  policy = list(mode = "evidence_only"))
do.call(neo_analyze_external, args) # first use refuses and shows pinned bytes to review
a <- do.call(neo_analyze_external, c(args, list(trust = TRUE))) # explicit consent
```

Approved pinned bytes can be reused without network/index access. Native or
direct imports never consult the registry. A deleted run loses the disposable
cache and requires trust again. Preserve raw results, converter bytes, original
inputs and the recipe outside runs to recreate the analysis.

## Executable offline example

Source builds containing `examples/external-predictions.R` also provide a
complete installed-package example; no checkout, licensed predictor or network
is needed after installation:

```r
library(Neoantimon2)
example_file <- system.file("examples", "external-predictions.R",
  package = "Neoantimon2")
stopifnot(nzchar(example_file))
source(example_file, local = TRUE)
external_results
external_explanation
external_example
```

This example is included in v0.3.0, but not in immutable v0.2.0 or the
historical Phase 3 setup commit. It builds on the bundled native Quick Start,
imports an evaluated zero, an explicit unsupported result
and a peptide-intrinsic score, and leaves two candidates without predictions.
A query with a different HLA remains unmatched. Candidate origins, native
results/explanations and original artifact hashes are retained; reopening the
external analysis returns identical results. Scores are synthetic technical
inputs, not prediction-accuracy evidence.

From the repository, run:

```sh
Rscript scripts/development/external_prediction_quickstart.R \
  --data-root ~/codex/data/Neoantimon2 \
  --run-root ~/codex/runs/Neoantimon2
```

The script preserves synthetic source/long-read inputs, X and the local converter
under data/external-interop-synthetic-v1, refuses to overwrite differing input,
and creates a new disposable run. It checks one shared candidate with three
origins, higher-is-better ranking with NOT_GRADED, direct normalized import,
reopening and unchanged native hashes. Small reference fixtures remain authored
repository inputs. Tests can select temporary data/run roots explicitly.

In-memory indices preserve exact matching, complete origin graphs and result
order for large inputs. Per-candidate explanations still include complete
global events/provenance; requesting many candidates can produce large JSON.

## Release scope

Both the immutable v0.2.0 release and v0.3.0 require a completed
local `neo_project`; an export bundle alone
cannot execute downstream analysis. External ranking/explanation and eligible
somatic grading are supported under their declared prerequisites. Other downstream
modules retain native evidence. Synthetic workflow tests do not certify licensed
predictor accuracy, clinical use or execution on other platforms. Community votes
are user reports, not official scientific validation. R/Python scripts run with
your permissions, not in a security sandbox.
