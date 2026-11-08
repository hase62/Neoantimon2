# First Real-Data Workflow

Use this page to move from the short Quick Start to a declared somatic VCF/HLA
analysis. The executable training input is **authored synthetic data**, processed
through the ordinary molecular reconstruction route. It is not a patient dataset,
human reference bundle, HLA experiment or real-predictor validation.

[Version choices](../../README.md#choose-a-version) · [Complete practical guide](getting-started.md) ·
[Reference setup](setup.md) · [Result fields](results-and-explain.md)

## 1. Check software and load the example helpers

Install a source revision containing `examples/first-real-run.R`, following the
[installation instructions](../../README.md#installation). For a locally reviewed
checkout, run this in a terminal (replace the example path):

```sh
R CMD INSTALL /path/to/Neoantimon2
```

Restart R after installation. Record the exact source commit and any uncommitted
changes; 0.3.0.9000 alone does not identify them. The immutable release/snapshot
do not contain this tutorial. Installing an older main/release cannot provide
helpers absent from that revision.

```r
library(Neoantimon2)
R.version.string
packageVersion("Neoantimon2")
length(getNamespaceExports("Neoantimon2"))  # 22
source(system.file("examples", "first-real-run.R", package = "Neoantimon2"), local = TRUE)
```

These sourceable example helpers are not additional exported package APIs.
Sourcing defines functions only: no download, provider command or project write.
Their static check uses the existing resource validator, length policy and HLA
resolver. It does not start `neo_run()` or a predictor, including help/version.

## 2. Prepare the reference and a new destination

For the first training run, use the included miniature fixture:

```r
run_root <- tempfile("first-real-run-")  # training only; temporary storage
a <- first_real_run_example(run_root)
```

This prepares declarations in memory, not a run directory. Its nine reference
files and one VCF are read-only package inputs. The GRCh38/GENCODE/MANE labels
are compatibility labels in a synthetic format fixture, not genuine resources.
Never mix its sequence files with real variants.

For real data, replace the **entire** reference configuration with an already
prepared compatible bundle: genome FASTA/FAI; GTF; transcript, coding-transcript
and protein catalogs; Sec metadata; MANE summary/GTF. Follow
[resource preparation](getting-started.md#use-r-based-setup) and keep original
versions/checksums. No download occurs in this tutorial. Optional acquisition is
a separate explicit several-GB operation. Use a durable run root outside source,
reference and model directories; this training temp directory disappears later.

## 3. Declare VCF, sample bindings and HLA

`a$samples`, `a$hla`, `a$inputs` and `a$config` are the ordinary inputs to
`neo_project()`. Replace them together for your own data, for example:

```r
# Real-data template ONLY: replace every placeholder with verified inputs.
# a$samples <- data.frame(sample_id="S1", patient_id="P1")
# a$hla <- list(P1=c("HLA-A*02:01", "HLA-B*07:02"))
# a$inputs <- list(list(input_id="S1-somatic", input_role="somatic_vcf",
#   sample_id="S1", path="/absolute/path/to/S1.vcf"))
# a$config <- your_compatible_reference_and_output_config
```

HLA is keyed by patient, inputs by sample. A sampled VCF must have exactly one
matching sample column. A supported sites-only VCF can instead use its explicit
input binding; that is **not** measured GT, AD or VAF. Do not rename headers,
invent genotype/counts, liftover or repair REF values to suppress an error.
Class II requires an explicit alpha~beta pair, such as
`HLA-DRA*01:01~HLA-DRB1*01:01`. The package retains unresolved HLA states;
this conservative introductory check stops on unresolved typing so that an
empty candidate set cannot be mistaken for a biological negative.

## 4. Choose predictor-free or configured prediction

The training configuration has `predictors=NULL`: generate candidates without
binding predictions. This is a supported mode. Its candidates remain `NOT_RUN`
and `NOT_GRADED`, not measured negatives, score zero or Grade E.

For native prediction, obtain an authorized compatible installation/models
separately, then declare provider/task/executable/options in a **new** config.
See [predictor setup](setup.md#configure-legitimately-obtained-prediction-tools).
Static executable/model-file presence does not establish version, working
dependencies, supported lengths/HLAs or model usability.

Only after reviewing the trusted installation, explicitly opt into capability
probes, for example:

```r
# Separate deliberate action, NOT part of automatic preflight:
# tools <- neo_check_predictors(config=your_predictor_config, smoke_test=FALSE)
# tools$predictors
```

Even `smoke_test=FALSE` executes help/version commands and may probe Python
dependencies. `smoke_test=TRUE` additionally runs synthetic inference. Neither
is silently invoked here. If an explicitly requested executable/model is absent,
fix it or deliberately choose a candidate-only config; nothing is auto-installed
or silently removed. External BA/EL/IM import is a separate
[documented workflow](external-predictions.md), preserving native results.

## 5. State lengths and run the read-only readiness check

The training config explicitly selects **Class I 8:12 and Class II 15**. The
historical package defaults remain **8:11 and 13:25**. The README length example
is an override, not a default change. An omitted length vector in a partial
policy is filled by `neo_config()`; the readiness report displays actual values.

```r
check <- do.call(first_real_run_check, a)
check$software
check$effective
check$checks
stopifnot(check$ok)
```

`FAIL` requires correction; `NOTE` records a limitation/optional absence.
Default output contains field/slot references, counts and controlled status
codes, not sample IDs, patient IDs, paths or raw errors. To inspect private
bindings locally, explicitly use `include_private=TRUE`; never paste that object
into a public issue. The defaults are safe for logging, but project results and
other native errors are not an automatic anonymization service.

The check reads reference files and honors declared checksums (which can read
large files), verifies VCF header bindings, resolves declared HLA syntax,
displays effective lengths, checks static provider files and output permissions.
It writes no directory or probe file. Permissions do not prove quota/free space
or future writes. It does not validate every variant body, REF allele, transcript
context or model. Full cross-resource checking is a separate explicit
`neo_check_resources(config=a$config, validate_reference=TRUE)` operation.

## 6. Execute after resolving FAIL rows

```r
p <- do.call(neo_project, a)  # declares the project; no synthetic predictor backend
p <- neo_run(p)             # ordinary route, on authored synthetic input here
```

This is the point that creates scientific outputs and directories. Do not run
it against your real inputs until reference/typing/provider choices are reviewed.
Keep the script, configuration, input/reference hashes and session information
in protected local storage. A run failure remains a failure; inspect the
identified field and never weaken a validation condition to force completion.

## 7. Read counts, missingness and the existing report

```r
events <- neo_events(p)
origins <- neo_origins(p)
results <- neo_results(p)
c(events=nrow(events), origins=nrow(origins), candidates=nrow(results))
table(results$prediction_status, useNA="ifany")
table(results$grading_status, useNA="ifany")
if (nrow(results)) table(unlist(results$not_evaluated_codes), useNA="ifany")
report <- neo_report_bundle(p)  # in-memory reporting; no second reconstruction
report$manifest$original_project_graph_counts
report$tables$event_summary
report$tables$pmhc_candidates
```

The included run produces **1 Event, 65 Origins, 65 candidates**, all `NOT_RUN`
and `NOT_GRADED`: Class I lengths 8/9/10/11/12 produce 8/9/10/11/12 candidates,
and Class II length 15 produces 15. This is a software example, not 65 validated
binders. No VAF/read counts were measured. Real datasets may produce different
counts, noncoding Events without peptide Origins or unsupported contexts.
Never substitute zero for missing scores/evidence. `not_evaluated_codes` gives
reasons; a failed provider is not simply an unrequested one.

`neo_report_bundle()` already contains Event representatives, unique-peptide,
all-pMHC, prediction and complete Origin/link tables. No new reporting API is
needed. Native status counts do not count separately imported External evidence;
include explicitly selected `analyses` for a full external-evidence report.
For large projects, use a bounded `candidate_id` report when appropriate;
do not automatically create an unbounded interoperability export just to count.

## 8. Trace the Origin rather than trusting the rank alone

```r
if (nrow(results)) {
  explanation <- neo_explain(p, candidate_id=results$candidate_id[[1L]])
  explanation$candidates[[1L]]$origins
  explanation$candidates[[1L]]$events
  neo_explain(p, candidate_id=results$candidate_id[[1L]], format="text")
}
```

One peptide can retain many molecular Origins and HLA alternatives; those are
not extra synthesized peptides or independent biological neoantigens. See
[result reporting](result-reporting.md) and [field interpretation](results-and-explain.md).

## 9. Persist, export deliberately and reopen

```r
run_path <- p$path  # private location; retain it in your protected run record
reopened <- neo_project(path=run_path)
stopifnot(identical(neo_results(reopened), results),
          identical(neo_origins(reopened), origins))
# Optional export to a NEW directory separate from inputs/references/project:
# exported <- neo_report_bundle(reopened, dir="/absolute/new/report-directory")
# exported$path  # manifest.json plus TSVs and complete linkage tables
```

Reopening verifies saved state without rerunning reconstruction. Training inputs
remain packaged/read-only; outputs are in the configured run root. A temp root
is not a permanent archive. Keep real projects and exports durable, protect
identifiers/private paths and review data terms before sharing. Readiness output
is diagnostic evidence, not proof of biological completeness or clinical utility.
