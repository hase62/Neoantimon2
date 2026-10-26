# Practical user guide

Follow this page from installation to a saved analysis, interpretation and optional
external prediction. The default route below is somatic analysis with separate
single-sample VCFs. Advanced input/analysis modes have their own contracts.

[Install and demo](#install-and-run-the-offline-example) →
[Prepare references](#use-r-based-setup) → [Declare and run](#configure-your-own-project) →
[Read results](#explore-the-completed-run) → [External prediction](#add-external-prediction) →
[Save and reopen](#reopen-and-inspect)

## Install and run the offline example

Use R ≥ 4.1.0. Install the released package, then run the small example:

```r
install.packages("remotes")
remotes::install_github("hase62/Neoantimon2", ref = "v0.3.0", upgrade = "never")
library(Neoantimon2)
source(system.file("examples", "quickstart.R", package = "Neoantimon2"), local = TRUE)
results[, c("sample_id", "ligand_sequence", "mhc_key",
  "prediction_status", "grading_status")]
neo_explain(p, candidate_id = results$candidate_id[[1L]], format = "text")
```

There are **1 event, 7 origins and 4 candidates**. The hand-authored references
and 3-/5-aa peptide windows are technical demonstrations, not real human-reference
or predictor settings. No network/model is needed after installation. The script
leaves `p`, `cfg`, `events`, `origins`, `results` and `reopened` in R and prints its
temporary input/run location. Temporary storage is not a permanent archive.

Runtime R dependencies are `jsonlite`, `digest`, `stringi`, `yaml`, `xml2` and
R's `tools`; source builds can require system libraries/compilers. Predictors,
models and human reference data are acquired separately. macOS arm64 / R 4.5.1
was validated; this is not a claim of validation on every R version/platform.

## Prepare your files and identities

Keep original inputs, references and generated runs in separate locations:

```text
~/samples/     S1.vcf, S2.vcf, S3.vcf
~/neo-data/    versioned references and acquisition receipts
~/neo-runs/    generated projects and results
```

These are example locations, not a mandatory layout. Absolute paths are allowed.
Declare who each sample belongs to, independently of filenames:

| Patient | Samples | HLA declaration applies to |
| --- | --- | --- |
| P1 | S1, S2 | Both S1 and S2 |
| P2 | S3 | S3 only |

Each supported VCF has a single sample column matching its `sample_id`. A minimal
synthetic VCF looks like this (**tabs**, not spaces, separate its fields):

```vcf
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	4	.	G	A	.	PASS	.	GT	0/1
```

This tiny record is illustrative; the assembly, contig and REF base must actually
agree with the selected genome. Do not silently fix mismatches or rewrite a
multi-sample VCF header. If you begin with a multi-sample file, select each sample
with an appropriate upstream tool, preserve provenance and bind the resulting
single-sample files separately.

Supply HLA from your typing/declaration. Class I uses e.g. `HLA-A*02:01`;
Class II requires a full alpha~beta pair, e.g.
`HLA-DRA*01:01~HLA-DRB1*15:01`. A beta chain alone does not define its partner.
Optional RNA/read-count, expression, CCF and MS evidence has explicit schemas;
missing evidence is not a measured negative. The basic recipe needs none of it.

## Use R-based setup

For the supported primary reference profile, provide all nine roles: genome
FASTA + FAI; GENCODE primary GTF, transcript FASTA, protein-coding transcript
FASTA, protein FASTA and selenocysteine metadata; MANE summary + GTF.
The catalog uses GRCh38.p14 / GENCODE 50 / MANE Select 1.5. Tiny demo files are
not substitutes. Matching coordinates, identifiers and versions matter.

Inspect a chosen destination without downloading:

```r
data_dir <- path.expand("~/neo-data")
refs <- neo_check_resources(profile = "minimum", dest = data_dir)
refs$files[, c("role", "status", "integrity")]
```

If you choose to obtain this profile, explicitly run the following **several-GB**
operation. Allow download/staging space; install `samtools` separately if a genome
FAI must be created. The three primary reference selections use R, not Python:

```r
refs <- neo_download_resources(profile = "minimum", dest = data_dir)
stopifnot(refs$ok)
```

Valid cached files are reused; invalid files are retained/reported rather than
silently overwritten. Individual MANE acquisition is not a complete analysis
bundle. The returned `resources` and `scientific` fragments feed `neo_config()`.
File checks do not by themselves prove cross-reference compatibility. If needed,
`neo_check_resources(config = cfg_user, validate_reference = TRUE)` below explicitly
runs the fuller compatibility check; it can take substantial time/memory.

Native prediction is optional. `neo_check_predictors()` inspects configured tools;
`neo_setup_predictor()` returns configuration fragments for legitimately obtained
installations. For example, after obtaining your NetMHCpan executable:

```r
# Configuration template: substitute an actual vendor installation.
setup <- neo_setup_predictor("netmhcpan", executable = "/your/tools/netMHCpan")
setup$checks$predictors
```

Missing tools produce diagnostics, not manufactured scores. A detected executable
is not proof its models work. Vendor software/models are not bundled and retain
licensing/registration restrictions. MHCflurry installation/environment changes
require explicit opt-in; this guide does not install any predictor.

## Configure your own project

With the complete profile already present, replace the paths, identities and
HLA in this recipe. No predictor is requested by this configuration:

```r
data_dir <- path.expand("~/neo-data")
run_dir <- path.expand("~/neo-runs")
vcf_path <- path.expand("~/samples/S1.vcf")

refs <- neo_check_resources(profile = "minimum", dest = data_dir)
stopifnot(refs$ok)
cfg_user <- neo_config(env = character(),
  paths = list(data_root = data_dir, run_root = run_dir),
  scientific = refs$scientific, resources = refs$resources)
p_user <- neo_project(
  samples = data.frame(sample_id = "S1", patient_id = "P1"),
  hla = list(P1 = c("HLA-A*02:01", "HLA-B*07:02")),
  config = cfg_user,
  inputs = list(list(input_id = "S1-somatic", input_role = "somatic_vcf",
    sample_id = "S1", path = vcf_path)))
p_user <- neo_run(p_user)
r_user <- neo_results(p_user, sample_id = "S1")
```

If you have reviewed a valid native predictor setup, supply
`predictors = setup$predictors` as an additional `neo_config()` argument **before**
declaring a new project. Do not silently alter a completed run to add predictions.
Keep your R recipe and original declarations with the inputs. `neo_project()`
validates and declares the project; `neo_run()` performs the analysis and saves its state.

### Several samples or patients in one project

For P1/S1/S2 and P2/S3, replace the corresponding declarations as follows.
Each file must exist and have the matching sample column; `cfg_user` is the
reference configuration from above.

```r
samples <- data.frame(sample_id = c("S1", "S2", "S3"),
  patient_id = c("P1", "P1", "P2"))
hla <- list(P1 = c("HLA-A*02:01", "HLA-B*07:02"), P2 = "HLA-A*03:01")
vcfs <- file.path(path.expand("~/samples"), paste0(samples$sample_id, ".vcf"))
inputs <- lapply(seq_len(nrow(samples)), function(i) list(
  input_id = paste0(samples$sample_id[[i]], "-somatic"),
  input_role = "somatic_vcf", sample_id = samples$sample_id[[i]], path = vcfs[[i]]))
p_cohort <- neo_project(samples = samples, hla = hla, inputs = inputs, config = cfg_user)
p_cohort <- neo_run(p_cohort)
neo_results(p_cohort, sample_id = "S2")
```

HLA is keyed by patient; input records are keyed by sample. Removing rows/records
and unused patient HLA entries gives the smaller scenarios. Several sample files
do not activate longitudinal analysis: that needs explicit sample graphs and
settings. Four complete offline cohort recipes are maintained separately as
[executable examples](multi-sample-examples.md), not as another introductory guide.

## Explore the completed run

For your real project, `r_user <- neo_results(p_user)` gives candidate assessments.
These commands inspect the completed **demo** `p` so they work without your data:

```r
events <- neo_events(p)
origins <- neo_origins(p)
results <- neo_results(p)
results[, c("sample_id", "ligand_sequence", "mhc_key",
  "prediction_status", "grading_status", "final_grade")]
neo_results(p, sample_id = "sample-1")
neo_explain(p, candidate_id = results$candidate_id[[1L]], format = "text")
```

An **event** is a normalized input observation. An **origin** records a particular
sequence/translation path. A **candidate** is a peptide and HLA presenter in one
sample. One event can support multiple origins; different origins can share a
peptide. Shared sequences keep every originating relationship and do not merge
sample occurrences. Explanations connect candidates to their origins, events,
reference/provider provenance and assessment reasons.

| Output / location | Meaning |
| --- | --- |
| `neo_events(p)` / `neo_origins(p)` | Input observations / peptide-producing provenance |
| `neo_results(p)` | Native candidate assessment; optional columns depend on evidence |
| `neo_explain(p, candidate_id = ...)` | Reasons and retained candidate → origins → events trace |
| `p$path` | Saved project; `manifest/` and `intermediate/` hold verified state |
| `results/export_bundle/` inside the run | `events.tsv`, `origins.tsv`, `candidates.tsv`, unique `peptides.tsv`, `pmhc_queries.tsv`, `candidate_origin_links.tsv`, full `annotations.json` and manifest |

The export tables are annotations/queries, not a replacement for the accessor's
native assessment. Accessors verify persisted artifacts. A filtered zero-row
view differs from a missing/unexecuted output. Some fields are lists: preserve
whole tables with `saveRDS(results, "my-results.rds")`; select scalar columns
explicitly for CSV rather than discarding provenance.

### Prediction, rank and grade are different

`NOT_RUN` means no successful evaluated primary prediction is available. In the
demo none was requested; other cases can be unsupported, unavailable or failed,
so inspect provider reasons. `NOT_GRADED` means grading is inapplicable or lacks
required evidence; it is **not grade E**. An explicit numerical zero and an
observed negative are not missing values.

Somatic A–E grades and ranks require technical eligibility, the applicable primary
metric/policy and sufficient evidence. They do not estimate a universal chance
of immune response. Advanced evidence overlays do not rewrite the underlying
native grade/rank. Optional source extensions can retain source-specific
`NOT_GRADED`. Detailed thresholds and column semantics belong in the result
reference after you understand this basic interpretation.

## Add external prediction

A completed project can export without native predictions or grades:

```r
bundle <- neo_export_bundle(p)
bundle$path
```

For your real analysis, use `p_user`. Run an external predictor on its compatible
peptide/HLA queries **outside Neoantimon2**, preserving the actual program/model
version and output. If necessary, review and run a converter to produce the
standard TSV. Converter scripts adapt formats; they do not install/run engines.

BA predicts affinity (e.g. nM); EL predicts eluted-ligand/presentation behavior;
IM models immunogenicity. Their raw scores/percentiles and directions differ.
Do not rescale one into another, equate IM to an EL percentile, or automatically
request somatic A–E grades from raw scores.

For a **BigMHC-EL** converted table from your real project, this is the complete
independent-ranking import recipe. Replace the path. Choose another provider and
metric explicitly for another method; do not combine incompatible channels:

```r
bundle <- neo_export_bundle(p_user)
# After separate prediction and conversion, set this to your standard TSV:
standard_tsv <- path.expand("~/predictions/standard.tsv")
scores <- read.delim(standard_tsv, colClasses = "character", check.names = FALSE)
channel <- scores[scores$provider == "bigmhc-el" & scores$metric == "score", ]
stopifnot(nrow(channel) > 0L)
fields <- c("provider", "run_id", "task", "metric", "score_direction")
stopifnot(nrow(unique(channel[, fields])) == 1L)
policy <- c(list(mode = "external_rank"), as.list(channel[1L, fields]))
a <- neo_analyze_external(p_user, standard_tsv, policy = policy)
r_external <- neo_external_results(a)
id <- r_external$candidate_id[r_external$external_status == "evaluated"][[1L]]
neo_external_explain(a, candidate_id = id)
```

Strict matching uses peptide or peptide+full HLA identity. All matched original
Origin links remain available. The original project is required; an export alone
is insufficient. Native results do not change. External rankings/explanations
are saved as separate derived analyses and do not automatically feed all native
longitudinal, vaccine or allo modules. The
[Converter Registry README](https://github.com/hase62/Neoantimon2-Converters#readme)
is the single practical selection/conversion walkthrough, including checksums,
consent and a synthetic before/after example.

To try external conversion/import offline without a predictor or long peptides:

```r
source(system.file("examples", "external-predictions.R", package = "Neoantimon2"),
  local = TRUE)
```

That example uses synthetic values and validates origin preservation and separate
results. It does not execute a prediction model or establish prediction accuracy.

## Reopen and inspect

```r
reopened <- neo_project(path = p$path)
stopifnot(identical(neo_results(reopened), neo_results(p)))
neo_explain(reopened, candidate_id = neo_results(reopened)$candidate_id[[1L]])
```

For your durable project, retain `p_user$path` and the original configuration,
input/reference/model versions and checksums outside the run. Use a new project
when changing scientific inputs/settings. Keep original predictor output if it
cannot be regenerated. A retained run can be verified/reopened; changing or
corrupting saved artifacts is not silently accepted.

## Reference shelf

The end-to-end path above is the recommended own-data guide. Open these only
for the exact contract you need: [input roles](inputs.md), [setup options](setup.md),
[resource/tool inventory](resources-and-tools.md), [result/grade rules](results-and-explain.md),
[external schema/policies](external-predictions.md), [VCF scope](vcf-compatibility.md),
[16 APIs](api-reference.md), [advanced workflows](workflows.md),
[limits](methods-and-limitations.md) and [troubleshooting](troubleshooting.md).

Longitudinal, vaccine-set and restricted indirect-HLA/allo workflows are distinct
opt-ins. They do not imply treatment recommendations; indirect-HLA is not official
PIRCHE. Do not infer biological absence from missing evidence.

## Install a source archive

As an alternative to GitHub installation, obtain the source archive and SHA256
manifest from the v0.3.0 release, verify the checksum, then run:

```r
install.packages("Neoantimon2_0.3.0.tar.gz", repos = NULL, type = "source")
```

Install its R dependencies first; restart R before replacing a loaded package.
This alternative is not an additional prerequisite for the workflow above.

## Agent-assisted analysis

Codex starts with [AGENTS.md](../../AGENTS.md); Claude Code starts with
[CLAUDE.md](../../CLAUDE.md). Both use the same [analysis playbook](../ai/ANALYSIS_PLAYBOOK.md).
The YAML manifest is an agent declaration, not a new package input format.
Agents still need your actual declarations, resources and authorization.

Research use only. [License](../../LICENSE), [NOTICE](../../NOTICE) and
[software citation](../../CITATION.cff) apply; third-party terms remain separate.
