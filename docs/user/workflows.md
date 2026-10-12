# Supported workflows

[Getting started](getting-started.md) · [Inputs](inputs.md) ·
[Results/explain](results-and-explain.md) · [Providers](resources-and-tools.md)

## Single-sample somatic analysis

Use the primary reference profile, explicit sample/subject metadata, supported
event inputs and subject-level HLA. `neo_run()` natively normalizes events,
reconstructs compatible coding/phase origins, translates altered sequences,
generates peptide windows and exact-self comparisons, and forms resolved
peptide–presenter candidates. Noncoding, unresolved or unsupported contexts retain
their limitations rather than guaranteeing candidates.

Optional configured binding/presentation predictors may execute locally. Supported
quantitative somatic tables enable technical eligibility/grade/rank. Other
evidence remains an independent overlay or source extension: input support does
not imply an executable provider. The [installed example](getting-started.md)
demonstrates a complete native run without any predictor. Its `NOT_RUN` and
`NOT_GRADED` values are expected.

Fusion/splice inputs import explicit structure/provider output. STAR-Fusion and
Arriba are not executed; SpliceAI/STAR junction support alone does not define a
splice structure. Long-read GTF/metadata inputs preserve isoform structures,
reference compatibility, end completeness and uncertainty. Their candidates
retain source-specific `NOT_GRADED` reasons. Use explanations for their extension
graphs as well as the main candidate view.

## Presentation and recognition evidence

Presentation input roles are `presentation_settings`, `processing_output`,
`immunopeptidomics_psm`, `immunopeptidomics_mzidentml`, `apm_evidence`,
`hla_expression`, and `arcas_hla_quant`. Settings use
`presentation-run-settings-v1` with explicit `version`, provider requests,
`max_queries` and experiment definitions as needed. PSMs retain run, experiment,
IP type/class, quality and peptide modification scope. Processing providers are
import adapters/command plans; no automatic real processing-provider execution
is supplied. A pan-HLA ligand observation does not identify a unique presenter.

Recognition roles are `tcr_rearrangement`, `tcr_airr`, `tcr_repertoire`,
`tcr_cell`, `tcr_receptor`, `tcr_reactivity`, `tcr_vdjdb`,
`recognition_settings`, and `recognition_output`. Settings use
`recognition-run-settings-v1` with versioned providers and explicit candidate,
profile and query/pairing limits. Default limits are 10,000 queries, 1,000 pairing
combinations and batch size 256. Import supported AIRR/repertoire/assay/provider
records; do not infer chain pairing, receptor functionality, recognition or
assay completeness from missing data. These providers do not have automatic real
execution through the current recognition interface.

Use installed versioned schemas for exact settings shapes and the
[provider table](resources-and-tools.md#external-tools-and-provider-inputs) for
accepted provider families. These independent dimensions are neither a universal
score nor a replacement for the somatic grade.

## Longitudinal and multi-sample analysis

An input with role `longitudinal_settings` loads a JSON document with
`version = "longitudinal-integration-1.0.0"` and explicit arrays `subjects`,
`contexts`, `relations`, optional `treatments`, `evaluability`, `providers`,
and `sample_runs`. Every project sample requires one context binding it to the
same declared subject. Relations remain within that subject: `same_subject`,
`same_lesion`, `different_region`, `primary_to_metastasis`, `pre_to_post`,
`serial_sample`, `matched_normal`, `technical_replicate`.

Contexts may declare one time axis: ISO date/UTC `collection_time`, numeric
`relative_time` with `time_origin` and `time_unit`, or `ordinal_time`. Comparable
explicit times establish order; absent/incompatible time stays unordered.
Neither a filename nor row order asserts chronology.

An aggregate project binds a completed sample run for **every** context. These
are explicit reproducible inputs, with individual graphs retained; regenerate
them from their original recipes if necessary. They are not accidental historic
development dependencies. `neo_results()` adds contextual evidence while keeping
sample-specific somatic priority, grade, rank and order. `neo_multisample()`
provides subject-scoped candidate trajectories and typed evaluability cells;
`neo_explain()` traces referenced sample graphs. Inspect original projects for
their individual event/origin tables. A single sample or matched-normal-only
comparison does not create a multi-tumor overlay.

This runnable example follows the installed quick start and repeats its same
**synthetic** event in a second sample. It demonstrates relationships, not
biological evolution or a measured absence:

```r
vcf2 <- file.path(input_root, "synthetic-2.vcf")
writeLines(sub("sample-1$", "sample-2", readLines(vcf)), vcf2)
p2 <- neo_run(neo_project(
  samples = data.frame(sample_id = "sample-2", patient_id = "synthetic-subject"),
  config = cfg, hla = list("synthetic-subject" = "HLA-A*02:01"),
  inputs = list(list(input_id = "synthetic-variant-2", input_role = "somatic_vcf",
    sample_id = "sample-2", path = vcf2))))
settings <- list(version = "longitudinal-integration-1.0.0",
  subjects = list(list(subject_id = "synthetic-subject")),
  contexts = list(
    list(sample_id = "sample-1", subject_id = "synthetic-subject",
      role = "tumor", ordinal_time = 1),
    list(sample_id = "sample-2", subject_id = "synthetic-subject",
      role = "tumor", ordinal_time = 2)),
  relations = list(list(type = "serial_sample", from = "sample-1", to = "sample-2")),
  sample_runs = list(list(sample_id = "sample-1", path = p$path),
    list(sample_id = "sample-2", path = p2$path)))
settings_file <- file.path(input_root, "longitudinal.json")
jsonlite::write_json(settings, settings_file, auto_unbox = TRUE, pretty = TRUE)
cfg_long <- neo_config(env = character(), paths = list(
  data_root = input_root, input_root = input_root,
  run_root = file.path(example_root, "longitudinal-runs")))
p_long <- neo_run(neo_project(config = cfg_long,
  samples = data.frame(sample_id = c("sample-1", "sample-2"),
    patient_id = rep("synthetic-subject", 2)),
  inputs = list(list(input_id = "longitudinal", input_role = "longitudinal_settings",
    sample_id = "sample-1", path = settings_file, version = "1"))))
long_results <- neo_results(p_long)
neo_multisample(p_long, patient_id = "synthetic-subject")
neo_explain(p_long, candidate_id = long_results$candidate_id[[1L]])
stopifnot(identical(neo_results(p), results))
```

Evaluability declarations distinguish `evaluated_positive`, `evaluated_negative`,
`not_evaluated`, `insufficient_evidence`, `not_applicable` and `ambiguous`.
Positive/negative assertions require explicit adequate opportunity, assay,
evidence reference and policy version. Gain/loss descriptions require explicit
chronology and supported absence; missing observations do not establish loss.

Generic `clonal_evidence`, `pyclone_vi` and `clonal_mapping` inputs import
provider results. Mutation/cluster IDs are scoped by provider, version, run,
subject and settings; exact mapping to events/origins is required. Provider
prevalence is separate from event CCF. No averaging, cross-run global clone,
de novo clustering or universal trajectory score is inferred.

## Declarative vaccine-set selection

`neo_config(vaccine_set = ...)` accepts a list or JSON/YAML policy filename.
Use `version = "vaccine-set-policy-1.0.0"` and explicit `max_peptides`.
Omitting the policy leaves ordinary outputs unchanged. Policy fields can declare
Class I/II minima/caps, required presenters, per-presenter/redundancy caps,
length limits, literal forbidden residues/motifs, candidate ID/key allowlists
and blocklists, confirmed exact-self exclusion, and required evidence.
Supported evidence requirements are `exact_self_evaluated`,
`longitudinal_evaluable`, `linked_clonal`, `ms_observed`. Missing required evidence
is insufficient evidence, not a satisfied negative constraint.

After the installed quick start, this creates a separate **synthetic** selection:

```r
cfg_set <- neo_config(env = character(), paths = cfg$paths,
  scientific = cfg$scientific, checksum = cfg$checksum,
  resources = cfg$resources, peptides = cfg$peptides,
  vaccine_set = list(version = "vaccine-set-policy-1.0.0",
    max_peptides = 2L, fill_to_capacity = TRUE))
p_set <- neo_run(neo_project(config = cfg_set,
  samples = data.frame(sample_id = "sample-1", patient_id = "synthetic-subject"),
  hla = list("synthetic-subject" = "HLA-A*02:01"),
  inputs = list(list(input_id = "synthetic-variant", input_role = "somatic_vcf",
    sample_id = "sample-1", path = vcf))))
set_results <- neo_results(p_set)
attr(set_results, "vaccine_set")
neo_explain(p_set, candidate_id = set_results$candidate_id[[1L]])$vaccine_set
stopifnot(identical(set_results$final_grade, results$final_grade))
```

The current `lexicographic-greedy-v1` algorithm is a deterministic heuristic,
not a global mathematical optimum. Capacity counts distinct ligand/presenter
keys, with separate subject sets. Eligibility, feasibility, coverage, redundancy,
escape evidence, marginal vectors and the decision trace are retained in the
`vaccine-set-selection-1.0.0` audit. Declarative minima may be infeasible; inspect
the feasibility/selection reasons rather than assuming a full set. Selection
does not mutate somatic grades/ranks or turn a selected ungraded candidate into
evidence of clinical suitability.

## Allo contexts and strict indirect-HLA

Allo support is enabled through configuration; there is no additional public
function or `mode = "allo"`. `neo_project()` still uses its supported `somatic`
mode. Construction validates and derives the separate opt-in graph, which
query-free `neo_explain()` can inspect before ordinary execution. A real
no-event allo-only `neo_run()` is not supported by the current run route.

Prepare `allo` with `version = "allo-foundation-1.0.0"`: explicit subjects,
transplant contexts, typing scenarios, a normalized literal reference JSON
(`path` and SHA256), and declared self-search universes. Roles/direction must
explicitly distinguish source, self, APC and responding subject. A source
sequence is not a presenter; a typing hypothesis is not a resolved choice.
Raw IMGTHLA acquisition supplies upstream data, not this normalized reference
or missing resolution/coverage attestations.

Add `indirect_hla` with `version = "indirect-hla-baseline-1.0.0"`,
`profile = "neo-literature-derived-indirect-hla-ii-dr-strict-v1"`, explicit
presenter records and a pinned external precomputed provider JSON. Presenter
records bind a context scenario (`scenario_id`), subject, literal presenter and
`resolution_policy = "explicit_typed_alpha_beta"`. Authorized external output
retains query/result/executable/model identities, license authorization and
execution provenance; synthetic evidence is labeled as such and never passed
off as executed vendor evidence.

The strict profile is Neoantimon2's **literature-derived implementation**, not
official PIRCHE, a PIRCHE reproduction or PIRCHE-equivalent software:

- The current executable strict scenario is solid-organ donor-to-recipient
  indirect analysis. Source/self loci are A, B, C, DRB1 and DQB1 with literal,
  complete precursor sequence/coverage requirements.
- Presenters are explicit DRA~DRB1 pairs. DQ, DP and DRB3/4/5 presenters are
  outside this strict profile; alpha chains are not inferred.
- Generate contiguous 15-mers with stride 1; no sequence padding or homology
  completion. Retain supported primary forward 9-aa core and offset mapping.
- Import NetMHCIIpan **3.2-family BA**, nM semantics; the qualifying threshold
  is strictly **< 1000 nM**. Exactly 1000 nM is excluded. EL ranks or the somatic
  4.3 execution adapter are not substitutes for this provider contract.
- Count unique qualifying non-self core/presenter keys after the supported
  source union and qualifying same-presenter self subtraction. Exhaustive
  declared self scope is required. Source subtotals can overlap; do not sum them
  to produce a context total.
- `complete_count` is `NULL` if required scope/evidence is incomplete. Known
  partial positives are retained without being presented as a complete count.

With your validated manifest objects bound to `allo_input` and `indirect_input`,
and your chosen roots bound to `data_dir` and `run_dir`, the public entry is:

```r
cfg_allo <- neo_config(env = character(),
  paths = list(data_root = data_dir, run_root = run_dir),
  allo = allo_input, indirect_hla = indirect_input)
p_allo <- neo_project(config = cfg_allo)
neo_explain(p_allo)$indirect_hla_explanation
```

These are user-supplied manifests, not hidden built-in patient inputs. To use
`neo_results()`, bind these opt-ins together with the ordinary event/reference
configuration, samples/HLA and inputs, then run that project. Its candidate table
retains separate `indirect_hla_result` and `indirect_hla_baseline` attributes;
somatic row counts do not count indirect members. Construction-only explanations
do not require a fabricated somatic event. Read completeness, comparison/reporting
scope and reasons; do not convert a null count to zero or clinical risk.

## Independent native allo evidence

Optional `allo_native_evidence` binds versioned observation files and explicit
context/sample/compartment scopes. The model is `neo-allo-hla-evidence-v1`, schema
`neo-allo-hla-evidence-1.0.0`. Its independent dimensions are source expression,
presenter expression, allele-specific HLA LOH, antigen-processing machinery,
MS observation and TCR recognition. Observation assignment requires qualified
context/subject/presenter/compartment applicability; wrong-scope evidence is not
automatically assigned.

Retained assignment states include `applicable_supported`, `applicable_negative`,
`ambiguous`, `insufficient_scope`, `not_evaluated`, `not_applicable`, `failed`.
Missing evidence is not negative. Inspect the `allo_native_evidence` result
attribute or `neo_explain(project)$allo_native_evidence`. This overlay does not
change strict indirect candidate membership, complete counts or somatic grade/rank,
and its independent evidence dimensions are not collapsed into a universal score.
