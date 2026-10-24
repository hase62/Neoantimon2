# Workflow overview

[Getting started](getting-started.md) · [Supported workflows](workflows.md) ·
[External predictions](external-predictions.md) · [Methods and limitations](methods-and-limitations.md)

![Neoantimon2 workflow: declared samples, patients, references and HLA lead to events, origins, peptides and peptide–HLA candidates. Native interpretation and independent external prediction retain origin links; specialized workflows have separate prerequisites.](../assets/neoantimon2-overview.svg)

[Full-size SVG](../assets/neoantimon2-overview.svg) ·
[High-resolution PNG](../assets/neoantimon2-overview.png) ·
[Vector PDF](../assets/neoantimon2-overview.pdf) ·
[Editable diagrams.net source](../assets/source/neoantimon2-overview.drawio)

## Read the figure

Read the blue **A — Native** row from left to right. Samples are explicitly
assigned to patients, and each patient's supplied HLA declarations apply to that
patient's samples. HLA typing is not performed automatically. VCF records are
bound separately to samples; a VCF with multiple sample columns is unsupported.
Reference resources are acquired separately and selected explicitly.

An event can contribute to several reconstruction and translation paths. A path
can contain several events, and several paths can yield the same peptide.
Peptide–HLA candidates retain all contributing origin links within their sample
and analysis context. A shared sequence does not merge sample occurrences.
Follow `neo_events()`, `neo_origins()`, `neo_results()` and `neo_explain()` to
inspect this trace, including uncertainty and contexts that cannot yield a
candidate.

Solid arrows show implemented processing or data flow. Dashed arrows show
optional evidence, configured predictors, or the conditional export from a
completed project boundary. Native prediction scores and somatic grades are not
prerequisites for export: runs with `NOT_RUN` / `NOT_GRADED` can still provide
a verified bundle. The reference arrow supplies the reconstruction pipeline;
it is not a separate variant caller. Boxes and lettered lane titles distinguish
the routes without relying on color.

The orange **B — External** row starts with a completed original project.
`neo_export_bundle()` preserves queries, events, origins and their associations.
Run the chosen model separately, then supply its output in the standard format
to `neo_analyze_external()`, directly or through an explicitly selected converter.
The optional [Converter Registry](external-predictions.md#local-and-registry-converters)
helps obtain reusable converters; it does not run prediction models. Strict
matching uses a literal peptide or literal peptide plus full presenter. Class II
requires its complete alpha~beta declaration; missing chains are not guessed.

External evidence and analyses persist separately. Inspect them with
`neo_external_results()` and `neo_external_explain()`; they do not replace native
results or automatically feed other evidence/selection overlays. External ranking
requires an explicit policy. External somatic grading additionally requires
truly declared binding/presentation percentiles and the original eligible core
somatic technical/evidence prerequisites. Arbitrary raw scores cannot produce
somatic A–E grades. See the [external guide](external-predictions.md) for the
executable example and exact contract.

## Prediction and evidence boundaries

Native prediction uses separately obtained, configured software and models.
These are the adapters summarized by the top-right box:

| Native adapter | Supported family | HLA class / lengths | Implemented tasks and returned metrics |
| --- | --- | --- | --- |
| NetMHCpan | 4.1 / 4.2 | I / 8–14 aa | Binding and presentation; BA nM/rank, EL score/rank |
| MHCflurry | 2.x | I / 8–15 aa | Binding and presentation; BA nM/rank, presentation score/percentile |
| NetMHCIIpan | 4.3 | II / 9–25 aa | Binding and presentation; BA nM/rank, EL score/rank, supported core mapping |

Presenter/residue, model and task checks still apply; the table is not a promise
that every allele or installation is supported. Predictor configuration checks
are distinct from successful scientific execution. NetMHCpan/NetMHCIIpan retain
the vendor's license/registration requirements. The R package does not bundle
them or human references. See [setup](setup.md) and
[resources and tools](resources-and-tools.md).

The optional evidence box summarizes supported observations such as RNA,
read-count and mass-spectrometry evidence. These have explicit input/provider and
scope requirements; not all are automatically computed. Presentation and
recognition observations retain independent interpretations. Ranking and grades
are conditional on technical eligibility and the applicable policy; a peptide
alone is not a validated neoantigen.

Missing/unavailable evidence is distinct from an evaluated zero or a measured
negative. `NOT_RUN` is not a biological negative, and `NOT_GRADED` is not grade E.
The figure groups evidence acquisition for readability without replacing the
[documented evidence states](methods-and-limitations.md).

## Specialized workflows are opt-ins

The three gray boxes show separate capabilities, not automatic branches from
external results:

- **Longitudinal analysis** requires explicitly related contexts, settings and
  sample graphs from the current project or explicitly associated completed
  sample runs. Informal timepoint labels alone do not activate it. It preserves
  sample-specific native grades/ranks and
  distinguishes missing observations from supported absence. Ordinary cohort
  input is not a longitudinal analysis.
- **Vaccine-set selection** uses native candidate/evidence projections and an
  explicit policy. Its deterministic heuristic does not change native grades or
  ranks, guarantee a global optimum, or establish treatment suitability.
- **Restricted indirect-HLA/allo analysis** has its own inputs, profile and
  evidence boundaries. The executable strict profile covers solid-organ
  donor-to-recipient analysis, source/self loci A/B/C/DRB1/DQB1 and explicitly
  typed DRA~DRB1 presenters. It uses 15-mers, primary forward 9-aa cores and pinned
  imported NetMHCIIpan 3.2-family BA evidence with a strict <1000 nM threshold.
  The native somatic 4.3 adapter is not a substitute. Incomplete required scope
  yields a null complete count, not zero. Native allo observations are independent
  overlays. This is not official PIRCHE, PIRCHE-equivalent software or a clinical
  risk model.

See [supported workflows](workflows.md) for these prerequisites, the separate
construction/query route and the absence of a no-event allo-only `neo_run()`.
The diagram does not imply that each workflow connects to every other workflow.

## Maintain or reuse the figure

The vector SVG is the primary artifact. Its text/shapes are editable and need no
network, scripts, remote fonts or external images to render. The PNG is
3600 × 2520 pixels (3×); the PDF retains vector objects and selectable text.
Open the uncompressed `.drawio` file in diagrams.net to edit individual groups,
labels, boxes and connectors. No account or hosted copy is required if using a
local diagrams.net installation. Edits are not sent to a service by the build
script.

The source uses groups, rectangles, plain text and polyline connectors with
explicit ports/waypoints. Preserve those primitives and the `neoLayer=background`
container setting when editing. The standard-library Python renderer implements
this small subset; it is not a general diagrams.net renderer. Text line breaks
are explicit. Recheck the output after edits, especially if another editor
changes geometry or styles.

From the repository root, rebuild the SVG with Python 3:

```sh
python3 scripts/development/build_workflow_figure.py \
  --source docs/assets/source/neoantimon2-overview.drawio \
  --output docs/assets/neoantimon2-overview.svg
```

Optional PNG/PDF exports and retained size previews require Chrome/Chromium.
Choose an output directory outside the repository for the previews:

```sh
python3 scripts/development/build_workflow_figure.py \
  --source docs/assets/source/neoantimon2-overview.drawio \
  --output docs/assets/neoantimon2-overview.svg \
  --png docs/assets/neoantimon2-overview.png \
  --pdf docs/assets/neoantimon2-overview.pdf \
  --work-dir ~/codex/runs/Neoantimon2/figure-preview
```

Use `--browser /path/to/chrome` if detection is insufficient. Browser export was
verified with macOS Chrome; other browsers/platforms are not claimed tested.
The script measures text bounds and retains 1200-, 900- and 600-pixel previews.
Visual inspection is still necessary. SVG output is deterministic; browser/PDF
metadata and font rendering can vary by environment. At narrow screen widths,
open the full-size SVG or PNG to read the detailed notes. The accompanying text
and README workflow description remain available to screen readers.

The earlier [workflow.svg](../assets/workflow.svg) remains accessible for existing
links. Use the current overview for the full Native/External boundaries. Reuse
remains subject to the repository [license](../../LICENSE) and
[attribution notice](../../NOTICE); third-party tools retain their own terms.

## Source-to-figure audit

These implementation anchors support the main connections. The figure does not
add public APIs or change an analysis contract.

| Figure connection or boundary | Implementation and maintained guide |
| --- | --- |
| Samples → patients; patient HLA → presenters | [`project.R`](../../R/project.R), [`peptide_records.R`](../../R/peptide_records.R): `peptide_binding()`, `peptide_construct()`; [cohort guide](multi-sample-examples.md) |
| References → normalized events and reconstruction/translation paths | [`events-import.R`](../../R/events-import.R), [`reconstruction_runtime.R`](../../R/reconstruction_runtime.R), [`translation_runtime.R`](../../R/translation_runtime.R); [input contract](inputs.md) |
| Origins → shared sequences → sample/context peptide–HLA candidates, with retained links | [`peptides.R`](../../R/peptides.R): `peptide_enumerate()`; [`peptide_records.R`](../../R/peptide_records.R): `peptide_origin()`, `peptide_candidate()`, `peptide_construct()` |
| Configured native tools → supported prediction evidence | [`predictor_registry.R`](../../R/predictor_registry.R): `predictor_registry()`, `predictor_candidate_status()`; [providers](resources-and-tools.md) |
| Evidence → eligible native priority/results/explanations | [`somatic_priority.R`](../../R/somatic_priority.R), [`results.R`](../../R/results.R), [`explain.R`](../../R/explain.R); [interpretation](results-and-explain.md) |
| Completed project → verified export graph with all origin links | [`external_export.R`](../../R/external_export.R): `external_completed()`, `external_graph()`, `external_merge_graph()`, `neo_export_bundle()` |
| Separate model output → optional conversion → exact matching | [`external_converter.R`](../../R/external_converter.R), [`external_predictions.R`](../../R/external_predictions.R); [external schema](external-predictions.md) |
| Imported evidence → separate derived view; eligible percentile grading only | [`external_analysis.R`](../../R/external_analysis.R): `external_policy()`, `external_derive()`, `neo_analyze_external()`; [policy boundaries](external-predictions.md#interpretation-policies) |
| Longitudinal / selection preserve native priority | [`longitudinal_runtime.R`](../../R/longitudinal_runtime.R), [`vaccine_runtime.R`](../../R/vaccine_runtime.R), [`vaccine_selector.R`](../../R/vaccine_selector.R); [workflows](workflows.md) |
| Indirect-HLA own profile, presenter/core/self scope and null-count semantics | [`indirect_contract.R`](../../R/indirect_contract.R), [`indirect_provider.R`](../../R/indirect_provider.R), [`indirect_baseline.R`](../../R/indirect_baseline.R); [restricted profile](workflows.md#allo-contexts-and-strict-indirect-hla) |
| Missingness and research limitations | [`peptide_records.R`](../../R/peptide_records.R), [`external_predictions.R`](../../R/external_predictions.R); [methods and limitations](methods-and-limitations.md) |
