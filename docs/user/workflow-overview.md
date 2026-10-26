# Workflow figure reference

The [README](../../README.md) is the beginner explanation; the
[practical guide](getting-started.md) is the continuous own-data workflow.
This page maintains the figure's meaning, exports and implementation mapping.

![Five stages: user inputs, sequence reconstruction, peptide–HLA candidates, separate native/external prediction and results with explanations.](../assets/neoantimon2-overview.svg)

[SVG](../assets/neoantimon2-overview.svg) · [PNG](../assets/neoantimon2-overview.png) ·
[Vector PDF](../assets/neoantimon2-overview.pdf) ·
[Editable source](../assets/source/neoantimon2-overview.drawio)

## Read the figure

Read the five stages left to right. Inputs bind sample variants to explicit
patients, patient-level HLA and compatible references. The reconstruction fork
illustrates multiple possible origins; it does not impose a one-event/one-path
relationship. Several events may contribute to a path, and different paths may
share a peptide. Candidates preserve every contributing Origin within their
sample/context. The peptide strip is a conceptual sequence, not a prediction.

Native prediction uses configured supported tools. External prediction instead
exports a **completed project**, runs a model separately, converts its output if
needed, and imports it strictly into that original project. Native scores or
grades are not prerequisites for export. The two labeled channels stay separate;
external evidence does not overwrite native results. Arrows show the conceptual
processing order, not automatic execution of every optional step.

Results expose available assessment and explanations. The small bars symbolize
prioritization only where eligible; they are not measured scores. Missing evidence
is not zero or a negative; raw IM scores do not automatically get A–E grades.
The README defines these terms before pointing to detailed contracts.

Advanced longitudinal, vaccine-set and restricted indirect-HLA/allo modules are
intentionally outside the beginner figure. They retain their separate opt-in
contracts in [advanced workflows](workflows.md), including original-source,
presenter/metric restrictions and lack of automatic integration with all external
analyses. Restricted indirect-HLA is not an official PIRCHE implementation.

## Maintain or reuse the figure

The vector SVG is the primary artifact. Its text/shapes are editable and need no
network, scripts, remote fonts or external images to render. The PNG is
3600 × 1320 pixels (3×); the PDF retains vector objects and selectable text.
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
  --work-dir ~/neo-figure-preview
```

Use `--browser /path/to/chrome` if detection is insufficient. Browser export was
verified with macOS Chrome; other browsers/platforms are not claimed tested.
The script measures text bounds and retains 1200-, 900- and 600-pixel previews.
Visual inspection is still necessary. SVG output is deterministic; browser/PDF
metadata and font rendering can vary by environment. The main figure is designed for 800–900 px README viewing. At narrow widths,
open the full-size asset for labels. The accompanying text
and README workflow description remain available to screen readers.

The earlier [workflow.svg](../assets/workflow.svg) remains accessible for existing
links. Use the current five-stage overview as the beginner figure. Reuse
remains subject to the repository [license](../../LICENSE) and
[attribution notice](../../NOTICE); third-party tools retain their own terms.


## Source-to-figure audit

These implementation anchors support the main connections and documented boundaries. The figure does not
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
| Advanced reference only (not shown): longitudinal / selection preserve native priority | [`longitudinal_runtime.R`](../../R/longitudinal_runtime.R), [`vaccine_runtime.R`](../../R/vaccine_runtime.R), [`vaccine_selector.R`](../../R/vaccine_selector.R); [workflows](workflows.md) |
| Advanced reference only (not shown): indirect-HLA own profile and null-count semantics | [`indirect_contract.R`](../../R/indirect_contract.R), [`indirect_provider.R`](../../R/indirect_provider.R), [`indirect_baseline.R`](../../R/indirect_baseline.R); [restricted profile](workflows.md#allo-contexts-and-strict-indirect-hla) |
| Missingness and research limitations | [`peptide_records.R`](../../R/peptide_records.R), [`external_predictions.R`](../../R/external_predictions.R); [methods and limitations](methods-and-limitations.md) |
