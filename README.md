# Neoantimon2

Neoantimon2 is an R package for reconstructing **neoantigen candidates** from
supported genomic inputs and interpreting them with prediction and experimental
evidence. It connects each candidate peptide to its genomic events, sample,
patient, sequence origins and declared HLA molecules, so you can trace how a
candidate was produced and why it received a particular assessment.

You can run configured predictors locally through the **native workflow**, or
export peptide–HLA queries for other methods and import their predictions through
the **external workflow**. Both preserve peptide origins and evidence provenance.
Ranking and grading depend on the selected workflow and available evidence.

Neoantimon2 is a **research tool**. Generated candidates, predictions and grades
do not establish clinical immunogenicity, treatment efficacy or clinical risk.

[Install](#installation) · [Quick Start](#quick-start) ·
[Your own data](#analyze-your-own-data) ·
[Multiple patients and samples](docs/user/multi-sample-examples.md) ·
[Documentation](#documentation)

## How it works

![Neoantimon2 workflow: sample-linked variants and declared patient HLA become traceable peptide–HLA candidates. Native prediction and a separate external export/import route support interpretation; specialized workflows require explicit settings.](docs/assets/neoantimon2-overview.svg)

Events lead to reconstruction/translation origins, peptide sequences and
peptide–HLA candidates. Several origins can yield the same peptide; each
candidate keeps its origin links for inspection.

**Native:** `neo_run()` uses configured local predictors and supplied evidence.
**External:** export a completed project, run the model separately, and import
standardized predictions by exact peptide or peptide–HLA matching. External
results and explanations remain separate from native results.

[Figure guide and accessible description](docs/user/workflow-overview.md) ·
[Full-size SVG](docs/assets/neoantimon2-overview.svg) ·
[PNG / PDF / editable source](docs/user/workflow-overview.md#maintain-or-reuse-the-figure) ·
[Concrete input/output examples](docs/user/input-output-examples.md)

For Codex or Claude Code, start with [agent instructions](docs/ai/AI_INSTRUCTIONS.md)
and the [shared analysis playbook](docs/ai/ANALYSIS_PLAYBOOK.md). These repository
guides use the released v0.3.0 APIs; they do not replace required inputs or
authorization.

## Installation

Use **R 4.1.0 or later**. Install the pinned
[v0.3.0 release](https://github.com/hase62/Neoantimon2/releases/tag/v0.3.0):

```r
install.packages("remotes")
remotes::install_github("hase62/Neoantimon2", ref = "v0.3.0", upgrade = "never")
```

- **Offline Quick Start:** the package and its R dependencies are sufficient.
  Tiny synthetic references are included for demonstration.
- **Your genomic data:** obtain a compatible reference bundle and prepare your
  input files, sample metadata and HLA declarations separately.
- **Optional prediction:** obtain and configure the chosen predictors and their
  models under their providers' terms. Licensed predictors and human reference
  genomes are **not bundled** with the package.

Runtime dependencies are `jsonlite`, `digest`, `stringi`, `yaml`, `xml2`, and R's
`tools`. Source builds of dependencies may require compilers/system libraries.
Suggested packages are for development/tests. Public tag and release-archive
installation were verified on macOS arm64 with R 4.5.1; this does not establish
validation on every platform or the minimum supported R version. See
[getting started](docs/user/getting-started.md) and
[resources and tools](docs/user/resources-and-tools.md).

Version 0.3.0 includes [R-based resource and predictor setup](docs/user/setup.md).
The earlier [v0.2.0 release](https://github.com/hase62/Neoantimon2/releases/tag/v0.2.0)
remains available with its original artifacts; it does not include the four setup APIs.

## Quick Start

After installation, run this complete offline example:

```r
library(Neoantimon2)
source(system.file("examples", "quickstart.R", package = "Neoantimon2"), local = TRUE)
neo_events(p)
neo_origins(p)
neo_results(p, prediction_status = "NOT_RUN")
neo_explain(p, candidate_id = results$candidate_id[[1L]], format = "text")
```

The example creates **1 event, 7 peptide origins and 4 peptide–HLA candidates**.
It checks that results remain identical after reopening the run. No predictor is
requested: prediction is `NOT_RUN`, and grading is `NOT_GRADED`. The candidates
are synthetic technical outputs, not biologically validated neoantigens.

The script prints the temporary directory containing its inputs and run and
leaves `p`, `events`, `origins` and `results` available for inspection. It needs no
network or licensed software after installation. Its tiny references and 3-aa/5-aa
peptide windows are test settings, not settings for real prediction. See
[results and explanations](docs/user/results-and-explain.md) for interpretation.

## Analyze your own data

Prepare the following before starting:

- A `samples` data frame with a unique `sample_id` for each sample and an explicit
  `patient_id` linking it to its patient.
- Patient-level `hla`, a named list keyed by `patient_id`, to construct
  peptide–HLA candidates. HLA typing is supplied by you, not inferred.
- An `inputs` list of separately declared records, each with `input_id`,
  `input_role`, `sample_id` and `path`. For VCFs, use one single-sample file per
  declared sample; its sample column, when present, must match `sample_id`.
- A complete, compatible reference bundle and writable run location, configured
  through `neo_config()`. The
  [input guide](docs/user/inputs.md) specifies the reference roles and formats.

Create a project with `neo_project(samples = ..., hla = ..., inputs = ...,
config = ...)`, execute it with `p <- neo_run(p)`, and inspect candidates with
`neo_results(p, sample_id = "S1")`. This describes the API sequence; the
[own-data recipe](docs/user/getting-started.md#configure-your-own-project)
provides the configuration details.

**Several patients or samples can share one project.** See the
[four complete, executable examples](docs/user/multi-sample-examples.md) for one
patient/one sample through two patients/three samples, with generated VCFs and
offline reference setup. Multiple separately bound VCFs are supported;
**multiple samples inside one VCF are not supported**.

Ordinary sample management does not imply a longitudinal comparison.
`neo_multisample()` exposes configured subject-scoped longitudinal summaries;
it is empty for the ordinary projects in those examples. See the
[longitudinal workflow](docs/user/workflows.md#longitudinal-and-multi-sample-analysis)
when you need trajectories, explicit time or evaluability across samples.

## Major capabilities

- **Somatic analysis:** reconstruct peptides from supported variants, preserving
  event, transcript, phase and normal-comparison provenance.
- **Native prediction:** execute configured NetMHCpan, NetMHCIIpan or MHCflurry
  adapters within their supported scope; software/models are acquired separately.
- **External interoperability:** export verified queries and origins, import
  standardized predictions, or explicitly select local R/Python converters.
  Arbitrary raw scores are not percentiles and do not establish scientific grades.
- **Evidence and explanations:** inspect provider provenance, missingness,
  evaluability, somatic assessments and supported evidence overlays.
- **Longitudinal analysis:** compare explicitly related samples within subjects,
  retaining sample-specific evidence and distinguishing missing from negative.
- **Vaccine-set selection:** select sets under declared constraints using the
  supported deterministic heuristic; selections are not treatment recommendations
  or a guaranteed global optimum.
- **Restricted indirect-HLA/allo analysis:** use the separate literature-derived
  profile and optional native allo evidence overlays. This is not an official
  PIRCHE implementation; its scientific restrictions remain explicit.

These workflows have separate contracts and prerequisites. External imports do
not automatically feed every vaccine, longitudinal, receptor, presentation or
allo module. See [supported workflows](docs/user/workflows.md) and
[methods and limitations](docs/user/methods-and-limitations.md).

## Documentation

Start with the [documentation index](docs/README.md), or go directly to:

- [Getting started](docs/user/getting-started.md) ·
  [Multiple patients and samples](docs/user/multi-sample-examples.md)
- [Inputs and configuration](docs/user/inputs.md) ·
  [Resources and tools](docs/user/resources-and-tools.md)
- [Workflow overview](docs/user/workflow-overview.md) ·
  [Workflows](docs/user/workflows.md) ·
  [External prediction interoperability](docs/user/external-predictions.md)
- [Public API reference](docs/user/api-reference.md) ·
  [Results and explanations](docs/user/results-and-explain.md)
- [Methods and limitations](docs/user/methods-and-limitations.md) ·
  [Reproducibility](docs/user/reproducibility.md) ·
  [Troubleshooting](docs/user/troubleshooting.md)
- [Contributing and bug reports](CONTRIBUTING.md)

Installed help is available through `help(package = "Neoantimon2")`.

The optional, separately maintained
[Community Converter Registry](https://www.hase62.jp/Neoantimon2-Converters/)
provides converter discovery. Versioned IDs, SHA256 verification, explicit
first-use execution consent and approved local cache reuse are supported.
Downloaded converters run with your privileges; checksums and community OK/NG
feedback are not security guarantees or official scientific validation. Native
analysis and direct standardized import do not require the Registry. See the
[external guide](docs/user/external-predictions.md) for its index and usage.

## Scientific limitations

Missing or unavailable evidence is distinct from a valid numerical zero or an
evaluated negative. `NOT_RUN` does not establish a biological negative;
`NOT_GRADED` is not grade E. No universal immunogenicity or clinical-risk score
is supplied. Real licensed predictors, clinical outcomes and other platforms
retain their documented validation limits. Read the
[methods and limitations](docs/user/methods-and-limitations.md) before interpreting
research outputs.

## Citation

Use [CITATION.cff](CITATION.cff) for the software citation and record the exact
version/commit, input/reference versions and external providers. A dedicated
Neoantimon2 paper is planned; the software citation is not a peer-reviewed
validation publication. The predecessor Neoantimon paper does not validate this
package.

## License

[PolyForm Noncommercial 1.0.0](LICENSE) governs this package; retain the required
[NOTICE](NOTICE). Commercial use requires separate permission or licensing from
Takanori Hasegawa: <hasegawa.t.d968@m.isct.ac.jp>. An inquiry does not grant permission.
Third-party tools, resources and converters retain their own license boundaries.
This GitHub release is not a CRAN submission.
