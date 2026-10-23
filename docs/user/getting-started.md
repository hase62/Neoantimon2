# Getting started

[README](../../README.md) · [Inputs/configuration](inputs.md) ·
[Resources/tools](resources-and-tools.md) · [Results](results-and-explain.md)

## Install and run the offline example

R >= 4.1.0 is required. Non-commercial academic and research use is welcome under
the [PolyForm Noncommercial License 1.0.0](../../LICENSE). See the
[README licensing policy](../../README.md#license) for commercial licensing
inquiries. **Neoantimon2 0.2.0 is released** on
[GitHub](https://github.com/hase62/Neoantimon2/releases/tag/v0.2.0). Install the
pinned public version:

```r
install.packages("remotes")
remotes::install_github("hase62/Neoantimon2@v0.2.0", upgrade = "never")
library(Neoantimon2)
stopifnot(as.character(packageVersion("Neoantimon2")) == "0.2.0")
```

The public tag and release-archive installation were verified on macOS arm64
with R 4.5.1. For a downloaded archive, use
`R CMD INSTALL --library=/your/test/library Neoantimon2_0.2.0.tar.gz`.
See the [installation notes](../../README.md#installation) for dependencies
and platform limitations.

The package's Imports are `jsonlite`, `digest`, `stringi`, `yaml`, `xml2`, and
`tools`. Suggested packages `testthat`, `jsonvalidate`, and `pkgload` are for
development/tests, not ordinary analysis. R's `tools` ships with R. Predictors,
their models and public references are separate resources. Source installation
of R dependencies may require system libraries/compilers. This documentation
does not establish support for untested operating systems.

Run the installed, synthetic/test-scale/non-clinical example:

```r
library(Neoantimon2)
source(system.file("examples", "quickstart.R", package = "Neoantimon2"), local = TRUE)
```

It requires no checkout, downloads, licensed software or network after package
installation. Installed fixture files are read through `system.file()`. The
ordinary pipeline processes a generated chr1:4 G>A synthetic VCF and prints
**1 event, 7 origins, 4 candidates, prediction NOT_RUN, reopen identical**.
The nine reference files are exact copies of committed reconstruction fixtures;
[installed data provenance](../../inst/extdata/quickstart/README.md) records their
hashes. Their assembly/catalog labels support the test contract; their sequences
are hand-authored, not human reference data. The 3-aa/5-aa peptide lengths are
demonstration settings, not production predictor lengths.

The example leaves objects `cfg`, `p`, `events`, `origins`, `results`,
`multisample`, `explanation` and `reopened` in your environment. It prints the
fresh temporary directory where it saved inputs and the run. Keep it during the
session to inspect/reopen; R temporary storage is not durable storage.

## Explore the completed run

See the [public API reference](api-reference.md) for the twelve v0.2.0 functions,
the four additional development setup functions, selectors and runnable examples.
The native workflow requires no
external converter or Registry connection.

## Use R-based setup in the development version

The published development source adds resource downloads/checks and predictor
configuration from R. These functions are absent from v0.2.0 even though the
development package still reports version 0.2.0. Use the pinned source commit
and availability check in the [R setup guide](setup.md#install-the-development-source),
then choose an explicit resource profile and destination. Existing analysis
configurations continue to work without these setup calls.

## Configure your own project

Choose original input, reference and output locations, for example
`~/neo-data/reference/`, `/analysis/sample01/inputs/` and
`/analysis/sample01/runs/`. Absolute paths may be anywhere. Supply your actual
paths; no Codex directory convention is required.

With the development version, [R setup](setup.md) returns reference configuration
fragments directly and does not require a checkout or YAML template. The
following existing file-based route remains available in both versions.

Acquire the primary references and obtain
[reference.primary.example.yaml](../../config/reference.primary.example.yaml)
from the GitHub checkout. Set `reference_config` to its saved filename,
`data_dir` to your acquisition destination, `run_dir` to your output root and
`somatic_vcf` to your VCF filename before running this recipe. Its VCF sample
column must be `sample-1`; replace metadata and HLA declarations with your own.

```r
cfg_user <- neo_config(file = reference_config, env = character(),
  paths = list(data_root = data_dir, run_root = run_dir))
p_user <- neo_project(
  samples = data.frame(sample_id = "sample-1", patient_id = "subject-1"),
  hla = list("subject-1" = "HLA-A*02:01"), config = cfg_user,
  inputs = list(list(input_id = "somatic-1", input_role = "somatic_vcf",
    sample_id = "sample-1", path = somatic_vcf)))
p_user <- neo_run(p_user)
neo_results(p_user)
```

For reproducibility, save the configuration, metadata, HLA and input declarations
outside runs. Retain original provider output if it cannot be recreated. Change
scientific inputs/settings by creating a new project, not mutating a completed
run. Add optional sources using the [input/configuration guide](inputs.md).

For a complete offline recipe with several patients or samples, see the
[four multi-patient/multi-sample examples](multi-sample-examples.md). Several
separately bound single-sample VCFs can share one project/run. This differs from
the [dedicated longitudinal workflow](workflows.md#longitudinal-and-multi-sample-analysis).

## Reopen and inspect

```r
reopened <- neo_project(path = p$path)
stopifnot(identical(neo_results(reopened), neo_results(p)))
reopened <- neo_run(reopened)
stopifnot(identical(neo_results(reopened), neo_results(p)))
```

`neo_project(path = ...)` uses the saved configuration and verifies identity,
input/resource integrity and safe paths. Retained intact output can be reused;
corrupt or changed output is not silently accepted. Run directories can remain
on disk for convenience. If a run was deleted, create a new project from the
original inputs and recipe. See [troubleshooting](troubleshooting.md) for recovery.
