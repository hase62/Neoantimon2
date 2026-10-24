# Reproducibility and source layout

[Documentation index](../README.md) · [Troubleshooting](troubleshooting.md)

Keep source code, original inputs, and generated runs separate. Users can choose
any locations; all paths are explicit. The maintainer convention is
`repos/Neoantimon2`, `data/Neoantimon2`, and `runs/Neoantimon2`.

Save the package version and source commit, configuration, sample/subject metadata,
HLA declarations, reference versions/checksums, original inputs, and provider
versions/outputs. Preserve manual curation and irreplaceable predictor output
outside generated runs. Licensing and access restrictions still apply to data
and third-party tools; the package does not redistribute these resources.

A completed run contains manifests, checkpoints, intermediates, caches, results,
and logs. Reopen with `neo_project(path = saved_run_path)` to verify saved identities
and integrity. Changed scientific inputs or configuration require a new project.
A deleted run cannot be reopened; regenerate it from the preserved recipe.
Different output locations and timestamps can change manifests without changing
biological identities. Compare scientific outputs separately from these fields.

## Public source map

| Location | Purpose |
| --- | --- |
| `R/`, `man/` | Package source and installed API help; R source remains flat |
| `inst/schema/` | Versioned scientific and interoperability contracts |
| `inst/extdata/`, `inst/examples/` | Authored synthetic examples, metadata, and installed quickstart |
| `tests/` | Regression tests and required synthetic fixtures |
| `docs/user/` | User documentation |
| `config/` | Current project and primary-reference configuration examples |
| `scripts/data/` | Public-resource acquisition and reproduction recipes, without downloaded data |
| `scripts/development/` | Public synthetic external quickstart and offline acquisition regression |

The two public scripts retain their paths because the regression harness and
user recipes reference them. The directory name does not imply that private
execution logs are distributed. Installed packages omit repository-only scripts
and guides; installed help and the offline quickstart remain available.
Essential historical scientific specifications are preserved separately from
public source. Frozen metadata may retain historical audit locators; these are
provenance fields, not runtime file dependencies. Internal audits, obsolete design
bundles, original research data, and old development Git history are excluded.

## Reproduce the synthetic external workflow

From the source checkout, use the [external quickstart](external-predictions.md#executable-offline-example).
It creates synthetic inputs outside runs and refuses to overwrite differing
existing inputs. The installed native quickstart uses temporary storage; copy its
recipe and inputs to durable storage if you need to retain the run.

For source regression, install the suggested dependencies and run:

```r
testthat::test_local(".", reporter = "summary")
```

```sh
R CMD build .
R CMD check --no-manual Neoantimon2_0.3.0.tar.gz
python3 scripts/development/test_public_data_acquisition.py
```

Licensed predictors, real reference bundles, clinical outcomes, and other
operating systems require separate validation. The standalone Community Converter
Registry is maintained at https://github.com/hase62/Neoantimon2-Converters and
served at https://www.hase62.jp/Neoantimon2-Converters/. Its infrastructure is not
part of this package; the Registry client and synthetic integration fixtures are.
