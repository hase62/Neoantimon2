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

## Snapshot saved-result reproduction

Choose this route to check the published **saved tables**, rather than to run a
new genomic analysis. Open the [reproducibility-v1 GitHub pre-release](https://github.com/hase62/Neoantimon2/releases/tag/reproducibility-v1)
and download the named [reproduction ZIP](https://github.com/hase62/Neoantimon2/releases/download/reproducibility-v1/Neoantimon2_Reproducibility_v1.zip)
and [SHA-256 file](https://github.com/hase62/Neoantimon2/releases/download/reproducibility-v1/Neoantimon2_Reproducibility_v1.zip.sha256).
The automatic GitHub **Source code** ZIP alone does not include the eight result
tables. This 20-API snapshot is distinct from 22-API main.

With Python 3 and its standard library, from the directory holding the download:

```sh
# Compare with the downloaded .sha256 file (macOS):
shasum -a 256 Neoantimon2_Reproducibility_v1.zip
unzip Neoantimon2_Reproducibility_v1.zip -d snapshot-reproduction
cd snapshot-reproduction
python3 -I scripts/verify_package.py
python3 -I scripts/verify_snapshot.py
```

Use a fresh extraction directory; keep Python assertions enabled (no `-O`). The
asset SHA-256 is `f4fb8eb9dff5af88b6cb73fc19aa582ef415753f45df3851cc60d555c527ef23`.
The package's own README/REPRODUCE describe every command and its scope.
Expected saved-table totals are **15 unique mutant peptides, 75 peptide–HLA
candidates, 36 Origins, 180 links and 300 external prediction measurements**.
The AGRN `AAPCRVLPQG` example has 5 Origins, 5 HLA candidates, 25 links and
20 external measurements. Native status rows are not external measurements.

The ZIP contains all eight input TSVs, frozen expected outputs, exact source
archives and scoped licenses. It needs no R, genome reference or model weights
for this read-only reaggregation. It does **not** rerun molecular reconstruction,
predictors, oracle cases, pVACseq or the full 41,072-record HCC1395 graph. Historical
qualifications remain unchanged, including the original driver/controller exit.
The original 15-sequence selector executable remains unrecovered; the saved
membership/reaggregation must not be presented as a complete VCF-to-pilot replay.
Independent scripts are MIT and author-owned explanatory materials CC BY 4.0;
Core/Converter remain PolyForm Noncommercial and scientific/third-party data
retain their own terms. See the package's LICENSES/SCOPE.md and NOTICE.md.

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
R CMD check --no-manual --install-tests Neoantimon2_0.3.0.9000.tar.gz
python3 scripts/development/test_public_data_acquisition.py
```

Licensed predictors, real reference bundles, clinical outcomes, and other
operating systems require separate validation. The standalone Community Converter
Registry is maintained at https://github.com/hase62/Neoantimon2-Converters and
served at https://www.hase62.jp/Neoantimon2-Converters/. Its infrastructure is not
part of this package; the Registry client and synthetic integration fixtures are.

## Development version and integration checks

The default branch uses 0.3.0.9000 as an explicit R development marker; the
published v0.3.0 tag, source archive and release assets are immutable. This does
not change scientific schema or identity versions. Record the Git revision as
well as DESCRIPTION because earlier unpublished snapshots still reported 0.3.0.
The sixteen released APIs remain compatible. The fixed snapshot adds four APIs
(20 total); main adds both indirect-HLA runners (22 total). The complete
[API reference](api-reference.md) and [version choices](../../README.md#choose-a-version)
distinguish the three editions.

After installing this source snapshot into an isolated library, run focused
`testthat::test_local()` filters `experimental-yield`, `result-reporting` and
relevant provenance tests, then the full source suite. Build and run
`R CMD check --no-manual --install-tests` for the full installed suite. The
repository-only [AI acceptance script](../../scripts/development/validate_ai_playbook.R)
checks the preserved baseline helper on the release or explicitly supported
current development version. [Historical evidence](../development/validation-evidence.md)
is distinct from a fresh integration test receipt.

The public R modules, tests/fixtures, man pages, runtime resources and installed
examples retain their conventional paths. Supplementary guides remain in docs;
there is no duplicate vignette hierarchy. Figure SVG/PNG/PDF/drawio files are
intentional documented source/exports, not disposable scientific outputs.
Generated caches, local libraries, raw datasets and internal audit reports stay
outside Git; `.gitignore` and `.Rbuildignore` already exclude their usual locations.
