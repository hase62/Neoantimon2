# Set up references and native predictors from R

[Resources and tools](resources-and-tools.md) · [Configuration](inputs.md) ·
[API reference](api-reference.md)

Version 0.3.0 includes these four R setup functions. They are absent from the
immutable v0.2.0 release. No development repository is needed at runtime.
Existing configurations and `neo_run()` work independently.

## Install v0.3.0

Follow the [pinned release installation](getting-started.md#install-and-run-the-offline-example)
or install the [verified source archive](getting-started.md#install-a-source-archive).
Restart R before replacing a loaded package, then verify its identity:

```r
library(Neoantimon2)
stopifnot(as.character(packageVersion("Neoantimon2")) == "0.3.0")
stopifnot(all(c("neo_download_resources", "neo_check_resources",
  "neo_check_predictors", "neo_setup_predictor") %in%
  getNamespaceExports("Neoantimon2")))
```

Earlier development snapshots exposed the setup APIs while reporting version
0.2.0. Record the source commit when reproducing those historical runs; the
immutable v0.2.0 release itself does not contain these APIs.

## Inspect resources before downloading

Choose your own durable data directory. This example respects an optional
environment setting; the directory does not need to follow a developer layout.

```r
library(Neoantimon2)
data_dir <- Sys.getenv("NEOANTIMON2_DATA_ROOT",
  unset = file.path(path.expand("~"), "neo-data"))
report <- neo_check_resources(profile = "minimum", dest = data_dir)
report$files[, c("resource_id", "role", "status", "integrity")]
```

Missing files produce diagnostics. A named file is checked for readability,
format and available integrity/identity information before reuse. The report
distinguishes missing resources, invalid content and missing external tools.

## Download public references

Start with an individual, manageable selection if you do not need a full genome:

```r
mane <- neo_download_resources(resources = "mane", dest = data_dir)
stopifnot(mane$ok)
mane$files[, c("role", "status", "integrity")]
```

`minimum` selects `genome`, `gencode` and `mane`: GRCh38.p14 primary FASTA/FAI,
GENCODE 50 primary annotation and all sequence catalogs, Sec metadata, and MANE
1.5 summary/GTF. It requires several GB of disk space plus staging. These three
selections use R's downloader and streaming decompression; Python is not needed.
Preparing a missing genome FAI requires `samtools` on PATH. A missing tool is
reported before starting that genome download.

For a complete reference setup, use
`refs <- neo_download_resources(profile = "minimum", dest = data_dir)`.
The returned fragments use the existing configuration contract:

```r
# After obtaining the complete minimum profile:
refs <- neo_download_resources(profile = "minimum", dest = data_dir)
cfg <- neo_config(env = character(),
  paths = list(data_root = data_dir, run_root = file.path(data_dir, "runs")),
  scientific = refs$scientific, resources = refs$resources)
check <- neo_check_resources(config = cfg, validate_reference = TRUE)
check$reference_validation
```

The full reference validator checks relationships between catalogs, annotations,
identifiers and declared assemblies/releases. It can take substantial time and
memory. Without `validate_reference = TRUE`, that check is explicitly `NOT_RUN`;
file integrity checks alone do not establish complete analysis compatibility.
An individual selection such as MANE supplies only its roles, not all nine
reference roles needed by the production pipeline.

`profile = "full"` preserves all ten existing selections. Individual IDs also
include `imgthla`, `1000g`, `hcc1395`, `processed-rna`, `longread`, `full-rna` and
`full-truth`. These use the packaged existing Python helper and selection-specific
Git, samtools/bcftools or SRA Toolkit commands. Supply `python = "/your/python3"`
when necessary. The [resource inventory](resources-and-tools.md#public-acquisition-is-separate-from-r-installation)
describes sizes, identities, usage conditions and additional tools. Full acquisition
is optional and may be very large; it is not a prerequisite for ordinary analyses.
Supporting datasets remain separate inputs. Successful acquisition does not make
every dataset compatible with every reference configuration; for example, IGSR
subsets retain their original contig names rather than silently adding `chr`.

## Existing files, checksums and interrupted operations

Valid existing files are reused. Invalid existing files are retained by default;
inspect the diagnostic, choose a separate destination, or explicitly request
`overwrite = TRUE` to replace invalid artifacts after validation. Owned staging
is cleaned on ordinary errors/interrupts. A process killed externally can leave
a staging directory or lock: check that no acquisition is running before dealing
with that owned state. Do not run different acquisition tools concurrently in
one destination.

A header-only `1000g` regional subset left by an earlier extraction is invalid,
even if its filename and index exist. The acquisition helper derives the region's
contig name from the source; it does not rename VCF contigs for analysis. Inspect
the diagnostic and recreate the subset in a separate destination, or explicitly
choose `overwrite = TRUE` after deciding which invalid files may be replaced.

SHA256 states distinguish declared hashes, R acquisition-receipt hashes, and
observed-only hashes. A local receipt binds observed bytes; it is not a publisher
authentication certificate. Reused files without source receipts do not acquire
an invented original download date or proof of origin. Complete gzip checks,
selected release/header checks and FASTA/index consistency complement checksums.
Full scientific reference validation remains a separate explicit check.

One packaged catalog is shared by R and Python. The existing standalone Python
and shell entry points remain available; no additional R setup call is required
to use their outputs in an existing configuration.

## Configure legitimately obtained prediction tools

Inspect the existing native integrations without installing anything:

```r
tools <- neo_check_predictors()
tools$predictors
```

`netmhcpan` supports the existing 4.1/4.2 adapter families; native `netmhciipan`
supports 4.3; `mhcflurry` uses the existing 2.x adapter. A detected version is
`CONFIGURED_UNTESTED`, not a promise that models or prediction execution work.
Missing or unconfigured tools return informative statuses and do not prevent
reference preparation or offline analysis.

After obtaining a vendor installation, pass its executable path. This complete
offline example intentionally demonstrates an unavailable installation:

```r
setup <- neo_setup_predictor("netmhcpan",
  executable = file.path(tempdir(), "not-installed-netmhcpan"))
setup$checks$predictors
cfg_with_predictor <- neo_config(env = character(), predictors = setup$predictors)
```

For your real installation, replace that path. `neo_setup_predictor()` returns
ordinary `predictors` request records; it never modifies an existing config or
vendor installation. `smoke_test = TRUE` explicitly runs a small synthetic query
through the existing native adapter and checks parsed evidence. `SMOKE_PASSED`
is a technical execution result, not validation of prediction accuracy.

Obtain NetMHCpan/NetMHCIIpan through the [official DTU routes](resources-and-tools.md#external-tools-and-provider-inputs),
including required licensing/registration and platform instructions. These tools
are not downloaded, redistributed, unpacked or installed by the setup API.
Indirect-HLA 3.2-family imports remain separate from the native 4.3 adapter.

## MHCflurry environments and models

Configure an existing MHCflurry installation with `executable`, explicit
`model_path`, and the matching `python` interpreter. Python/package/backend checks
do not prove that a model tree is compatible; use an explicit smoke check when
appropriate. Older 2.x installations require TensorFlow even when their package
metadata omits it; newer 2.3 releases use PyTorch. The runtime probe checks the
chosen environment rather than assuming the active Python is the correct one.

Optional installation requires **`install = TRUE`**, a **new** `environment`
directory and an explicit MHCflurry `version`. Additional backend dependencies
must be explicit `package==version` pins in `dependencies`. No existing or active
Python environment is changed. Failed new environments remain marked incomplete;
successful package preparation records a requirements snapshot and still reports
missing models separately. Models are not fetched automatically. Choose Python,
package/backend versions and models using the
[provider's documentation](https://openvax.github.io/mhcflurry/).

Only macOS arm64 with R 4.5.1 was validated in this development cycle. Mock-tool
tests do not establish compatibility with every vendor build or platform.
