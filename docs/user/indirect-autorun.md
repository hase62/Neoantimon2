# Literature-derived indirect-HLA autorun v0.1

This optional **local development** API automates the existing strict molecular
baseline. It is not official PIRCHE, Frost/PIRCHE v4 numerical reproduction,
clinical rejection-risk prediction or immune-response validation. The frozen
`cmpb-code-v1` source has 20 APIs; this branch adds one optional API. No existing
Native Grade, somatic/external score, Yield mathematics or Origin ID changes.

## One command after explicit setup

Supply either the existing `allo-foundation-1.0.0` declaration or a compact list
of `donor` and `recipient` typing. Compact typing requires **two explicit copies**
per locus A/B/C/DRB1/DQB1 and recipient DRA; repeat an allele only for explicitly
known homozygosity. One allele is not silently completed into a genotype.

```r
run <- neo_indirect_autorun(
  allo = list(donor = donor_typing, recipient = recipient_typing),
  presenters = c("HLA-DRA*01:01~HLA-DRB1*01:01"),
  reference = list(path = "hla-reference.json", sha256 = reference_sha256),
  executable = "netMHCIIpan",
  model = list(path = "licensed-netMHCIIpan-3.2", identity = model_release_id),
  authorization = list(status = "authorized_local_installation",
                       terms = "My legitimately installed local provider; applicable vendor terms"),
  output_dir = "indirect-run")
run$result$contexts[[1]]$complete_count
```

The typing variables are named lists such as `A = c("A*01:01", "A*02:01")`;
all required loci must be supplied. Declare all typed recipient DRB1 presenters.
The alpha/beta pair is never inferred. The first autorun profile accepts explicit
DRA*01:01 with two-field DRB1 only. Other DRA alleles, extra-resolution beta names,
DRB3/4/5, DQ/DP presentation, HSCT/reverse directions and PIRCHE-I are not silently
mapped into this profile. The underlying import/foundation workflows remain
available independently with their existing uncertainty semantics.

## Reference setup

Use the existing normalized `hla-reference-1.0.0` manifest with exact allele to
literal precursor-protein mappings. Acquire an appropriate versioned
[IPD-IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/) resource under its terms and
retain source/release, sequence identifiers and normalization provenance. A raw
FASTA or alignment download alone is not a complete reference manifest. Do not
strip gaps, complete partial sequences, select a nearest allele or assert full
coverage without evidence.

The manifest records `reference_source`, `reference_release`, nomenclature
source/release, species, provider/parser versions, `sequence_policy="precursor"`,
`completion_policy="none"`, coordinate convention and per-allele sequence,
coverage and domain intervals. The compact input convention requires a
`full_precursor` interval spanning each supported complete sequence. Keep partial
or unavailable entries explicit. An example of the **format only**, using
synthetic sequences, is installed as
`extdata/indirect-autorun/synthetic-reference.json`; it is not a real IPD dataset.
The original field definitions are in `inst/schema/allo-foundation-1.0.0.json`.

```r
reference_sha256 <- digest::digest(file = "hla-reference.json",
                                   algo = "sha256", serialize = FALSE)
```

Typing is conservatively normalized, then resolved by the existing foundation.
Autorun requires exact reference allele mappings for evaluable calls; it rejects
an otherwise sequence-equivalent but differently named allele substitution.
Full foundation declarations can preserve incomplete/unknown typing. The existing
baseline then returns an unavailable complete count rather than a definitive zero.

## Licensed NetMHCIIpan 3.2 setup

Obtain the executable and matching model data through the official
[DTU NetMHCIIpan 3.2](https://services.healthtech.dtu.dk/services/NetMHCIIpan-3.2/)
registration/download process. Neoantimon2 never downloads vendor binaries,
weights or accepts the license. `authorization` is an explicit user declaration,
not an automated legal-license check.

Supply the physical installed `-rdir` root containing `data/`, a versioned model
identity and optionally its expected tree SHA-256. A computed tree pin is always
saved; a supplied pin is checked. Internal model-tree symlinks are rejected
rather than followed outside the declared tree. Use a physical installation
root, not a directory of unrelated tools. The executable and every model file
are hashed before prediction and checked again afterward.

Actual `-h` output must identify the 3.2 family and advertise `-f`, `-a`,
`-inptype`, `-rdir` plus `-list` or `-listMHC`. The supported allele list is
executed and each declared beta checked exactly (e.g. DRB1_0101). Prediction uses
peptide input (`-inptype 1`), one explicit recipient DR per invocation and sorted
bounded batches (default 256). A supported `-tdir` is used if advertised. BA is
3.2's affinity output; the modern 4.3 `-BA` switch is not blindly passed. Threads
are bounded with OMP/BLAS environment values 1. Each command has a timeout and
separate retained stdout/stderr. Argument values are individually quoted through
base R `system2`; no `sh -c` or user-interpolated command string is used.

The parser accepts the legacy `Seq/Allele/Pos/Core/Affinity(nM)` BA layout and
the `Pos/MHC/Of/Core/Affinity(nM)` BA layout, with optional percentile output.
It requires a 3.2 banner, literal 15mer, exact requested allele, explicit zero-based
core offset 0–6 and a nine-residue forward core matching that ligand register.
No missing or ambiguous core is inferred. Unknown columns, changed versions,
invalid numbers/cores, unexpected queries, conflicting duplicates, missing self
rows and failed/truncated output stop with retained diagnostics. Harmless identical
duplicate evidence is normalized by the existing provider contract. Rank, when
absent, remains NULL. Raw numeric strings are preserved in raw_output.

DTU's 3.2 webpage documents the legacy BA columns but labels its displayed example
as 3.1. The mock tests exercise these explicit layouts; they are not an authentic
3.2 execution receipt. Until a licensed 3.2 engine is tested locally, real vendor
CLI/model integration remains **REAL_PREDICTOR_E2E_NOT_TESTED**. An incompatible
wrapper fails with an actionable setup message rather than falling back to 4.3.

## What is generated and counted

The runner calls existing `allo_foundation`, `indirect_generate` and
`indirect_self_catalog`. It retains donor peptide Occurrences, including repeated
sequences and coordinates, and every self-window mapping. Only exact
`(15mer, recipient presenter)` predictor queries are deduplicated. A peptide shared
by donor and self is queried once but retains both roles. All eligible recipient
self 15mers are included, not just donor ligands.

The output is converted to the unchanged `indirect-binding-import-1.0.0`, then
passed through existing `indirect_config`, `indirect_baseline` and public reporting.
The frozen strict policy uses **affinity <1000 nM** (exactly 1000 is excluded),
same-presenter designated primary 9mer self comparison, and a unique
`(core, DR presenter)` union. Missing/incomplete self evidence cannot establish
nonself or a definitive zero. An inspectable conditional result is separate from
the complete count; binding is not presentation or clinical risk.

## Dry run, repeat and failure handling

```r
plan <- neo_indirect_autorun(allo, presenters, reference,
                            output_dir = "indirect-plan", dry_run = TRUE)
# No executable, authorization or models are needed; no provider command runs.

same <- neo_indirect_autorun(allo, presenters, reference, executable = tool,
  model = model, authorization = authorization,
  output_dir = "indirect-run", resume = TRUE)
```

Use a fresh directory initially. Completed runs/plans may be reopened only with
`resume=TRUE`, identical input/reference/tool/model/code identity and verified
saved artifact/checksum files. Ordinary repeated calls refuse an existing target;
they never overwrite it or repeat inference. Changing ordering alone does not
change semantic query identities. Changing a pin or artifact is rejected.
Failed/interrupted runs retain logs and a FAILED manifest. Partial-run resumption
is deliberately not supported in v0.1; choose a new directory after diagnosis.
No failed query is reported as a zero-epitope analysis.

## Outputs

- `run.json`, `SHA256SUMS.txt`: request/tool/data/code identity, process state and
  all persisted artifact hashes.
- `allo-config.json`, `reference.json`, `foundation.json`, `plan.json`: frozen
  supplied inputs and existing foundation/generation semantics.
- `queries.tsv`, `query-origin-links.tsv`, `donor-occurrences.tsv`, `self-windows.tsv`,
  `generation-coverage.tsv`: all molecular and self links, independent of query dedup.
- `raw/`, `execution/`: exact argv, bounded environment, exit status, stdout/stderr,
  executable identity, model-tree and raw-result hashes. Vendor assets are not copied.
- `binding.json`, `indirect-config.json`: normalized evidence and existing import config.
- `query-coverage.tsv`, `summary.tsv`, `result.json`, `explanation.json`, `result.txt`:
  imported status, complete/conditional counts and existing evidence graph.

List-valued TSV cells are JSON. Scientific result references are portable relative
artifact links. Execution argv and original reference bindings may retain local
installation locations as historical metadata; inspect them before sharing results.
The JSON/RDS-like in-memory list can be reopened through this API's completed-run
resume without executing the provider again. This is a standalone workflow, not
an unsupported no-event somatic `neo_run()` invocation.

## Installed synthetic end-to-end demonstration

```sh
Rscript --vanilla -e 'source(system.file("examples", "indirect-autorun.R", package="Neoantimon2"))'
```

It uses only installed toy reference/typing and a clearly labeled Python 3 mock,
creates temporary files and prints a five-unit synthetic molecular result. The
fixture is not NetMHCIIpan, an IPD dataset, an authentic predictor score or a clinical
example. `execution_kind="MOCK_EXECUTION"` and `real_predictor_executed=FALSE`
are explicit. The existing synthetic-provider schema's `execution="not_executed"`
means no genuine predictor execution; metadata separately records the actual mock
process. Both source and installed package tests cover the orchestration.


## Practical use and scope checklist

1. Install this local development source (`R CMD INSTALL /path/to/Neoantimon2`)
   with its DESCRIPTION dependencies. Frozen `cmpb-code-v1` does not include this API.
2. Supply known donor/recipient HLA copies, every typed recipient DR presenter,
   and a complete versioned normalized precursor manifest; compute its SHA-256.
3. Obtain and validate the legitimate **3.2** binary and data as described in
   [NetMHCIIpan setup](netmhciipan32-setup.md). Set the setup variables used above.
4. Call `neo_indirect_autorun()` once. The API writes the donor/self query files,
   invokes the tool, parses raw BA/core output, creates `binding.json`, performs
   self comparison and writes the final result and checksums automatically.
   **No manual predictor-result creation, score editing or conversion is required.**
5. Inspect `run$status`, `run$result$contexts[[1]]$completeness_status`,
   `complete_count`, coverage and reason codes. Operational COMPLETED with an
   unavailable complete_count is not a complete scientific estimate.

The installed synthetic example shown above is runnable without a licensed tool;
it requires Python 3 and clearly reports mock-only execution. Production requires
real HLA sequences and a working licensed tool. Do not reuse toy reference
proteins as biological HLA resources.

This profile supports **SOT donor-to-recipient, DRA*01:01–two-field DRB1** and
NetMHCIIpan **3.2 BA** only. Official PIRCHE v4/Frost has different prediction
semantics; numerical equivalence and clinical performance are unverified.
[Scope and future candidates](indirect-autorun-scope.md) lists the remaining areas;
they are not silently enabled by a successful run. Self coverage concerns the
explicit HLA precursor universe, not the entire human proteome.
