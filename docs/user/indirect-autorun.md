# Indirect-HLA autorun: NetMHCIIpan 4.3k standard, explicit 3.2 legacy

`neo_indirect_autorun()` executes the configured local predictor, imports its raw
output and reports molecular counts in one call. Users do not manually create,
edit or convert prediction result files. This is research software, not official
PIRCHE/Frost or a clinical risk model. Version alone is insufficient: record the
source commit, exact predictor/model identity and analysis profile.

## One command after setup

```r
library(Neoantimon2)
run <- neo_indirect_autorun(
  allo = list(donor = donor_typing, recipient = recipient_typing),
  presenters = c("HLA-DRA*01:01~HLA-DRB1*01:01"),
  reference = list(path = "hla-reference.json", sha256 = reference_sha256),
  executable = file.path(netmhciipan_root, "netMHCIIpan"),
  model = list(path = netmhciipan_root, identity = model_release_id),
  authorization = list(status = "authorized_local_installation",
                       terms = accepted_local_license_description),
  predictor_version = "4.3k", output_dir = "indirect-run")
run$result$contexts[[1]]$complete_count
```

The variables are explicit setup values, not inferred patient data.
`donor_typing` and `recipient_typing` are named lists: A, B, C, DRB1, DQB1, plus
recipient DRA. Each locus has **two explicit character allele copies** (e.g.
`A=c("A*01:01","A*02:01")`). Repeat only for known homozygosity. Declare every
typed recipient DR presenter as an explicit alpha/beta pair. The existing
`allo-foundation-1.0.0` configuration is also accepted.

The scope remains human SOT donor-to-recipient, source/self A/B/C/DRB1/DQB1,
15mer ligands and DRA*01:01 with DRB1 names having two two-digit fields. DQ/DP,
DRB3/4/5 presentation, other alpha chains and HSCT/reverse directions are not
silently enabled. The internal parser's DP documentation example is not an
advertisement of DP autorun support.

## Reference and vendor setup

Supply a versioned normalized `hla-reference-1.0.0` manifest with exact allele
mappings, literal precursor sequences, release/source/nomenclature provenance,
coverage/domain intervals and `completion_policy="none"`. No nearest allele,
partial-sequence completion or genotype imputation is performed. Raw FASTA alone
is not a complete manifest. Acquire appropriate [IPD-IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/)
resources under their terms. Compute the reference file SHA-256 explicitly:

```r
reference_sha256 <- digest::digest(file="hla-reference.json",
                                   algo="sha256", serialize=FALSE)
```

Use the legitimately acquired [DTU 4.3 distribution](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/)
for the host architecture (this verification uses **4.3k Darwin_arm64**).
Follow its README's site settings for the launcher NMHOME/TMPDIR. Neoantimon2
never obtains a license or downloads vendor assets. Supply the **full installation
root**, with physical `data/` and its host-platform directory. The 4.3k adapter
automatically selects `Darwin_arm64`/`Darwin_x86_64` (or matching installed host
platform) when a platform `bin/` is present; `model$runtime_subdir` can state it
explicitly. The resolved runtime path is pinned. Internal data links are recorded
without traversal; their physical files are hashed once. Any link escaping the
installation root is rejected. The legacy model validator remains strict about
symlinks.

Executables and the complete supplied model tree are hashed before/after use.
An expected `model$sha256` can be supplied. Actual help and allele-list execution
must establish the exact release and requested DR support. 4.3k help does not
print a version banner; the executed list and prediction outputs do. Unpatched
4.3, 4.3j or other releases are not silently accepted as 4.3k.

## BA and EL stay separate

The 4.3k adapter always explicitly requests **`-BA`** for this BA analysis.
Plain peptide input uses `-inptype 1` with context encoding off (actual help's
0 default is checked). `inversion="vendor_default"` keeps the vendor default
(DP-only inversion); `inversion="all"` explicitly adds `-inv_all` for 4.3k.
Filtering/sorting/unique-core options are not enabled. Source/self query coverage
is checked in full. Batches and OMP/BLAS threads are bounded; every argv, exit,
stdout and stderr is retained. 4.3k retains its official TMPDIR-based temporary-name template, recorded in
raw stdout. An output-derived long `-tdir` caused a native buffer overflow in
the Darwin arm64 verification; no such override is passed.
A zero process exit is insufficient: errors/no rows, missing columns/queries,
conflicting duplicates, invalid numbers/cores and wrong alleles fail closed.

Normalized affinity and percentile_rank mean **BA nM and %Rank_BA**. Raw
`Score_BA`, `Score_EL`, `%Rank_EL`, `Core_Rel`, `Inverted`, peptide/MHC/Core/Of
are retained as original strings and exported in `model-metrics.tsv`. EL-only
output cannot complete this BA analysis. EL, binding affinity and experimental
positivity are different quantities. Neither rank nor Score_EL becomes nM.

## Explicit model-specific counting policy

4.3k uses `neo-netmhciipan43k-indirect-hla-ii-dr-ba1000-orientation-v1` and its
own versioned evidence/result schemas. Its **affinity <1000 nM** cutoff is an
explicit exploratory assumption carried from a familiar numerical cutoff; it
is not cross-model calibration or numerical compatibility with classical
PIRCHE-II v3. Exactly 1000 is excluded. Official PIRCHE v4/Frost is different.

A reported inverted core is retained in **model binding direction**. When
`Inverted=1`, Of indexes the ligand **after reversal**. The forward physical
zero-based start is `15 - 9 - Of`; for forward cores it is Of. The reported
nine residues must match that orientation/register exactly. No core is guessed.
The binding reliability is retained, without imposing an invented reliability
threshold or interpreting it as an immune-response probability.

For this new profile, self matches require exact **core sequence + orientation +
presenter**, and the count is the union of those triples. Forward and inverted
binding contexts are not silently collapsed. This is an explicit conservative
engineering policy, not a validated TCR-recognition rule. Each qualifying self
15mer is predicted against the same DR. Incomplete source/self/presenter evidence
makes complete_count unavailable, never a falsely definitive zero. Self means
the supplied HLA precursor universe, not the whole human proteome.

Only exact `(15mer, presenter)` predictor queries are deduplicated. All donor
Occurrences, allele/source coordinates and self-window links remain. Molecular
generation policy/IDs are reused independently of model choice. Selecting a new
analysis profile does not replace Event/Origin semantics, Native Grade/rank or
Experimental Yield calculations.

## Legacy compatibility and reopening

4.3k is the new default. For the unchanged 3.2 adapter, explicitly add
`predictor_version="3.2"`. Its profile, strictly <1000 nM rule, forward-only
primary-core self policy and unique `(core,presenter)` count are preserved.
`predictor_version="auto"` opts into executed 3.2/4.3k detection, with the actual
profile/version recorded; it does not relabel either model as the other.

```r
plan <- neo_indirect_autorun(allo, presenters, reference,
  predictor_version="4.3k", output_dir="indirect-plan", dry_run=TRUE)
# No predictor/model/authorization is required for this plan.
# Reopen a completed run with the identical full argument list plus resume=TRUE.
```

Existing target directories are never overwritten. Completed-run/plan resume
requires identical code, inputs, tool/model and all artifact hashes, and never
re-executes inference. Failed/interrupted batches are retained but cannot resume
in place. Previously frozen 113f801 runs retain their old code fingerprint; use
that revision for their exact resume, or start a fresh directory under this update.
Do not rewrite old manifests to bypass an identity mismatch.

## Outputs and executable examples

`run.json` and SHA256SUMS record resolved model/profile, actual BA settings,
cutoff and counting policy. `queries.tsv`, `query-origin-links.tsv`,
`donor-occurrences.tsv`, `self-windows.tsv` and coverage tables preserve all paths.
`raw/` and `execution/` retain original evidence and exact argv/tool/model hashes.
`binding.json` and `indirect-config.json` declare the selected versioned contract;
`summary.tsv`, `result.json`, `explanation.json` and `result.txt` expose the result.
List-valued TSV cells are JSON. Physical argv metadata may contain local paths;
review before sharing. Vendor assets are never copied into reporting output.

After installing this local source and dependencies (including Python 3 for mocks):

```sh
Rscript --vanilla -e 'source(system.file("examples","indirect-autorun43.R",package="Neoantimon2"))'
# The old indirect-autorun.R example remains an explicitly selected 3.2 mock.
```

The mock is a controlled software fixture, not an actual predictor or real HLA
protein dataset. Its five units are an artificial golden expectation. Genuine
model execution on the existing synthetic HLA fixture is a separate technical
E2E; its result is not forced to match a mock. Neither is clinical validation.
See [scope and future candidates](indirect-autorun-scope.md) and the
[legacy 3.2 acquisition guide](netmhciipan32-setup.md).
