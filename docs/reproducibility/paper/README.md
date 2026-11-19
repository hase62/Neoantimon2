# Neoantimon2 paper replay

This page provides the frozen inputs, independent expectations, and execution
and verification scripts for three bounded studies. Use the exact commit
cited in the paper. The unchanged historical
[eight-table snapshot](https://github.com/hase62/Neoantimon2/releases/tag/reproducibility-v1)
has a separate scope.

| Study key | Scope and recorded outcome |
| --- | --- |
| `mr11` | Eleven synthetic conditions. Nine sequence expectations, one unresolved-contract case and one intentional invalid-annotation rejection. Ten zero exits and one expected exit 1, 2,500 checks and four separate-process readers. |
| `study6` | Six additional synthetic conditions, including conditional Origin support and phase or boundary behavior. Six successful cases, separate-process reopening and 946 checks. |
| `r09` | Controlled post-prediction comparison using native pVACtools report/filter commands and external queries. The saved comparison had 45 initial and 49 fresh-process checks. Predictor-like scores are controlled input. This scope does not run pVACseq reconstruction or model inference. |

The history of the original 18-condition study and the full HCC1395
analysis are separate from these three runnable scopes. The counts above
refer to recorded validation. A new run produces its own written evidence.

## Install prerequisites

Use Python 3.12 or newer and R 4.1 or newer with `jsonlite`, `digest`,
`stringi`, `yaml`, and `xml2` in R's normal library path. Retain at least
25 GiB of free disk space. The frozen drivers use isolated R libraries, a
6 GiB R vector heap limit and single-threaded numerical-library settings.
Python 3.11 can run the read-only hash check. R09 additionally requires a
separately installed pVACtools 7.1.5 interpreter with the exact
[969-file source/resource inventory](payload/r09/r09/PVACTOOLS_EXPECTED_SOURCE.json).
The recorded R09 environment used Python 3.10.19 for pVACtools, with
numpy 1.26.4, pandas 1.5.3, scipy 1.15.3 and pysam 0.24.1. Use a separate
Python 3.12+ interpreter for the preparation and replay drivers.
Preparation checks the reader's environment and uses local dependencies.
No predictor assets or pVACtools runtime are bundled.

The shared Core archive is the restored 502-file **0.3.0.9000, 20 API** source,
SHA-256 `af58b053bc57fcb1f46fb8378ee16173621579e4aec156cdb54767d26792486a`.
The 22-API development source in the repository is a separate implementation
scope. Each run installs the verified 20-API source in its own library.
[The payload lock](payload-lock.json) identifies the three donor candidates and
all retained file bytes. Preparation recreates packaging manifests in the run
directory and copies the one shared source archive under the original paths
used by the drivers. Input, expected-value and execution file bytes stay fixed.
R09's [protocol](payload/r09/r09/R09_CANDIDATE/ROUTE_AND_SCOPE.md)
records the original study design. Its historical Core digest identifies the
original execution material. The active freeze and shared-source SHA-256
above identify the runnable replay.

## Clone and verify

In a POSIX-compatible shell, replace `PAPER_COMMIT` with the full commit SHA
cited in the paper.

```sh
git clone https://github.com/hase62/Neoantimon2.git
cd Neoantimon2
git checkout --detach PAPER_COMMIT
paper="$PWD/docs/reproducibility/paper"
run_root="$PWD/runs/paper-replay"
mkdir -p "$run_root"
python3 -I "$paper/prepare.py" mr11 --verify-only
```

The read-only command verifies the payload file set and SHA-256 values for
all three scopes, the fixed Core archive, all 502 source entries, its version
and its 20-export interface. Choose one route below. Use a new run directory
for every invocation. Set `run_root` to another absolute local `runs` path
if preferred.

## Verify the reported numbers

The [released snapshot](https://github.com/hase62/Neoantimon2/releases/tag/reproducibility-v1)
is already public. This page adds a separate route for reaggregating observations
cited by the paper. After checking out the full paper commit as above, use a
fresh directory under `runs` and Python 3.12 or newer. No R packages, model
runtime, project files, completed manuscript tables or private paths are needed
for the number check.

```sh
verification="$run_root/number-check-1"
mkdir -p "$verification"
asset="$verification/Neoantimon2_Reproducibility_v1.zip"
curl --fail --location --silent --show-error \
  https://github.com/hase62/Neoantimon2/releases/download/reproducibility-v1/Neoantimon2_Reproducibility_v1.zip \
  --output "$asset"
python3 -I "$paper/check_numbers.py" --snapshot "$asset" > "$verification/numbers.json"
```

The checker verifies the release SHA-256
`f4fb8eb9dff5af88b6cb73fc19aa582ef415753f45df3851cc60d555c527ef23`
and [every included observation's size and SHA-256](observations-sha256.json).
It reads the eight-table archive where needed and prints JSON. Review the
top-level keys listed below. The [source ledger](observations/source-ledger.json)
lists descriptive source labels, file names, byte counts and SHA-256 values.
The compact records carry the fields used by the comparisons. Original private
workspaces, review archives, bulky references, trained weights and vendor
runtimes are outside this page.

| JSON key and paper scope | Input records, filters and calculation |
| --- | --- |
| `oracle` (SM2, S3, S6–S8) | `observations/oracle/cases.json`, all 54 case IDs from the three archived fixture sets at Core `491d72f746399f5fc8aa88bb9f985ef2ce477e34`. Compare `expected` with `observed` by protein, spliced DNA, peptide length/start/sequence where applicable, and available WT correspondence. Count the three qualified cases separately and count empty expected window sets. The separate public pipeline records count distinct `peptide_sequence`, `origin_id` and candidate-to-Origin links. The low-level germline comparison uses sequence set differences. The synthetic isoform check translates the saved exons and compares the three retained source links. |
| `pvacseq` (SM3, S4–S5, S9) | `observations/pvac/`, all 24 selected Event rows, 168 versioned Native and exact VEP-matched contexts, 445 unfiltered VEP rows, the pVACtools 7.1.5 converted rows and saved FASTA. Join allele coordinates and versioned transcript, with the recorded SNV/deletion coordinate convention. Compare literal VEP 116.3 `FrameshiftSequence` proteins or apply the saved `Protein_position`/`Amino_acids` edit to `WildtypeProtein`. Separately join full converted `index` to `MT.` FASTA headers and test fragment containment. The missing PDE3B fragment and the four unconverted contexts stay separate. These are checks of archived sequence views. They do not execute VEP or pVACseq. |
| `full_graph` (Table 4, SM1) | `observations/full-graph/`, all Event keys (`input_id`, `line`, `selected_alt_index`), Origin IDs, candidate IDs and their indexed endpoints in the original audit record order. The compressed link fields are zero-based index lists into the corresponding data rows, excluding headers. Check unique keys and endpoint bounds, sum links and read the saved process exit records. The graph came from the HCC1395 run at Core `11c68160609fc08f7303e07425f5cebfaf6af4a9`. The R driver exited -15 and its controller exited 143. The saved graph was separately verified. |
| `pilot` (Table 4, SM4, S10) | `observations/pilot/`, all five Event IDs, 3,696 Origin rows, their sequence and source Event fields, all candidate peptide–HLA keys and the full link table from `run-c1cbdd292231765567abf64b57c44bfc`. Count distinct sequences over Origins, split lengths 8–12 and 15, join exact link keys, then select the three lexical first sequences per length 8–12. `source-vcf-rows.tsv` identifies the five matching SEQC2 v1.2.1 source lines by file and one-based line number (headers included). Pinned upstream URLs and SHA-256 values are in `source-vcf-inputs.json`. The pilot was selected directly from the two source VCFs. Its original final 15-sequence selector was not recovered. |
| `snapshot` (Figure 1, Table 4, SM4–SM6, S1–S2) | Public ZIP `inputs/report/{events,origins,pmhc_candidates,candidate_origin_links,origin_event_links,pmhc_predictions}.tsv`. Compare natural variant keys to the five source rows and peptide/HLA membership to the pilot. Count IDs, source-link roles and records of kind `measurement` with source `External`. For AGRN, filter sequence `AAPCRVLPQG` before counting links, Events and measurements, and read the exact versioned transcript IDs. `observations/authentic/` retains original MHCflurry 2.0.6 BA and BigMHC EL/IM raw outputs and normalized records. The checker follows SHA-256/row/column locators and verifies raw score strings against converted and saved imported values. Best-HLA score per sequence within one provider, task, metric, version and run gives the sequence ranking. Use ascending affinity for BA and descending score for EL/IM. Break HLA ties by `mhc_key` then `candidate_id`, and rank sequence ties by sequence then HLA. Intersect top-5 and top-10 sequence sets and divide by the union for Jaccard. |
| `saved_replays.study18*`, `integrated18` (Tables 2–3, SM7–SM12, SM17–SM18) | `observations/saved/study18-*` provides all 36 case-version execution/check records, the original support statuses and separate reader evidence. Count distinct condition IDs, process exits and per-case check objects. Read the four original graph counts from their per-execution `counts` object and compare patched support status and selected counts. Count raw import and omitted-query check labels separately. The later integration ledger reuses 18 baseline records, so its counts must not be added to the original run. The complete old independent inputs, drivers and full output projects are unavailable here. |
| `saved_replays.mr11`, `study6` (SM13, SM19–SM20, S15–S17, S22) | Read archived per-check rows and exit/reader records, and compare MR11 exit codes to `payload/mr11/expected/independent_expected.json`. Count positive sequence, unresolved-contract and rejected inputs separately. `observations/study6-graph/` has both complete matched graph projections. In `observations/saved/study6-*`, group decisions by case and candidate, and join per-Origin support. Check observed boundary protein, first stop and Origin windows. |
| `s18` (SM14, S18) | All 36 candidate–Origin rows in [Neoantimon2_S18_observations.tsv](observations/saved/Neoantimon2_S18_observations.tsv). Recompute single-Origin eligibility from Event-level ALT and phase, then compare with the saved full-Origin status. See the procedure below. |
| `saved_replays.r09-*` (SM15, S20–S21) | `observations/saved/r09-*` contains nine command exits, the native binding/aggregate row projections, all 24 fresh normalized context–HLA relations across four output routes, and first/fresh per-query observations. Test agreement by exact natural keys including peptide, transcript, register identifier, HLA and score. Recompute the 2/4/0 relation counts after external Event exclusions and compare first and fresh answers. The two native aggregate modes are separate routes. These scores are declared controls. No Neoantimon2 or pVACtools native eligibility change is inferred from the external query. |
| Other `saved_replays` fields (SM7, SM17–SM19) | Inventory and count the 199 archive manifest entries and match before/after SHA-256 rows for 175 integrated-distribution, 55 restored-study and 429 development files. Count or quote the package-test receipts and R CMD check statuses at their respective software versions. This confirms the saved receipts, not a fresh test-suite execution or the availability of the complete historic distribution archive. |

The three runnable routes below are separate from this saved-record check.
Their inputs, independently frozen expectations and driver scripts are included.
Other scopes use archived observations. In particular, this page cannot repeat
the old 18-condition execution, the full HCC1395 analysis, VEP/pVACseq upstream
reconstruction, the original selector execution or trained-model inference.

### S18 single-Origin queries

The denominator is **6 conditions × 2 HLA candidates × 3 Origins = 36
queries**, comprising 24 queries from the four earlier conditions and
12 from the two additional conditions. Group the supplied saved observations
by `case_id` and `HLA`. Within each group, retain one `origin_id` at a time
for the same candidate. Each row records `origin_dna_evidence` for its required
Events and `origin_phase_condition` for the Origin.

1. If any required Event has tumor ALT count below 5 or normal ALT count
   above 1, mark the Origin ineligible (`ineligible`).
2. If every required Event passes the DNA criteria and the recorded Origin
   phase is conditional, mark it pending evaluation (`not_evaluated`).
3. Otherwise, with complete evidence and resolved phase, mark it eligible
   (`eligible`).

Use the single retained Origin's status as the query result. For the
full-Origin comparison in each group, any eligible Origin gives an eligible
candidate, all Origins ineligible gives an ineligible candidate, and the
remaining mixtures give `not_evaluated`. Compare the recomputed statuses to
the saved observations, and cross-check all 12 additional-study rows against
the earlier query evidence. Native selection is read from the saved full-Origin
records. These are read-only queries over archived observations.

To check S18 without downloading the public snapshot, run
`python3 -I "$paper/check_numbers.py"` and inspect the `s18` output. The
checker reports 36 queries over 12 candidate-condition groups, with a status
different from the full-Origin status in 14 queries and agreement in 22.

## Execute the bounded replays

### MR11

```sh
run="$run_root/mr11-1"
python3 -I "$paper/prepare.py" mr11 --runs-dir "$run"
python3 -I "$run/candidate/run_replay.py" --out "$run/output"
```

Verify the per-case exits and independent checker status in
`output/FINAL_EXECUTION_RECEIPT.json`, the separate reader in
`output/SEPARATE_READER.json`, and the case and comparison tables in
`output/verification/`. Reserved invalid-annotation case S03 should exit 1
inside an overall successful driver run. The unresolved case is a contract
check, outside the nine sequence expectations.

### Additional six conditions

```sh
run="$run_root/study6-1"
python3 -I "$paper/prepare.py" study6 --runs-dir "$run"
python3 -I "$run/candidate/study6/scripts/execute.py" --work "$run/output"
```

Preparation runs the original byte and freeze-binding verifier. The driver
runs the separate R reader and independent saved-result checker. Inspect
`output/EXECUTION_RECEIPT.json` and
`output/verification/INDEPENDENT_REVIEW.json` for case exits, reopen status
and the independent comparisons.

### R09

Set `PVAC_PYTHON` to the absolute path to the pVACtools 7.1.5 interpreter
in the environment you installed.

```sh
run="$run_root/r09-1"
python3 -I "$paper/prepare.py" r09 --runs-dir "$run" --pvac-python "$PVAC_PYTHON"
candidate="$run/candidate/r09/R09_CANDIDATE"
python3 -I "$candidate/scripts/run_r09.py" --work "$run/output" --pvac-python "$PVAC_PYTHON"
python3 -I "$candidate/scripts/query_saved.py" --root "$candidate" --execution "$run/output" --out "$run/queries-first"
PYTHONPATH="$run/candidate/r09/pvac-isolated" PYTHONDONTWRITEBYTECODE=1 \
  "$PVAC_PYTHON" -B "$candidate/scripts/query_saved.py" --root "$candidate" \
  --execution "$run/output" --out "$run/queries-fresh" --reread \
  --baseline "$run/queries-first/FIRST_QUERY_RECEIPT.json"
```

Preparation checks the original freeze and copies the verified local pVACtools
package into the run directory. Check `output/EXECUTION_RECEIPT.json` and the
two query receipts under `queries-first/` and `queries-fresh/`. Examine the
native binding-filter and aggregation outputs separately from the controlled
input table. The fresh query and separate R reader test retrieval after a new
process starts.

## Licenses and attribution

See [license scope and sources](LICENSE-SCOPE.md) for the applicable notices
and terms. The archived Core keeps its own license and notice. Third-party
tools, model assets, and original private audit records are obtained separately
where applicable.
