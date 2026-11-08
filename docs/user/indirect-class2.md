# Indirect HLA / PIRCHE-inspired class II research workflow

This development branch adds `neo_indirect_class2()` without altering the frozen
`neo_indirect_autorun()` 3.2/DRB1 contracts, somatic reconstruction, Origin
identities, Native Grade/rank or Experimental Yield. It is not official PIRCHE,
Frost, PIRCHE-I, B-cell matching or a clinical risk model.

## Install this Draft development branch

```r
remotes::install_github("hase62/Neoantimon2", ref="feat/indirect-hla-expanded")
```

For an exact scientific reproduction, replace ref with the reviewed commit SHA
reported by the Draft PR. Package version alone is insufficient: verify that
`neo_indirect_class2` is exported. The default main and historical tags are not
changed by this Draft PR. No predictor or IPD data is installed by this command.

## One command after typing and licensed installation

```r
library(Neoantimon2)
# Fictional research declarations, not patient observations or inferred haplotypes.
donor <- list(A=rep("A*02:01",2), B=rep("B*07:02",2), C=rep("C*07:02",2),
  DRB1=rep("DRB1*01:01",2), DRB3=list(status="absent",copy_number=0),
  DRB4=list(status="absent",copy_number=0), DRB5=list(status="absent",copy_number=0),
  DQA1=rep("DQA1*01:01",2), DQB1=rep("DQB1*05:01",2),
  DPA1=rep("DPA1*02:02",2), DPB1=rep("DPB1*19:01",2),
  MICA=rep("MICA*001",2), MICB=rep("MICB*002",2))
recipient <- donor
recipient$A <- rep("A*01:01",2)
recipient$DRA <- rep("DRA*01:01",2)
run <- neo_indirect_class2(
  donor=donor, recipient=recipient, profile="dr-expanded-v1",
  installation=list(path=netmhciipan_root,
    authorization=list(status="authorized_local_installation",
                       terms=accepted_local_license_description)),
  reference=list(release="3.66.0",cache="ipd-cache"),
  output_dir="indirect-expanded-run")
run$result$group_summary
```

The two installation variables are explicit local setup values. No user creates
or converts prediction TSVs. The 13 configurable donor sources are A/B/C, DRB1,
DRB3/4/5, DQA1/DQB1, DPA1/DPB1 and MICA/MICB. MICA/B are source proteins, not
class II presenters. Character values declare observed allele copies; repeated
values explicitly declare homozygosity, never inferred by the software. A record
may instead state `status="typed"`, `"absent"`, `"untyped"` or `"uncertain"`,
with `alleles`, optional `copy_number` and `evidence`. Missing loci are untyped,
not absent. DRB1 never supplies DRB3/4/5 presence or copy number.

## Reference preparation and exact identity

An explicit release/cache specification retrieves about 17 MB of official
[ANHIG/IMGTHLA](https://github.com/ANHIG/IMGTHLA) protein FASTA, Allele_status,
release metadata and LICENCE, or verifies an existing immutable cache. The
3.66.0 default source pin is `b2055ea33ef464781b9e8ecab20b275b6704f878`
(October 2026). For another explicit release, its official revision is resolved
and frozen; an explicit 40-character `commit` can also be supplied. Every cached
file is SHA-256 checked, and release metadata must agree. `download=FALSE`
requires the cache and fails rather than downloading or repairing changed files.
An incomplete cache is retained; choose a new cache after a failed acquisition.

Reference preparation is local: it writes the existing hla-reference-1.0.0
record layout with full precursor coordinates, source alleles, release, commit,
file hashes and sequence identity. Expanded MICA/B records are consumed only by
this wrapper, not silently enabled in the frozen allo reader. Full status
requires the official Full CDS annotation, Start=1, End=3*(protein length+1),
literal amino acids and full coverage. An exact partial allele is never replaced
by a complete relative. Partial records have no guessed precursor coordinates
and generate no eligible windows. Missing sequences remain unavailable.

A two-field HLA protein name (or the one-field three-digit MICA/B protein name)
may resolve through fully sequenced synonymous members of that exact named
protein group, only if their complete proteins are identical. All supporting
accession/allele names and this protein-reference operation are retained. This
is not genotype completion or nearest-allele substitution. Full allele calls
require their exact records. Ambiguous/uncertain typing must be declared as such;
expression suffixes are not silently regarded as ordinary expressed proteins.
A supplied local `reference=list(path=...,sha256=...)` is also accepted.

The default self universe is the explicitly selected 13 HLA/MIC precursor loci,
not the entire human proteome. Both source and self scopes are written in the
result. Incompleteness is judged against those declared scopes, not against all
possible biology. DRA is validated as a presenter chain; it is not a donor
source in the 13-locus profile. No sequence or missing terminal region is padded.

## DR versus exploratory DQ/DP

`dr-expanded-v1` supports DRA*01:01–DRB1/3/4/5. With one explicitly typed DRA
protein, automatic DR declarations use only explicitly typed beta alleles.
All DR beta loci must be characterized as typed/absent for complete automatic
DR coverage. Automatic DR pairing assumes expression from typing, not from an
expression assay. Its presenter records carry
`expression_basis="typing_based_assumption"`; explicit declarations instead carry
`expression_basis="user_declaration"`. Both record
`experimental_expression_status="not_verified_by_workflow"`: the workflow does
not independently verify experimental expression, including a user's assay claim.
The compatibility field `expression_confirmed` remains the analysis-admission
flag (TRUE for automatic DR), not a measured-expression result. These fields are
saved in the request/plan and public presenter summaries, including JSON/TSV
exports and explanations. They do not change eligibility, scores or counts.
An untyped additional beta locus leaves the overall DR count
unavailable, even if known presenters were evaluated. Alternatively, explicit
`presenters` selects a clearly labeled subset rather than claiming all expressed
DR molecules. Two-field presenter protein names are required; higher-resolution
names are never truncated to make a model match.

`dqdp-exploratory-v1` requires explicit experimentally/otherwise established
pairing declarations; the software does not establish their expression:

```r
pairs <- list(
  list(pair="DQA1*01:01~DQB1*05:01",pairing="confirmed_cis",
       expression_confirmed=TRUE,evidence="explicit fictional research declaration"),
  list(pair="DPA1*02:02~DPB1*19:01",pairing="explicit_expressed",
       expression_confirmed=TRUE,evidence="explicit fictional research declaration"))
# Same call as above, with profile="dqdp-exploratory-v1", presenters=pairs,
# and a NEW output_dir. Both chains must be explicitly typed in recipient.
```

Allowed pairing provenance includes confirmed_cis, confirmed_trans and
explicit_expressed. Unknown pairing is rejected; unconfirmed expression remains
unevaluable. No DQ/DP Cartesian product is constructed. DR and DQ/DP have separate
profile IDs and result schemas. DQ and DP retain separate group totals; no pooled
clinical risk score is returned.

## Predictor and score semantics

Reuse a legitimately acquired [NetMHCIIpan 4.3k installation](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/)
for the host architecture; follow its launcher NMHOME setup. No predictor or
weights are downloaded. The entire physical installation and internal data links
are pinned. Actual -h/-list and every prediction banner must establish 4.3k;
requested exact allele names must exist in the executed model list. Unsupported
pairs remain visible, never fuzzy substitutes. All commands are local CPU batches
with one OMP/BLAS thread, maximum 1024 peptides and a 25 GiB free-space reserve.

Every invocation explicitly requests `-BA`, peptide input, no context, no output
filtering. Vendor-default inversion is retained, or `inversion="all"` explicitly
requests it for all loci. A long output-derived temporary path is not passed.
The previous Darwin buffer-overflow diagnosis and legacy adapter remain intact.

Affinity nM and %Rank_BA are separate from Score_EL/%Rank_EL. Score_BA,
Core_Rel, Inverted, raw strings, exact version, settings and stdout are retained.
EL-only/missing/invalid/truncated/contradictory output fails, even if the vendor
launcher exits zero. BA <1000 nM is the default **explicit assumption** for each
isotype, not cross-model calibration. `thresholds=list(DR=...,DQ=...,DP=...)`
keeps them separate. An exactly equal value is excluded.

An inverted Core maps into the reversed full 15mer at the reported Of; its
forward physical start0 is 15-9-Of. Counts are unique(core,orientation,presenter)
within each group. Self uses the same triple from BA-qualifying self 15mers.
This conservative orientation-aware rule is explicit engineering policy, not a
validated TCR self-recognition rule. Core_Rel is alignment reliability, not an
immune-response probability.

## Results, missingness and reopening

`COMPLETED` means complete only within the declared input scopes and pairing
assumptions. `INCOMPLETE` retains actual evidence and reasons; unknown counts
are NULL/blank, never measured zero. Conditional counts are provisional nonself
relative to observed self, not a complete score. Explicit source absence or a
known complete source with no eligible windows can legitimately yield zero.

- inputs.json, donor_typing.tsv, recipient_typing.tsv, coverage.tsv: all declarations
  and missing/absent/partial states.
- reference.json / reference-provenance.json: locally generated pinned reference.
- origins.tsv: every donor and self window, allele, source occurrence, coordinates
  and sequence identity; no Origin is synthesized from uncertain sequence.
- queries.tsv / query_origin_links.tsv: only exact sequence+pair queries deduplicate;
  all source observations remain.
- raw/ / execution/: original outputs, exact argv/exit and executable/model pins.
- predictions.json/.tsv, model_metrics.tsv: all donor AND self predictions.
- candidates.tsv / candidate_origin_links.tsv: every donor pMHC and full links,
  including unsupported/unscored combinations.
- units.tsv, presenter_summary.tsv, locus_summary.tsv, group_summary.tsv:
  counted contexts and overlapping locus attributions. Do not sum locus rows.
- result.json, explanation.json, run.json, SHA256SUMS.txt: typed results,
  assumptions and complete artifact identity. TSV list cells use lossless JSON;
  zero-row tables may be empty and their authoritative JSON arrays are empty.

Use the exact same call with `resume=TRUE` to reopen a completed/incomplete/dry
run, without inference or window regeneration. Inputs, code, tool, reference and
all artifact hashes must match. Failed runs remain intact and require a fresh
output directory. `dry_run=TRUE` prepares the reference and windows but invokes
no predictor; use a local cache with download=FALSE for a fully offline plan.
Old 3.2/DRB1 usage continues through [neo_indirect_autorun](indirect-autorun.md),
with its unchanged five-locus contract and explicit predictor_version="3.2".

## Rights and reproducible examples

Core code retains PolyForm Noncommercial. IPD's [LICENCE](https://github.com/ANHIG/IMGTHLA/blob/3660/LICENCE.md)
is Attribution-NoDerivs and requires permission before distributing modified
data. Accordingly no processed IPD reference, full protein/derived peptide table,
or model assets are committed. Official source retrieval and local regeneration
are reproducible; factual versions/checksums and independently authored synthetic
software fixtures are included. Cite IPD-IMGT/HLA and the predictor's original
work; upstream terms and possible other rights remain applicable.

The installed `examples/indirect-class2.R` is a runnable MOCK-only example.
The repository `scripts/development/indirect_class2_ipd_validation.R` is an
opt-in genuine IPD/model validation workflow; it requires the separate licensed
installation/cache and writes only to a new local destination. Real model
execution on fictional typing is technical interoperability evidence, not
experimental immunogenicity or clinical validation.


### Independently reproduce the audited technical E2E

After installing this development revision and obtaining the licensed model:

```sh
Rscript --vanilla scripts/development/indirect_class2_ipd_validation.R "/path/to/netMHCIIpan-4.3" "$WORK/reference-cache" "$WORK/real-final"
python3 scripts/development/verify_indirect_class2_ipd.py "$WORK"
```

Set WORK to a new local working directory first. The second command needs Python
3.9+ standard library only. It invokes no predictor: it checks official-source
hashes and full-CDS annotations, every expected precursor window, exact model
alleles, original BA/EL strings, oriented core registers, self subtraction and
all linkage/group/locus counts. It writes audit receipts outside the saved case
folders. The benchmark-specific verifier pins 3.66.0 and the five declared cases;
it is not a generic verifier for arbitrary cohorts/releases. No runtime benchmark
or clinical experiment is performed. Processed IPD/model-output sequences remain
local and are not bundled in the repository. See the validation receipts in
[the development evidence summary](../development/indirect-class2-validation.md).
