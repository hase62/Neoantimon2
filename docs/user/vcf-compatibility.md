# VCF input scope and compatibility validation

[Inputs](inputs.md#genomic-and-hla-input) · [Multiple samples](multi-sample-examples.md) ·
[Methods and limitations](methods-and-limitations.md)

This guide describes the implemented VCF contract and the development validation
suite. **It is an input-compatibility test, not a peptide reconstruction or
neoantigen prediction benchmark.** No HLA typing, known peptides, licensed
predictors or immunogenicity measurements are required.

## Implemented input contract

- Use tab-separated, site-only VCFs (eight fixed columns) or separately bound
  single-sample VCFs (the fixed columns, `FORMAT`, and one sample). The fixed
  column names/order are required. Any sample column must match the input's
  authoritative `sample_id`; filenames do not assign identity. Multiple sample
  columns are rejected, including tumor/normal pairs. Several separate files can
  share one ordinary project as described in the [cohort guide](multi-sample-examples.md).
- Plain VCF and gzip/BGZF files ending in `.gz` are accepted. Extension and magic
  bytes must agree. BCF is not an input format. A readable file is not a promise
  of full VCF specification conformance.
- The tested sources contain VCF 4.1, 4.2 and 4.3 structures. The parser does not
  enforce a `##fileformat` version gate; unrecognized metadata lines are skipped.
  A missing `#CHROM` header, malformed record width, invalid positive-integer POS,
  duplicated column-header line or inconsistent sample binding fails explicitly.
- Small-variant REF/ALT must be uppercase forward-strand A/C/G/T sequences. REF
  must exactly match the selected FASTA/FAI at the declared contig and position.
  There is no implicit `22`/`chr22` alias repair or REF correction. Recognized
  GRCh37/hg19/b37 declarations conflict with the supported modern reference route;
  declared GRCh38 patch versions and contig lengths are checked where supplied.
- Each selected ALT is normalized separately. Common suffix/prefix trimming
  retains an anchor for unequal-length alleles; repeat indels are left aligned
  against the selected reference. No-op REF=ALT and unsupported sequence syntax
  fail. This phase does not change the existing normalization or event ID rules.
- `somatic_vcf` retains asserted ALTs regardless of GT and FILTER. It does not
  silently apply a PASS-only or genotype filter. `germline_vcf` selects unique
  non-reference GT indices: `0/0` contributes no events, `1/1` contributes one
  occurrence, and `1/2` selects both ALTs. Missing GT fields leave allele selection
  unknown and retain ALTs. An explicitly unresolved/malformed germline GT such as
  `./.` or an out-of-range index is rejected under the existing contract.
- FORMAT fields are looked up by name. Valid omitted trailing sample values are
  treated as missing while original FORMAT text is retained. Empty/duplicate
  keys and excess sample values are rejected. Do not replace missing values with
  numerical zeros.
- Header-declared `Number=A` and `Number=R` INFO/FORMAT arrays are projected by
  the original ALT index; R arrays retain separate REF and selected ALT values.
  Incorrect cardinalities in projected arrays and duplicate INFO keys are rejected. Raw strings,
  literal missing dots and valid zero measurements remain distinct. Whole-array
  `.` missingness retains the released projection behavior; unknown, scalar,
  fixed-length and Number=G fields remain in raw provenance rather than being
  reinterpreted as per-ALT values. Arrays on GT-filtered germline records are not
  projected or fully validated. General Type validation is not promised.
- GT/PS and PID/PGT assertions remain source/sample/contig-scoped. Usable diploid
  phased assertions can establish local phase; missing, non-diploid, invalid or
  conflicting assertions retain explicit states. In particular, phased GT without
  PS is record-scoped in the existing implementation; it does not establish a
  chromosome-wide phase block. No phase-orientation algorithm was changed.
- Recognized symbolic SV and breakend assertions are retained with unsupported
  reconstruction status when their REF, END and mate coordinates are valid.
  Retention does not mean a resolved peptide or general gVCF/SV workflow. Invalid
  mate coordinates, spanning-deletion `*` and unsupported ALT syntax fail.
- Repeated normalized small variants can merge into one event while retaining
  every distinct source occurrence. Fewer unique events than ALT occurrences
  therefore does not necessarily indicate data loss. Import errors are fatal to
  that file operation; no partial event result is returned by the suite.

The ordinary production pipeline still requires its complete compatible
[reference bundle](inputs.md#primary-reference-bundle). The component harness
uses only a verified FASTA/FAI and declared sample context to isolate import,
normalization, merge and per-source phase assertions. It does not establish that
an arbitrary reference configuration can run the complete package workflow.

## Two test levels

**Level 1: routine offline regression.** `tests/testthat/test-vcf-robustness.R`
uses tiny hand-authored fixtures. It checks boundaries, indels/MNV trimming,
multiallelic selection, A/R projections, zero/missing values, genotype and phase,
FORMAT order, malformed declarations, identity, compression and deduplication.
It requires no internet or external VCF tools and runs with the ordinary R tests:

```sh
Rscript -e 'testthat::test_local(filter = "vcf-robustness")'
```

**Level 2: opt-in public-data compatibility.** The source-checkout script
[`validate_vcf_compatibility.py`](../../scripts/development/validate_vcf_compatibility.py)
reuses validation, atomic acquisition and metadata helpers from the existing
[`download_public_data.py`](../../scripts/data/download_public_data.py).
Its [manifest](../../scripts/data/vcf_compatibility_manifest.json) pins original
URLs, checksums, sizes, assemblies, generation methods, provenance and usage
conditions. Payloads remain outside Git. Acquisition/reuse dates are recorded;
unknown original transfer dates for old caches are explicitly null.

Requirements are Python 3.9+, R with package Imports plus `pkgload`, `bcftools`,
`samtools`, `gzip`, and `/usr/bin/time` for measured peak RSS. The measured
versions/platform and every case command are recorded in the run artifacts.
Use a new run directory for each invocation:

```sh
neo_data_root="$HOME/codex/data/Neoantimon2"
neo_run_root="$HOME/codex/runs/Neoantimon2/vcf-compatibility-$(date +%Y%m%dT%H%M%S)"
python3 scripts/development/validate_vcf_compatibility.py \
  --data-root "$neo_data_root" --run-root "$neo_run_root"
```

The default run verifies/reuses the complete pinned cache, prepares small
representative cases and executes the small suite. To acquire missing pinned
inputs, add `--acquire`; large acquisitions additionally require
`--allow-large-downloads`. The complete reference/panel cache is substantial
(about 3.15 GB decompressed genome plus 0.85 GB compressed genome and 0.19 GB
panel), so reuse existing verified data. `--case` selects execution after common
preparation; it does not currently reduce cache requirements. An unknown cache
checksum stops execution rather than silently selecting a different dataset.
No licensed predictor acquisition is performed.

Run the larger case separately with a fresh run directory:

```sh
python3 scripts/development/validate_vcf_compatibility.py \
  --data-root "$neo_data_root" --run-root "${neo_run_root}-large" \
  --large --case seqc2-full-large
```

This imports all 41,072 original SEQC2 records. A per-import timeout defaults to
600 seconds and can be set with `--timeout`; a timeout is a resource-limited
unresolved result, not automatically a parser defect. Peak RSS includes R and
package loading. It is an observation for that environment, not a guaranteed
production resource bound.

## Selected public inputs and explicit preprocessing

| Source | What it exercises | Important boundary |
| --- | --- | --- |
| IGSR/1000 Genomes chr22 integrated phased panel | VCF 4.3, GT selection/phase, INFO allele annotations, sample handling | Original 2,548-sample input is rejected. `bcftools view --no-update -s NA12878 -r 22:20000000-20050000` is separate recorded preprocessing; original bytes are preserved. |
| SEQC2/HCC1395 v1.2.1 through the DeepSomatic public mirror | Site-only VCF 4.1, real SNVs/indels, compression, order, duplicates, larger files | Used for import facts, never peptide or neoantigen truth. Controlled REF/build/corruption/truncation mutations are labeled negative derivatives. |
| Pinned Illumina Strelka2 demo SNV/indel outputs | Caller-specific INFO/FORMAT, no GT, normal/tumor header | Original two-sample files are rejected; TUMOR extraction is explicit. Upstream demo20 is hg19 chr20:862001–867000 remapped to local coordinates, with its own supplied reference; component checks do not introduce GRCh37 production support. |
| Pinned Broad GATK Mutect2 public test resources | Multi-sample and malformed allele-cardinality rejection | The mitochondrial test file declares MMQ Number=A but provides two values for one ALT; it is retained unchanged as an expected rejection, not repaired to create a passing case. |

For the positive IGSR component case, the verified primary chr22 sequence is
explicitly extracted and named `22` to match the source VCF. The derived reference
and its parent checksum are recorded. A separate case checks that the same VCF
fails against the unchanged `chr22` reference. This is a declared reference
choice, not an implicit runtime contig conversion or a complete annotation bundle.

Consult the [IGSR data-use statement](https://internationalgenome.org/IGSR_disclaimer/),
[DeepSomatic source documentation](https://github.com/google/deepsomatic/blob/r1.9/docs/deepsomatic-quick-start.md),
[Strelka demo provenance](https://github.com/Illumina/strelka/blob/0b1a13376e6b49ef9147668bd5e1509b11f7e495/src/demo/README.md)
and [GPL-3.0 terms](https://github.com/Illumina/strelka/blob/0b1a13376e6b49ef9147668bd5e1509b11f7e495/LICENSE.txt),
and [GATK Apache-2.0 terms](https://github.com/broadinstitute/gatk/blob/0cde69eed30339f5978cbb1ac6e5cf3662f9e1f8/LICENSE.TXT).
No upstream payload is committed or redistributed by this repository, and public
access is not treated as a blanket third-party redistribution grant.

## Oracles, outcomes and denominators

An independent Python inspection establishes original record/ALT counts, sample
columns, GT-selected alleles and header-declared A/R attributes. `bcftools norm`
then normalizes independently selected single-ALT queries against the exact
reference. Every emitted source occurrence must match the expected ALT index,
normalized REF/ALT/position, sample/patient, raw attributes and projected values.
Gzip/plain and reversed/duplicated inputs also check stable ordered event IDs.
`bcftools` alone is not used as a full conformance validator: some malformed
cardinalities require independent header/value checks.

Each file-level case is classified as `SUPPORTED_PASS`, `EXPECTED_REJECTION`,
`UNEXPECTED_ERROR`, `INCORRECT_INTERPRETATION` or `UNRESOLVED`. Successful process
exit alone is insufficient. Structured JSON reports include raw records, ALT
alleles, selected occurrences, homozygous-reference filtering, unique emitted
events, diagnostics, runtime and peak RSS. Unreadable corrupt files have unknown
record denominators, not invented zeros. An intended duplicate is accounted for
separately from missing data. Cases often reuse the same underlying input; do not
add their denominators as independent biological evidence.

The suite writes `verified-sources.json`, `prepared-cases.json`,
`environment.json`, `compatibility-results.json`, `compatibility-summary.json`
and per-case oracle/import/process artifacts. Source verification records include
all original sample identifiers, including the complete population-panel header.
Generated derivations also have
an acquisition record under the chosen data root. The manifest/source pins and
transformation descriptions support reproduction without storing payloads in Git.

## Limits and next coverage

This finite collection does not establish universal VCF compatibility. Positive
public multiallelic/Number=R examples are limited; deterministic fixtures cover
those associations more directly. A current, valid single-sample Mutect2 GRCh38
output should be added without altering source headers. GIAB, additional callers,
gVCF blocks, complex SVs/BNDs, more ploidies and multi-gigabyte imports remain
future coverage. General VCF Type/specification conformance and detection of a
plain file truncated exactly at a valid record boundary are not promised. Use
pinned transfer checksums to detect transport-level truncation.
