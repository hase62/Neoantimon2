# Inputs and configuration

For the v0.3.0 release candidate, [R-based setup](setup.md) can supply
reference configuration fragments directly. The declarations below remain the
analysis contract, including when resources were obtained by Python or shell.

[Getting started](getting-started.md) · [Resources/tools](resources-and-tools.md) ·
[Workflows](workflows.md)

## Samples, subjects and paths

`samples` is a data frame or array of records. Every record needs unique
`sample_id` and explicit `patient_id` (the subject identity in somatic projects).
Optional `timepoint`, `tissue` and `role` are strings. Do not infer subject identity
from filenames. Subject-level HLA is a named list keyed by `patient_id`; declare
one genotype per subject across its samples.

See the [four executable multi-patient/multi-sample examples](multi-sample-examples.md)
for complete declarations, offline setup and sample-specific result selection.
Multiple separately bound single-sample VCFs can be processed in one project;
this does not enable multiple sample columns inside a VCF or automatically
create longitudinal summaries.

Each `inputs` record needs `input_id`, `input_role`, `sample_id` and `path`.
`sample_id` must belong to the project. `input_id` is a unique instance label;
`input_role` chooses the actual parser/meaning. Optional `version` records the
provider/parser release; fusion/splice and long-read inputs require it. Optional
file checksum is `list(algorithm = "sha256", scope = "file", value = ...)`, with
64 lowercase hexadecimal characters. Source-specific configuration records may
instead accept a plain SHA256 string; follow the relevant format.

Absolute and home-relative paths can point anywhere. Relative inputs resolve
against `input_root`; relative references against `data_root`; relative logical
`tools.*`/`models.*` resources against `tools_root`. Input/resource roots are
conveniences, not an obligatory directory layout. Predictor request
`executable`/`model_path` fields are used by their adapter: prefer absolute paths
there rather than assuming logical-resource root resolution. Allo manifest and
vaccine-policy filenames supplied in a config file are anchored at that file's
directory. Use absolute manifest paths in direct R overrides.

Keep original inputs/resources outside run roots; overlapping trees are rejected.
Missing declared inputs are errors. An omitted optional evidence source stays
unavailable/not evaluated; it is not created with zeros.

## Primary reference bundle

The current public reference example declares GRCh38.p14 primary assembly,
GENCODE 50 and MANE Select 1.5. Nine logical roles are required for this modern
reference route. Catalogs are retained; MANE supplements selection rather than
replacing the primary GTF's coordinate authority.

| Logical resource | Format / role |
| --- | --- |
| `reference.genome_fasta` | Uncompressed genomic FASTA for indexed sequence access |
| `reference.genome_fai` | Matching five-column FAI: contig, length, offset, line bases, line width |
| `reference.annotation_gtf` | Primary-assembly GENCODE GTF; `resource_id = "gencode_primary_assembly_all_catalogs_mane_supplemental_v1"` |
| `reference.transcript_fasta` | Transcript FASTA catalog |
| `reference.protein_coding_transcript_fasta` | Protein-coding transcript FASTA catalog |
| `reference.protein_fasta` | Protein translation catalog and exact-self lookup scope |
| `reference.selenocysteine_metadata` | GENCODE Sec metadata; literal U/translation exceptions stay explicit |
| `reference.mane_summary` | MANE summary with transcript mapping/status and GRCh38 fields |
| `reference.mane_gtf` | MANE Ensembl genomic GTF supplement |

GENCODE/MANE catalog files may be plain or gzip according to their actual filename
and bytes; a compression/extension mismatch fails. Do not replace a missing
protein catalog with an arbitrary proteome. Selected assembly, release, transcript
versions, coverage and scope must agree. Genome FASTA/FAI geometry and input REF
alleles are checked. Contig aliases such as `1` and `chr1` are not silently repaired
by the runtime. Selected unsupported or incomplete sequence contexts retain
uncertainty; they are not completed by homology or an arbitrary isoform.

The legacy reference route uses `reference.refFlat` and `reference.refMrna`.
It must not be mixed with modern reference roles. A predecessor reference version
or old directory is not a default for Neoantimon2. Prefer the documented current
primary profile for a new analysis.

## Genomic and HLA input

| Input | Requirement / format | Computation and missingness |
| --- | --- | --- |
| Somatic variants | For the ordinary somatic example: `somatic_vcf`; VCF plain/gzip, explicit sample binding | Native REF validation, small-variant normalization, phase/reconstruction; no event input yields a failed ordinary production run |
| Germline variants | Optional `germline_vcf`; same reference/sample semantics, genotype/phase retained | Germline context can affect reconstruction and comparison; absence does not mean a fully observed germline genome |
| MAF | `somatic_maf` / `germline_maf`; named coordinates/alleles, optional declared build/sample | Imported alleles are reference-validated; imported consequence labels do not replace selected reference reconstruction |
| ANNOVAR output | `somatic_annovar` / `germline_annovar`; TSV `Chr`, `Start`, `End`, `Ref`, `Alt` | Import only; no ANNOVAR installation/execution. Indels require REF-validatable anchor/interval semantics; `Ref=0` is rejected |
| VEP output | `somatic_vep` / `germline_vep`; supported named-coordinate/allele table | Import only; not general automatic VEP execution |
| Other event representations | `canonical_events`, `germline_snp_table`, `somatic_bnd_table`, `sequence_fragment`, `sequence_entire_region` | Dedicated parsers, not arbitrary free-form inputs. BND/SV retention does not guarantee peptide reconstruction; sequence provenance/coverage limits remain explicit |
| HLA declaration | Needed to create peptide–presenter candidates; `hla = list("subject-1" = c(...))` | Literal declared alleles normalized to supported presenter identities; no typing inference. Omitted/unresolved HLA can yield no candidates despite retained events/origins |

The supported VCF has standard `#CHROM POS ID REF ALT QUAL FILTER INFO` fields;
when a sample column exists it must match the bound `sample_id`. Multi-sample
VCFs are rejected: supply a separately selected single-sample input. Number=A/R
fields preserve their selected ALT relation. See the [VCF scope and compatibility suite](vcf-compatibility.md)
for genotype selection, missingness, parser limits and reproducible public-input tests. GT/PS and PID/PGT phase evidence are
source/sample/contig-scoped; missing/conflicting phase is not resolved by guessing.
Declarations of reference/contig build and contig length must match the selected
reference. Small-variant REF/ALT must be supported genomic alleles; symbolic SVs
are not treated as resolved substitutions.

Class I example: `HLA-A*02:01`. Class II requires an **explicit alpha~beta pair**,
for example `HLA-DRA*01:01~HLA-DRB1*04:01` or a supplied DQ/DP pair. A beta chain
alone does not infer an alpha chain or cis/trans pairing. Expression suffixes and
typing ambiguity must not be silently erased. A molecularly resolved presenter
can still be unsupported by a particular predictor's narrower scope.

## Optional evidence and source extensions

All optional paths are supplied by the user. None is a mandatory public bootstrap
download. Tables preserve observed/raw values and versioned provenance.

| Category | Current entry / format | Native versus imported; absent input |
| --- | --- | --- |
| DNA/RNA read evidence | `neo_config(somatic_evidence = list(sources = ...))`, kinds `dna_reads` / `rna_reads`; CSV or TSV keyed by exact `event_id` | Imports counts; derives VAF/technical criteria. No raw BAM pileup in the R input API; missing counts remain unavailable/not evaluated |
| Raw DNA/RNA BAM, FASTQ | Upstream provider inputs; optional repository shell workflows use explicit paths | Not accepted as ordinary `dna_bam`/`rna_bam` event input roles. Generate supported evidence externally; no automatic alignment/caller execution through `neo_run()` |
| Expression | Somatic source kind `expression`; transcript/gene IDs, `transcript_tpm` / `gene_tpm` | Import; matching transcript links to origins. Missing source/match is distinct from measured zero TPM |
| Purity/CN/CCF/generic LOH | Source kinds `purity`, `facets_cn`, `ccf`, `loh` | FACETS result import; conditional CCF derivation only when prerequisites exist. No automatic CN calling or clonal inference |
| Fusion | Inputs `star_fusion` / `arriba`; native provider named-column output and explicit `version` | Import breakpoints/support; reconstruct compatible transcript structures/ORFs. No STAR-Fusion/Arriba execution; candidates stay source-specific `NOT_GRADED` |
| Splice | `splice_structure` plus optional `spliceai_vcf`, `star_junction` evidence | Supplied structure enables reconstruction; SpliceAI/STAR support alone is not a fabricated splice structure |
| Long-read | `long_read_gtf`, `sqanti3_gtf`, `flair_gtf`, `isoquant_gtf`; matching metadata roles, optional `long_read_fasta` / `long_read_settings` | Import explicit exon structure and version, reconstruct compatible isoforms/ORFs; no aligner/isoform-caller execution. Incomplete ends/evidence stay explicit; source-specific `NOT_GRADED` |
| HLA LOH | `neo_config(hla_loh = ...)`; named-column LOHHLA output | Import by default; optional explicit execution wrapper has additional dependencies. Allele loss is separate from generic CN LOH and does not rewrite somatic grade |
| Processing/stability | `presentation_settings` + `processing_output` | Imports supported MHCflurry/NetChop/NetCTLpan/NetMHCstabpan output. Independent origin/flank or peptide–presenter scope; no automatic real-provider runner |
| MS ligands | `immunopeptidomics_psm` TSV or `immunopeptidomics_mzidentml` XML, with experiment scope in presentation settings | Imports PSMs/assay observations; no search engine execution. Missing MS is not evidence of absent presentation |
| HLA/APM expression | `hla_expression`, `arcas_hla_quant`, `apm_evidence` with named measurements | Import only; genotype, expression, allele loss and APM remain independent |
| TCR/recognition | AIRR/receptor/reactivity/VDJdb roles and `recognition_settings` / `recognition_output` | Imports repertoire/assay/provider evidence; no real recognition predictor execution. Pairing, functionality and assay completeness stay explicit |
| Longitudinal/clonal prevalence | `longitudinal_settings`, `clonal_evidence`, `pyclone_vi`, `clonal_mapping` | Subject-scoped overlay/import; exact mapping required. Provider prevalence is separate from event CCF; no averaging or inferred global clone |
| Allo/indirect-HLA | `allo`, `indirect_hla`, optional `allo_native_evidence` config | Explicit context, typing, normalized literal HLA reference and pinned provider JSON. Raw IMGTHLA alone is insufficient; see workflows |

For somatic evidence each source has `kind`, `sample_id`, `path`, and optional
`source`, `source_version`, `reference_build`, `checksum`. Do not declare two
sources of the same kind/sample. Event tables join exact event IDs, not gene names
or row order. You can first materialize/import events to inspect their IDs, then
create a new project with a saved, explicitly mapped evidence table.

DNA columns include `event_id`, `tumor_ref_reads`, `tumor_alt_reads`, `tumor_depth`,
`normal_ref_reads`, `normal_alt_reads`, `normal_depth`; corresponding RNA fields
use `rna_`. Missing values are null/NA, not invented zeros. Purity is a fraction;
CCF is not VAF. FACETS requires declared `reference_build` and named `chrom`,
`loc.start`, `loc.end`, `tcn.em`, `lcn.em`; optional `normal_cn`, `purity`, `ploidy`
support scoped derivation. Only a unique segment fully containing the event can
be assigned. No assumed diploidy, purity of one or mutated-copy multiplicity.

Generic PSM TSV requires `psm_id`, `spectrum_id`, `peptide_sequence`; retain
quality/decoy/modification fields and explicit experiment/IP/HLA/sample scope.
A pan-HLA observation does not assign a unique presenter. AIRR Rearrangement TSV
requires `sequence_id`, `locus`, `junction_aa`; paired chains require actual
pairing evidence. AIRR JSON and other receptor tables have distinct parsers.
The full accepted role names and versioned settings shapes are in the
[workflow guide](workflows.md); installed schemas can be located with
`system.file("schema", package = "Neoantimon2")`. They are technical format
references, not permission to populate unknown fields with convenient defaults.

## Configuration precedence and templates

`neo_config()` applies defaults → JSON/YAML file → environment → explicit named
arguments. Unknown keys/types and duplicate keys are rejected. `env = character()`
disables environment overrides for a saved recipe. Top-level settings include
`paths`, `scientific`, `runtime`, `checksum`, `resources`, `peptides`, `predictors`,
`somatic_evidence`, `somatic_priority`, `hla_loh`, `vaccine_set`, `allo`,
`indirect_hla`, `allo_native_evidence`.

| Setting | Meaning |
| --- | --- |
| `paths` | `data_root`, `input_root`, `tools_root`, `run_root`, optional `repository`. File-relative roots are anchored at the config file; direct argument/environment roots must resolve absolutely |
| `scientific` | Species taxon (default `NCBITaxon:9606`), selected `genome_assembly`, optional `policy_id`; changing these changes scientific configuration |
| `runtime` | Positive integral `threads`; `verbosity = "quiet"`, `"info"`, or `"debug"` |
| `resources` | Named logical records with `path`, `type`, `required`, `resource_id`, `version`, `checksum`, `source`. Types: file/directory/executable |
| `checksum` | `metadata` default observes path/stat identity; `sha256` observes/verifies content; `none` is not content verification. Declaring a checksum is not the same as universally enabling verification |
| Optional policies | Supply only when their evidence/workflow is intended; defaults do not install a predictor or enable every overlay |

Directory SHA256 hashing through logical resource integrity is rejected; predictor
model directories have their own adapter-specific file identity. Reference-bound
declared checksums are additionally verified by reference validation. Missing
`required = TRUE` resources fail; optional missing resources do not establish
biological absence. Standard roots default to R's per-user data directory, not
the source checkout. Default `input_root` is `data_root`; `tools_root` is its
`tools` child. These defaults can be overridden.

Environment suffixes are `DATA_ROOT`, `RUN_ROOT`, `INPUT_ROOT`, `TOOLS_ROOT`,
`REPOSITORY`, `THREADS`, `VERBOSITY`, `CHECKSUM`, under `NEOANTIMON2_`.
`neo_project(reference = ...)` overrides matching config resources using role
names without the `reference.` prefix.

The current [primary template](../../config/reference.primary.example.yaml)
contains all nine required reference roles. The
[project loader example](../../config/project.example.yaml) demonstrates loader
syntax but does not by itself supply a complete biological reference bundle.
Use the current project and primary-reference templates for analysis; historical
resource inventories are not loader configurations.

## Optional predictor and quantitative-evidence configuration

Add this fragment to a **complete** primary-reference YAML configuration when
you intend to use an installed predictor and your own measured evidence. Replace
the executable, evidence paths, sample identity and source versions with your
actual values. This is a configuration template, not bundled provider software
or synthetic measurements:

```yaml
predictors:
  - provider: netmhcpan
    task: mhc1_presentation
    executable: /opt/providers/netMHCpan
    options:
      timeout_seconds: 600
somatic_evidence:
  sources:
    - kind: dna_reads
      sample_id: sample-1
      path: /analysis/sample01/evidence/dna-counts.tsv
      source: upstream-read-counts
      source_version: your-pipeline-version
    - kind: rna_reads
      sample_id: sample-1
      path: /analysis/sample01/evidence/rna-counts.tsv
      source: upstream-read-counts
      source_version: your-pipeline-version
    - kind: expression
      sample_id: sample-1
      path: /analysis/sample01/evidence/expression.tsv
      source: upstream-expression
      source_version: your-pipeline-version
```

`predictors` is an array of requests, not a named provider map. Request fields
are `provider`, `task`, optional `executable`, `model_path`, `options`. Current
provider identifiers are `netmhcpan`, `netmhciipan`, `mhcflurry`; their respective
Class I/II tasks are `mhc1_binding`, `mhc1_presentation`, `mhc2_binding`,
`mhc2_presentation` within each adapter's supported class. Version is detected
from the installed provider, not forced by an unsupported `version` request key.
MHCflurry needs its real model directory in `model_path`. Options include
`timeout_seconds`, explicit `retry_token`, and MHCflurry-only `affinity_only`
(TRUE requires its binding task). A changed retry token creates a new provider
request identity; it does not make failed output successful.

Use the [getting-started recipe](getting-started.md#configure-your-own-project)
to load the completed configuration and run. Set ordinary production peptide
lengths rather than the demo's 3-aa/5-aa windows. Omit unused optional source
records instead of pointing at nonexistent placeholder files. Keep raw provider
measurements and source versions; do not invent counts to obtain a grade.
