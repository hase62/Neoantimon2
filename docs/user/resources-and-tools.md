# Public resources and external tools

[README](../../README.md) · [Inputs/configuration](inputs.md) ·
[Workflows](workflows.md)

## R setup in v0.3.0

Version 0.3.0 includes R-based reference acquisition, resource checks and
native-predictor configuration. See the [R setup guide](setup.md). These APIs
are absent from the immutable v0.2.0 release. The standalone Python/shell
route below is preserved, and existing analysis configurations remain valid.

## Public acquisition is separate from R installation

Large references and supporting datasets are not included in the package.
Clone the GitHub repository, or download its standalone acquisition script:

```sh
git clone https://github.com/hase62/Neoantimon2.git
cd Neoantimon2
python3 scripts/data/download_public_data.py --help
```

Alternatively, from your chosen script directory:

```sh
curl --fail --location \
  https://raw.githubusercontent.com/hase62/Neoantimon2/main/scripts/data/download_public_data.py \
  --output download_public_data.py
python3 download_public_data.py --dest "$HOME/neo-data" --profile minimum
python3 download_public_data.py --dest "$HOME/neo-data" --resource imgthla,1000g
python3 download_public_data.py --dest "$HOME/neo-data" --profile full --record-sha256
```

Choose a writable destination on a disk with sufficient space. The script's
relative output names are conveniences; absolute resource/input paths may point
to other locations in your analysis configuration. Users have no mandatory
Neoantimon2 physical data-root layout. `NEOANTIMON2_DATA_ROOT` can replace `--dest`.
The default invocation still requires an explicit profile or resource selection.

`--profile minimum` selects `genome,gencode,mane`. `--profile full` selects all
ten resources below. `--resource` takes a comma-separated list using **these
exact identifiers**, without spaces. Profile and resource options may be combined;
duplicates are processed once. There is no `dev` profile. A full acquisition is
optional supporting/development data, not a prerequisite for every analysis.

| Resource ID | Selected identity / content | Additional tools |
| --- | --- | --- |
| `genome` | GRCh38.p14 primary assembly; compressed source, uncompressed genomic FASTA and matching FAI | samtools |
| `gencode` | GENCODE 50 primary GTF, transcripts, protein-coding transcripts/translations, Sec metadata | Base tools |
| `mane` | MANE Select 1.5 GRCh38 summary and Ensembl genomic GTF | Base tools |
| `imgthla` | IMGTHLA 3.65.0 at commit `5b915f27f7f620361cf83cb626eeac8a03c0247c`; original upstream checkout | Git |
| `1000g` | GRCh38 chr22 release `20190312_biallelic_SNV_and_INDEL`, `v2a_27022019`; NA12878 chr22 and 20–25 Mb subsets | bcftools; installed VCF/index support |
| `hcc1395` | Named DeepSomatic HCC1395/SEQC2 chr1 quick-start BAM/reference/truth artifacts | samtools, bcftools |
| `processed-rna` | Named HCC1395 tumor/normal HISAT tutorial BAMs; no immutable published release claimed | samtools |
| `longread` | Original SRA accession `SRR23456709` | SRA Toolkit `prefetch`, `vdb-validate` |
| `full-rna` | Original SRA accession `SRR37849525` | SRA Toolkit `prefetch`, `vdb-validate` |
| `full-truth` | SEQC2 truth VCF 1.2.1 and high-confidence BED 1.2, GRCh38 | bcftools |

Base tools are Python 3, curl and gzip; selection-specific validators also use
the tools shown. The optional repository RNA/alignment workflows may need
`fasterq-dump`, `sam-dump`, STAR, samtools, bcftools, bgzip/tabix as appropriate.
They are separate shell workflows, not a requirement to load the R package.
Original SRA acquisition does not automatically generate FASTQ/aligned BAM.

The primary genome alone is roughly 3 GB uncompressed, with its compressed
download and staging adding more space. Catalogs add further space. Population
VCF, BAM and SRA selections can be multi-GB; full-RNA conversion/alignment may
need substantially more space for FASTQ, indexes and temporary files. No fixed
total size is promised for `full`: inspect selected upstream sizes and allow
space for staging plus completed data before downloading. Prefer individual
selections when you do not need all supporting datasets.

Download the matching primary reference configuration separately when using the
standalone script:

```sh
curl --fail --location \
  https://raw.githubusercontent.com/hase62/Neoantimon2/main/config/reference.primary.example.yaml \
  --output reference.primary.yaml
```

Set your `data_root` to the acquisition destination and your `run_root` to a
separate generated-output location. See [configuration](inputs.md#configuration-precedence-and-templates).

## Completion, reuse and provenance

New downloads use owned staging paths on the destination filesystem. The
downloader must exit successfully and the file must pass the appropriate
format/readability checks before promotion. A nonempty partial file is not a
completed download. Checks include gzip integrity, genome/FAI geometry,
VCF/index queries, BAM quickcheck/index queries, and SRA validation as applicable.
Each promoted file is replaced atomically; do not run concurrent acquisitions
into the same destination because index/data promotion is not a multi-file
database transaction.

Existing validated files can be reused. An invalid existing artifact is not
accepted as successful merely because it exists. An existing IMGTHLA checkout
with a conflicting commit/release or local modification is rejected without
changing it: select a separate destination rather than overwriting your changes.

`acquisition/<resource>.json` records selected URLs, release/accession, semantic
scope, artifact locations/sizes and observed tool versions. `--record-sha256`
adds observed file SHA256 values where file artifacts are recorded. Universal
predetermined SHA256 values are **not** required for every public resource.
An observed checksum establishes the bytes used; it does not by itself prove
assembly, catalog scope or scientific compatibility.

Do not let a mutable `latest` release silently replace a scientific reference.
Use the selected release/scope, retain provenance, and create a reviewed new
configuration for changed references. Named unversioned tutorial/quick-start
URLs are recorded honestly and do not guarantee immutable future bytes.
Raw IMGTHLA is not an already normalized literal HLA reference manifest for allo
analysis; typing, sequence content, coverage and release policy must be supplied
under that separate contract.

Upstream data remain subject to their terms. The script records applicable
source/terms references, including [GENCODE](https://www.gencodegenes.org/pages/data_access.html)
and [IPD-IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/licence/); acquisition does
not grant redistribution rights.

## External tools and provider inputs

The table describes **current Neoantimon2 adapters**, not everything a provider
can do. No tool/model is installed by `neo_run()`. Real external execution is
optional; obtain and test your local installation separately. Local-provider
execution has additional environment requirements beyond the ordinary R Imports.
Provider values preserve their own units/tasks/scopes; they are not interchangeable.

| Tool/provider | Purpose and current integration | Version/scope and user-supplied paths | Acquisition boundary |
| --- | --- | --- | --- |
| NetMHCpan | Local Class I binding/EL execution via `predictors` | Standard 4.1/4.2 EL/BA output families, including detected patch suffixes; explicit `executable` or PATH. Adapter lengths 8–14, supported A/B/C syntax. Optional newer transfer-learning modes are not requested | Manual official DTU acquisition/registration/license; no vendor redistribution |
| NetMHCIIpan, somatic adapter | Local Class II binding/EL execution | 4.3 family; explicit alpha~beta presenter; `executable` or PATH, adapter lengths 9–25 | Manual official DTU acquisition/license |
| NetMHCIIpan, indirect-HLA import | Import externally executed BA evidence; **no local indirect predictor runner** | **3.2-family BA**, pinned build/model/query/result identity, `indirect_hla$provider$path`. 4.x EL is not a substitute | Authorized manual provider acquisition/execution; supply provenance, not redistributed binaries/models |
| MHCflurry | Local Class I binding/presentation CSV adapter | 2.x; `executable` (`mhcflurry-predict`) and explicit nonempty `model_path`; lengths 8–15. Binding-only option restricted to binding task | User installs package/models via provider instructions and terms; not downloaded by the public-resource script |
| ANNOVAR / VEP | Event-table imports only | `somatic_annovar`/`germline_annovar`, or supported VEP roles; user supplies output path/version/build. No annotation runner | ANNOVAR registration/license/personal links remain manual. Follow VEP's official installation/data terms externally |
| FACETS | Copy-number `cncf` import | `somatic_evidence` kind `facets_cn`, output `path`, source version and `reference_build`; no FACETS execution | User-generated output; provider software acquired separately under its terms |
| LOHHLA | Named-column output import; optional explicit execution wrapper | `hla_loh$sources` or `hla_loh$execution`. Wrapper needs script, BAMs, HLA calls/FASTA, purity/ploidy and configured novoalign/GATK paths plus LOHHLA's R dependencies | User supplies official software and any restricted dependencies; do not bypass their acquisition requirements |
| STAR-Fusion / Arriba / SpliceAI / STAR junctions | Source/support imports; native reconstruction from supplied structures | Exact supported native named fields, versioned input roles/output paths; no caller execution through R | User runs/acquires providers separately under their terms |
| SQANTI3 / FLAIR / IsoQuant | Long-read exon-structure/metadata imports | Supported GTF/GFF representation plus explicit version; optional supplied transcript FASTA/settings | User runs providers externally; raw long reads are not inferred isoforms |
| MHCflurry processing / NetChop / NetCTLpan / NetMHCstabpan | Processing/stability **output imports** via presentation settings | MHCflurry 2/class1_processing; NetChop 3.1 C-term-3.0 or 20S-3.0; NetCTLpan 1.1; NetMHCstabpan 1.0 family. `processing_output` and `output_input_id`; configured executable/model detection does not execute them | DTU tools require official manual acquisition/license as applicable; user manages MHCflurry models. No guessed processing CLI |
| PRIME / pMTnet / ERGO-II / NetTCR | Recognition **output imports**, independent of pMHC/priority | PRIME 2.1 (MixMHCpred >=3.0 dependency); configured pMTnet/ERGO-II model/features; NetTCR-2.0, Class I 9-mer/A*02:01 initial scope. `recognition_output` + settings; no real recognition runner | User-managed software/models/backgrounds and applicable provider terms; manual registration/license steps are not automated |
| AIRR / VDJdb / MS / arcasHLA | Repertoire, curated annotation, measured ligand or expression imports | Supported named tables/AIRR JSON/mzIdentML and explicit assay/sample/version scope; user-supplied output paths | Consult selected dataset/provider rights; no credentials/private or restricted datasets bundled |
| PyClone-VI / generic clonal provider | Imported prevalence/cluster evidence | Supported exported table, subject/run/version, exact event mapping; no de novo clustering | User runs provider externally and supplies output/mapping under its terms |

Other historical inventory names do not establish current adapter support.
For example, listing an ASCAT/Sequenza/CNVkit/PureCN directory does not make raw
outputs accepted by the current FACETS parser. Do not relabel another format as
FACETS. Fixture providers are for explicit synthetic tests, not production models.

Official starting points: [NetMHCpan 4.2](https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/),
[NetMHCIIpan 4.3](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/),
[MHCflurry](https://openvax.github.io/mhcflurry/commandline_tools.html),
[ANNOVAR acquisition](https://www.openbioinformatics.org/annovar/annovar_download_form.php).
These links are not substitute licenses or personalized download links. If a
provider needs an account, email registration, license acceptance or credentials,
complete that process manually. Keep vendor binaries/models intact; Neoantimon2
compatibility belongs in its own wrappers. Never put credentials in committed
configuration, shared reports or acquisition metadata.
