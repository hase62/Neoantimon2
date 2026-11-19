# License scope and sources

Copyright Takanori Hasegawa and Neoantimon2 contributors, for their respective
contributions.

- The independently authored reproduction scripts and this page's explanatory
  text are distributed under the repository's
  [PolyForm Noncommercial License 1.0.0](../../../LICENSE) and
  [copyright notice](../../../NOTICE).
- The archived 20-API Neoantimon2 source in `payload/shared/` retains its
  original [license](payload/shared/LICENSE) and [notice](payload/shared/NOTICE).
  Its archive and source-member hashes are pinned in the replay.
- `payload/` contains the small study inputs, independent expectations and
  replay scripts for the synthetic scopes. `observations/` contains saved
  research observations and projections used for the numerical checks.
  These records are published for research verification with the attributions
  below. Publication does not grant a new blanket license over third-party
  data, output records, or their upstream sources.
- The existing [Reproducibility Snapshot v1](https://github.com/hase62/Neoantimon2/releases/tag/reproducibility-v1)
  remains separate. Its bundled license, notices and file-specific terms
  apply to those release assets.

## Upstream attribution

| Source | Scope in this page |
| --- | --- |
| [SEQC2 / NCBI v1.2.1](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/v1.2.1/) and [NCBI molecular data policy](https://www.ncbi.nlm.nih.gov/home/about/policies/) | Public HCC1395 cell-line records underlie the saved graph and pilot. The page includes derived identifiers, links and selected records, with upstream VCF hashes in [the input ledger](observations/pilot/source-vcf-inputs.json). See [Fang et al. (2021)](https://doi.org/10.1038/s41587-021-00993-6). |
| [GENCODE](https://www.gencodegenes.org/pages/data_access.html) and [Ensembl](https://www.ensembl.org/info/about/legal/disclaimer.html) | Reference-derived annotations and selected sequence views in the saved comparison. Obtain the full reference and annotation resources from their providers. |
| [MHCflurry 2.0.6](https://pypi.org/project/mhcflurry/2.0.6/) | Saved prediction measurements and derived tables. The Apache 2.0 software license applies to the upstream software where independently obtained. See [O'Donnell et al. (2020)](https://doi.org/10.1016/j.cels.2020.06.010). |
| [BigMHC at `199e890`](https://github.com/KarchinLab/bigmhc/blob/199e890fa69a63646730bcecea0a520e08ca4f02/LICENSE) | Saved EL/IM scores and derived tables. The upstream Academic License governs the independently obtained software. The page includes prediction output records, with no upstream code, binary or weights. See [Albert et al. (2023)](https://doi.org/10.1038/s42256-023-00694-6). |
| [pVACtools 7.1.5](https://pypi.org/project/pvactools/7.1.5/) | Saved comparison results and a separate synthetic R09 replay. R09 reuses a separately installed copy of pVACtools and its dependencies under their upstream terms. The runtime inventory consists of file names and hashes. |

The published tree excludes the original private audit archives, machine-local
records, full genomic reference bundles, predictor runtimes and model weights.
