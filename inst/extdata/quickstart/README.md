# Synthetic quick start reference

SYNTHETIC / TEST-SCALE / NON-CLINICAL. These nine files are byte-identical
copies of the committed `tests/fixtures/reconstruction/` reference at commit
`635823fc4a41308d010396fdababb2ca5d88878b`. They contain hand-authored sequences
and annotation for independent reconstruction tests, not human reference data
or patient observations. The original fixture README documents the oracles.

The installed `examples/quickstart.R` uses these files through `system.file()`.
Its VCF is generated from the existing chr1:4 G>A synthetic oracle. Reference
version/assembly declarations are compatibility labels used by the original
fixture; they do not turn these tiny sequences into GRCh38/GENCODE/MANE data.
The 3-aa/5-aa peptide windows are for demonstration only and are unsuitable for
production prediction. Default real-pipeline execution processes these synthetic
inputs; the package's synthetic/mock execution modes are separate test transports.

Only the example reads this directory. It is not a default resource bundle,
reference lookup, predictor model or scientific selection policy.

| Copied filename | SHA256 |
| --- | --- |
| `genome.fa` | `f1c53bd8b0e05c50e9435aa73fda4c7f1cca07805f734113736da8c68abd4172` |
| `genome.fa.fai` | `4342a7d0dcbe59c7ce3af3f2ec77578e2c5f4be623d58ecd26b4b1be13709a04` |
| `annotation.gtf` | `d5f5ca6ce6d9c1cdd5d1f92bb1a2ed1a4989f643584de8b519e1039e181a69b2` |
| `transcripts.fa` | `164f4bf066cdfdac1829663c77cd8b39536e6607f4ce8ed1584d5039922ce438` |
| `coding.fa` | `164f4bf066cdfdac1829663c77cd8b39536e6607f4ce8ed1584d5039922ce438` |
| `proteins.fa` | `781fb048401564c1b6a21c22e9c4c5c20540171d68547f82e86d5e5c4460d890` |
| `sec.tsv` | `c58768c4acea903be01ed4f291be6edf78dafd4caf1ff64c35d80b41894c5c23` |
| `mane.tsv` | `091e62a223aeb20f3e4075536a10f39b2b4026ad7a2dd033e91ae2efd63a5a01` |
| `mane.gtf` | `71b651954f87e751d8ebc84dc5797ea9c8130753c5b6b27c201a08003e40c0c5` |
