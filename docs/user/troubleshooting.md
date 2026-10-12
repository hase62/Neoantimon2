# Troubleshooting

[Getting started](getting-started.md) · [Inputs](inputs.md) ·
[Resources/tools](resources-and-tools.md) · [Results](results-and-explain.md)

Start with the reported error, the saved configuration, run manifest/logs and a
candidate's `neo_explain()` output. Execution errors and unavailable biological
evidence need different remedies. Preserve original inputs and irreplaceable
provider outputs before attempting recovery.

| Symptom | What to check / supported remedy |
| --- | --- |
| Wrong/missing file | Confirm spelling, case, permissions and actual file type. Relative input paths use `input_root`; references use `data_root`; use an absolute path to remove ambiguity |
| Config rejected | Use a JSON/YAML analysis template with current loader keys; Use `config/project.example.yaml` or `config/reference.primary.example.yaml`. Unknown/duplicate keys and invalid types fail |
| Resource absent/invalid | Supply the complete selected reference route. Modern primary references need all nine roles; do not mix modern and legacy bundles. A nonempty file can still be invalid |
| Reference/assembly mismatch | Match genome, annotation/catalogs, event build and provider build. Use a new project with a consistent reference; the package does not perform liftover to repair a mismatch |
| Contig missing or REF mismatch | Match literal contigs (`chr1` versus `1`), genome/FAI and variant REF. Normalize/correct upstream against the selected assembly; do not rename just one file or overwrite a run's original input |
| FASTA access/compression error | Genome must be uncompressed with its matching FAI geometry. Catalog gzip extensions must match actual bytes. Reacquire/decompress correctly and retain provenance |
| Partial download/reuse failure | Rerun the acquisition script with the same destination/selection; validated files can be reused. Failed staging is not a completed resource. Do not bypass validators or rename partial files as final |
| Full profile too large | Choose `minimum` or explicit `--resource` selections instead. Plan for multi-GB genomes, population variants, RNA/FASTQ and SRA; verify space on destination/staging volumes before starting |
| Checksum mismatch | Compare the exact original resource/input and declaration. Corruption or scientific input changes require repair to original bytes or a new project with honestly updated provenance; never weaken checking just to reopen |
| HLA malformed/ambiguous | Supply supported literal subject-level declarations. Class II needs explicit alpha~beta pairing. Do not erase expression suffixes, guess an alpha chain, or equate provider-supported alleles with all molecularly resolved alleles |
| No candidates | Inspect events, origins, coding effects, valid peptide lengths, HLA resolution and reasons. Retained noncoding/SV events or unresolved sequences need not yield candidates; do not invent a candidate to ensure a nonempty table |
| Predictor not installed | Either omit it and retain `NOT_RUN`, or obtain the supported executable/models from the provider and configure their actual paths. Installing the R package does not install licensed predictors |
| Provider unavailable/license required | Complete registration/license/manual download with the owner/provider. No credential bypass, vendor redistribution or vendor-file modification is part of this workflow |
| Imported provider output missing | Supply the declared file with the supported provider version, metric, candidate/input/assay scope and provenance. The package does not execute every import-only provider to generate it |
| `NOT_RUN` after a configured attempt | Read primary/provider status records and explain. Somatic priority can project unavailable/unsupported/failed primary evidence to `NOT_RUN`; the run log distinguishes execution errors |
| `NOT_GRADED` | Inspect technical eligibility, primary prediction and source policy. Missing evidence is not grade E; source-extension candidates may intentionally have no somatic grading policy |
| Optional evidence unavailable | Omit the source honestly or supply a supported observed artifact; do not fill missing measurements with zeros or assert negative assays without adequate evidence |
| Multi-sample output empty/unordered | Declare subjects, contexts, relationship/time axes and completed sample runs explicitly. One sample or only matched-normal integration does not create a multi-tumor summary; filenames do not establish time |
| Indirect complete count null | Inspect provider completeness, presenter/sequence coverage and exhaustive qualifying self-search. Partial known positives cannot be reported as a complete count; raw HLA downloads alone do not satisfy the normalized reference contract |
| Reopen/resume mismatch | Reopen with `neo_project(path = ...)` using the saved configuration. Restore original inputs/artifacts if available; changed scientific settings require a new project. A missing run must be regenerated from original inputs and recipe |
| Interrupted lock | Inspect the failed process and run before manual recovery. Locks are not automatically stolen; do not delete an active process's lock or rewrite a checkpoint to mark it successful |

Run directories may be retained. There is no required deletion policy. A completed
run can reuse verified checkpoints; it must not silently accept changed required
inputs, reference bytes or corrupted outputs. Scientific reproducibility comes
from package/source, original resources/inputs and saved configuration/provider
inputs, not an arbitrary historical development run.

For a support report, include package/R version, operating environment, sanitized
configuration, exact error/status/reason and a small synthetic reproducer where
possible. Do not include patient files, private identifiers, credentials, private
links, licensed binaries or personalized vendor download URLs.
