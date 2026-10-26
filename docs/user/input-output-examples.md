# Input and output examples

[Workflow overview](workflow-overview.md) · [Executable cohort examples](multi-sample-examples.md) ·
[Agent playbook](../ai/ANALYSIS_PLAYBOOK.md)

![Illustrative Neoantimon2 inputs and outputs: explicit patients/samples and single-sample VCFs; patient HLA and full Class II pairs; versioned references and optional observations; normalized events, origins, peptide–HLA candidates, distinct prediction metrics and traceable results.](../assets/neoantimon2-input-output-examples.svg)

[SVG](../assets/neoantimon2-input-output-examples.svg) ·
[PNG](../assets/neoantimon2-input-output-examples.png) ·
[Vector PDF](../assets/neoantimon2-input-output-examples.pdf) ·
[Editable source](../assets/source/neoantimon2-input-output-examples.drawio)

The panels are independent illustrative snippets, **not one reconstructed
biological case**. In particular, the illustrative nine-mer SIINFEKLA is not a
claim about translation of chr1:4 G>A in a real human transcript. Transcript TA.1
is a technical-fixture identifier. The prediction numbers are invented format
examples; no upstream model generated them.

The VCF panel abbreviates columns for readability. A complete runnable VCF also
needs its metadata and full header/fields, as provided in the cohort guide.
Each supported VCF has one sample column; separate input records associate
separate files with S1/S2/S3. Patient HLA applies to that patient's samples.
Class II requires both presenter chains; the figure's DRA~DRB1 example is not
interchangeable with a DR beta-only or an unspecified DP/DQ heterodimer.

GRCh38.p14, GENCODE 50 and MANE Select 1.5 name explicit, compatible reference
resources. They are not inferred from a VCF or patient ID. Tiny installed
synthetic references carry contract-test labels and are not human reference
files. RNA, CCF/clonality and MS observations have their own supported evidence
schemas and optional settings; their absence is not an observed negative.

Read results through the supported accessors and explanations. Several origins
may share a peptide and each candidate retains their links. Native NOT_RUN and
NOT_GRADED do not mean zero, biological negative or grade E. BA affinity_nM, EL
presentation_percentile and IM score represent separate metric identities;
IM evidence does not automatically generate a somatic A–E grade. External
rankings and explanations remain separate from native results.

Complete tested R declarations for all four cohorts are maintained in the
[cohort guide](multi-sample-examples.md) and the installed
`examples/multi-sample.R`. The [AI playbook](../ai/ANALYSIS_PLAYBOOK.md) provides
reproducible agent-assisted modes without inventing a new package input format.

Rebuild this supplementary figure with the [same renderer and instructions](workflow-overview.md#maintain-or-reuse-the-figure),
substituting the `neoantimon2-input-output-examples` source/output basename.
