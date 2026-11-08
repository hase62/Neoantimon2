# Documentation: choose the depth you need

**Choose a version first:** [v0.3.0 (16 APIs), reproducibility-v1 (20), or main (22)](../README.md#choose-a-version).
Main includes the integrated indirect-HLA workflows; the immutable tags do not.

**1. Understand and start:** the [five-minute offline Quick Start](../README.md#quick-start)
needs no predictor or genome download after installation. For indirect-HLA work,
go directly to the [class II one-command guide](user/indirect-class2.md); for the
snapshot's saved tables, use [Reproducibility Snapshot](user/reproducibility.md#snapshot-saved-result-reproduction).

**2. First real-data run:** use the [continuous training and preflight tutorial](user/first-real-data-workflow.md) before switching to real references/VCF/HLA.

For collaborators, the [offline HTML viewer](user/result-reporting.md#offline-html-viewer) exposes saved Candidates, every linked Origin and separate model evidence without R.

**3. Analyze your data:** follow the [practical user guide](user/getting-started.md)
from top to bottom: files/identities → references → configuration → execution →
interpretation → external prediction → reopening. This is the one recommended
own-data entrypoint. [Cohort recipes](user/multi-sample-examples.md) and
[file examples](user/input-output-examples.md) are worked supplements.

**4. Look up a contract:** use the shelf below for details, not as a prerequisite
reading chain.

| Reference question | Page |
| --- | --- |
| Accepted files, sample/HLA declarations and configuration | [Inputs](user/inputs.md) |
| Which VCF features are supported/rejected? | [VCF compatibility](user/vcf-compatibility.md) |
| Resource/tool identities and advanced acquisition | [Resources](user/resources-and-tools.md), [R setup](user/setup.md) |
| Result fields, missingness, grades and ranks | [Results](user/results-and-explain.md) |
| Standard external TSV, policies and execution consent | [External contract](user/external-predictions.md) |
| Indirect HLA: 13 sources, DRB1/3/4/5 and explicit exploratory DQ/DP | [Class II workflow](user/indirect-class2.md) |
| Existing five-source DRB1 runner, 4.3k or explicit 3.2 | [Standalone autorun](user/indirect-autorun.md) |
| Snapshot eight-table saved-result reproduction | [Package and commands](user/reproducibility.md#snapshot-saved-result-reproduction) |
| Longitudinal, vaccine and other extensions | [Advanced workflows](user/workflows.md) |
| Source VAF/BAM, expression/ASCAT and numerical priorities | [Source evidence](user/source-evidence-and-priorities.md) |
| Experimental planning and independent evidence categories | [Experimental Yield](user/experimental-yield.md) |
| Multi-model unique-peptide selection and scenario composition | [Selection/planning](user/model-selection-and-planning.md) |
| Event representatives and complete four-level export | [Reporting](user/result-reporting.md) |
| Measured validation and its limits | [Validation evidence](development/validation-evidence.md) |
| Arguments and return values | [Public API reference](user/api-reference.md) |
| Scientific scope and reproducibility | [Limitations](user/methods-and-limitations.md), [reproducibility](user/reproducibility.md) |
| An error or an unclear result | [Troubleshooting](user/troubleshooting.md) |
| Figure meaning, exports and rebuild instructions | [Figure reference](user/workflow-overview.md) |
| Agent-assisted analysis | [Shared playbook](ai/ANALYSIS_PLAYBOOK.md) |

[Contributing](../CONTRIBUTING.md), [license](../LICENSE), [notice](../NOTICE) and
[citation](../CITATION.cff). The independent Converter Registry handles format
adapter selection/conversion; the user guide explains that boundary before
linking its walkthrough.
