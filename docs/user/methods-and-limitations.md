# Scientific methods and limitations

[Documentation index](../README.md) · [Workflow details](workflows.md) ·
[Result interpretation](results-and-explain.md)

Neoantimon2 connects canonical input events to reconstructed sequence origins,
then to literal peptide–HLA candidates and scoped evidence. Events, origins, and
candidates are different entities: multiple origins may share one peptide. Exact
sequence and presenter matching preserves every origin association.

Missing evidence is not zero or a biological negative. `NOT_RUN` does not mean a
predictor failed biologically; `NOT_GRADED` does not mean grade E. Provider scope,
provenance, evaluability, and uncertainty remain explicit. There is no universal
immunogenicity or clinical-risk score. Research grades and set selections are not
treatment recommendations. See [grading and explanation](results-and-explain.md).

Native sample-specific analysis (including multiple samples in one project),
subject-scoped longitudinal summaries,
deterministic lexicographic vaccine-set selection, and the strict
literature-derived indirect-HLA profile have separate contracts. The selection
heuristic is not a global optimum. Indirect-HLA is not an official PIRCHE
implementation: its strict profile retains DRA~DRB1 presenters, contiguous 15-mers,
9-aa core self comparison, and strict `<1000 nM` BA evidence. Incomplete scope
has a null complete count. See [supported workflows](workflows.md) for details.

External analysis requires a completed local `neo_project`; an exported bundle
alone cannot execute downstream analysis. Standardized predictions support
separate evidence, user-selected ranking, and eligible somatic grading with the
required semantics and prerequisites. Arbitrary raw scores are not percentiles
and cannot become grades. External imports do not automatically feed all vaccine,
longitudinal, receptor, presentation, or allo modules. See
[external prediction scope](external-predictions.md).

Community OK/NG feedback and converter availability are not official scientific
validation, peer review, or a security guarantee. Downloaded converter scripts
require explicit first-use consent and checksum verification; local scripts are
executed only when explicitly selected. Scripts run with the user's privileges,
not in a sandbox. Licenses for third-party tools and converters remain separate.

The synthetic tests verify implementation contracts and interoperability, not
clinical utility. Real licensed predictors, clinical use, other platforms, and
the minimum supported R version have not all been independently validated.
Acquire licensed tools and reference data separately under their provider terms;
see [resources and tools](resources-and-tools.md). The software citation is not a
peer-reviewed Neoantimon2 validation publication.
