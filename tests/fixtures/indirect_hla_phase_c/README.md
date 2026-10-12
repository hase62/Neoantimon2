# Phase C public reporting fixtures

These files project the frozen, hand-authored Phase B mechanics fixture. They
contain synthetic sequences and imported synthetic binding rows, not real
patient data, historical predictor execution, or official PIRCHE results.

The five result fixtures cover complete strict-profile analysis (5 union units,
1 self-excluded unit), incomplete self search (null complete count), incomplete
recipient DRB1 scope (null complete count, 5 conditional units), an explicitly
declared accessory DRB3 presenter (unsupported, no counted contribution), and a
failed provider row (unevaluable rather than nonbinding).

`compatibility-expectations.json` independently records configuration
distinctions from the frozen reconstruction matrix: S7's inclusive boundary,
S4's unknown equality comparator and changed SOT self domains, S9 Frost's metric
and missing artifacts, S10's separately unpinned build, and the later-v4.x audit
boundary. Every target remains non-official.

Detailed complete/incomplete explanations and a two-ligand/one-core union
example are emitted by the Phase C E2E script into the external review evidence
directory. The tests verify these against independent count expectations,
namespace/reference guards, versioned schemas, and a separate Git-baseline
capture. The runtime still permits one strict SOT scenario per project;
collection tests combine separately computed context projections without
expanding the frozen constructor.
