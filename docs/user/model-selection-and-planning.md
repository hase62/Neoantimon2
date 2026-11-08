# Model selection and assumption-aware experimental planning

[User guide](getting-started.md) · [External evidence](external-predictions.md) ·
[Yield](experimental-yield.md) · [Complete reporting](result-reporting.md)

This workflow uses existing public APIs and a sourceable example. It adds no
predictor, probability model or exported API. It neither reconstructs variants
nor writes to a saved project. Use a bounded saved `neo_priorities()` table or a
small existing project's explicitly selected external analyses. Never export or
load a whole genomic project merely to demonstrate a small selection.

```r
library(Neoantimon2)
source(system.file("examples", "publication-planning.R", package = "Neoantimon2"))
# scores: an existing neo_priorities-compatible long table, with explicit hla_class and prediction_scope.
# Select one sample/class; retain candidate_id and every Origin/Event association.
# Each named selector must resolve one exact provider/version/model/run/metric.
# comparison <- publication_model_selection(
#   scores, channels = my_explicit_selectors, N = my_requested_N,
#   profile = my_explicit_planning_profile, pool = "common_evaluated")
```

The example functions `publication_model_selection()` and
`publication_scenarios()` are ordinary sourced example helpers, not package APIs.
Their inputs are explicit; they do not download models, infer missing scores,
create an official interoperability bundle, or automatically choose assumptions.

## Compare selections, not incompatible scores

`neo_experimental_yield()` already selects the best-scoring HLA proxy for each
unique peptide and retains all candidate IDs, Origin IDs and alternative pairs.
The example reuses that ordering. A tie uses peptide sequence, HLA and candidate
identity according to the existing Yield policy. Probability rates do not affect
selection. No BA/EL/IM scores are rescaled or summed across channels.

The supplied profile is needed by the existing API; probabilities in the
selection plans are **not** an assessment of model accuracy and are not used to
rank or compare models. Do not interpret a cross-model planning profile as
predictor-specific calibration.

`common_evaluated` restricts to exact candidate IDs evaluated by **every** selected
channel before choosing the HLA proxy. `model_available` instead ranks each
model's evaluated pool. The declared pool is the union of candidate IDs in the supplied records; it is
not a claim about candidates outside that explicit input scope. Coverage distinguishes evaluated, recorded unavailable
and absent records. Full input records and full-pool plans retain unavailable
statuses/reasons and alternative associations, even when the common pool is
empty. The common-pool proxy can differ from the full-pool proxy.

Outputs include requested/actual N, available pool, coverage, selected sequences,
exact predictor identities, set intersections, model-specific selections and
Jaccard similarity. N is user-specified, including zero; fewer available peptides
are reported without extra invented candidates. An empty union has undefined
Jaccard similarity, not perfect agreement. Samples and MHC classes must be
analyzed separately. This particular composition requires peptide–HLA scope;
peptide-only evidence remains a separate Core-supported workflow.

Prefer independent models of the same task and intended use. Different tasks
require `allow_different_tasks = TRUE`. Comparing BA, EL and IM then describes
**model-and-objective choice**, not an isolated effect of model architecture,
calibration, sensitivity, specificity or clinical accuracy. Two metrics from
one predictor are not two independent models.

## Compare explicit planning scenarios and budgets

```r
# profiles: a named list of user-authored versioned profile objects/JSON paths.
# Do not silently add default low/high rates.
# plans <- publication_scenarios(
#   scores, profiles = my_explicit_profiles, budgets = my_budget_vector,
#   channel = my_explicit_selector)
# plans$budgets; plans$targets; plans$curves
```

Each call passes the **same full pool** to `neo_experimental_yield()`; budget N is
looked up in the resulting cumulative curve. Truncating the input before
cohort-relative binning would change the bin reference population. The helper
requires identical binning definitions, peptide order and HLA proxies across
scenarios. It delegates all expectation/probability mathematics and target-N
calculation to the public API.

Expected positives and the probability of at least one positive are different.
The latter assumes independent outcomes; shared tumor biology and TCR/assay
conditions can violate that assumption. A–D/B1/B2 evidence does not multiply
probabilities. Unknown predictions remain unknown; observed zero and assumed
probability zero remain distinct. Rates 0 or 1 are valid mathematical boundary
assumptions, not automatically plausible biological scenarios.

The four IMPROVE-derived rates (0.0139, 0.0230, 0.0306, 0.0639) are the existing
`improve-exploratory-v0.1` scenario. Applying them to BigMHC-IM rankings is an
explicit cross-model assumption. They are not validated BigMHC probabilities,
confidence intervals or forecasts of an individual patient's immune response.
Additional scenario rates and their provenance must be supplied explicitly.

## Reproducible demonstration and its limits

The retained HCC1395 example has 15 unique sequences, 75 peptide–HLA pairs and
180 candidate–Origin links. Genuine MHCflurry BA and BigMHC EL/IM measurements
can demonstrate selection changes using the same 75 evaluated pairs. They are
three different tasks, so this is not a same-task model-accuracy comparison.
The original scenario yields expected positives 0.5346 and independence-based
P(at least one) 0.4208215262598308. All four target probabilities (50%, 80%, 90%,
95%) are unreachable with those 15 peptides under that scenario.

Keep original raw score strings, model/weight provenance, Converter pins,
analysis IDs, linkage tables and input hashes alongside the derived tables.
The 15-significant-digit public numeric boundary is distinct from raw score
preservation. Figure/table reproduction should use explicit local input paths
and record their hashes. Do not redistribute licensed weights or experimental
reference databases merely because a saved local example used them.

No new full HCC1395 run or predictor execution is required. See the
[versioned evidence inventory](../development/validation-evidence.md) for the
historical full-run/controller distinction and qualified molecular-oracle cases.
An interoperability or reconstruction match is not experimental immunogenicity
validation. Before submission, authors must archive legally redistributable
reproduction inputs/results and verify the target journal's current instructions.
