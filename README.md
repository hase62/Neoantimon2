# Neoantimon2

Neoantimon2 is an R package for **somatic-variant-derived neoantigen reconstruction
and prioritization with molecular provenance**. It connects declared variants to
transcript/protein paths, mutant peptides, peptide–HLA candidates and independent
prediction/evidence records.

It succeeds the original [Neoantimon](https://github.com/hase62/Neoantimon), described
by Hasegawa et al., *Bioinformatics* (2020),
[doi:10.1093/bioinformatics/btaa616](https://doi.org/10.1093/bioinformatics/btaa616).
That paper describes the predecessor; it does not validate every Neoantimon2
extension. Candidates and model scores are research outputs, not clinical or
experimental proof of immunogenicity.

**This local autorun development branch is version 0.3.0.9000, with 21 public APIs.**
The frozen `cmpb-code-v1` snapshot has 20 APIs; this work does not modify main.
The immutable published [v0.3.0](https://github.com/hase62/Neoantimon2/releases/tag/v0.3.0)
has 16 APIs and does not contain the new numerical-priority, Experimental Yield
or reporting APIs. Record both package version and commit.

![Inputs, reconstruction, peptide–HLA candidates, independent prediction channels and traceable results.](docs/assets/neoantimon2-overview.svg)

[Quick start](#quick-start) · [Inputs](#inputs-and-reconstruction) ·
[External predictors](#external-predictor-workflow) · [Experimental Yield](#experimental-yield-estimation) ·
[Reporting](#complete-result-reporting) · [Validation](#validation-and-scientific-boundaries)

## Key capabilities

| Capability | Implemented scope |
| --- | --- |
| Reconstruction | Somatic SNVs/indels; supported normalized compound/phased contexts, transcript-aware CDS/protein reconstruction and explicit unsupported states |
| Peptides and HLA | Configurable Class I lengths, including 8–12mers; supported Class II workflows with explicit alpha~beta presenters; every eligible window under the documented reconstruction policy |
| Provenance | Event → Origin → peptide–HLA links; sequence/query deduplication preserves all biological Origins |
| Native assessment | Existing A–E Grade and ranking policies; optional configured native predictors |
| Biological evidence | VCF-first DNA VAF, optional DNA/RNA BAM allele evidence, externally quantified expression, generic CN/purity and supported native ASCAT imports |
| External evidence | Separately identified BA/EL/IM observations, exact peptide/HLA matching and model/run/metric provenance |
| Numerical priorities | `neo_priorities()` keeps prediction-only, exploratory adjusted ordering and original-paper PI/PR methods separate |
| Experiment planning | Optional evidence A–D/B1/B2, independent score-bin rates, yield curves and target candidate counts |
| Reporting | Complete Event, unique-peptide, pMHC and Origin views, independent channel representatives and explicit linkage tables |

Specialized structural/source, long-read, longitudinal and indirect-HLA workflows
have explicit prerequisites and limited contracts. This is not arbitrary SV,
splice-event or complex-mutation support. See [advanced workflows](docs/user/workflows.md)
and [scientific limitations](docs/user/methods-and-limitations.md).

## Installation

Use R ≥ 4.1.0. Install the current default-branch development package:

```r
install.packages("remotes")
remotes::install_github("hase62/Neoantimon2", ref = "main", upgrade = "never")
```

For a reproducible published release instead:

```r
remotes::install_github("hase62/Neoantimon2", ref = "v0.3.0", upgrade = "never")
```

Restart R before replacing a loaded package. The commands use normal R package
installation; local source build, isolated installation and examples are checked
in the [integration validation recipe](docs/user/reproducibility.md). Network and
system-library/compiler requirements depend on your platform. For exact
reproduction, pin the audited commit rather than a moving branch.

Human references, predictor engines and model weights are not bundled or installed
automatically. Native predictors and BAM readers have optional dependencies;
see [setup](docs/user/setup.md) and [resource/tool contracts](docs/user/resources-and-tools.md).
Acquire licensed tools separately under their terms.

<a id="try-it-in-five-minutes"></a>
## Quick start

This is an executable **synthetic, offline** example. It downloads no references
and invokes no predictor:

```r
library(Neoantimon2)
source(system.file("examples", "quickstart.R", package = "Neoantimon2"), local = TRUE)
results[, c("sample_id", "ligand_sequence", "mhc_key", "prediction_status", "grading_status")]
neo_explain(p, candidate_id = results$candidate_id[[1L]], format = "text")
reopened <- neo_project(path = p$path)
stopifnot(identical(neo_results(reopened), results))
```

The hand-authored example produces **1 Event, 7 Origins and 4 candidates** using
tiny 3-/5-aa test windows. Prediction is `NOT_RUN`, Grade is `NOT_GRADED`; these
are not negative predictions or grade E. These windows are not human-analysis
settings. The script leaves `p`, `events`, `origins` and `results` in R.

<a id="analyze-your-own-data"></a>
## Inputs and reconstruction

Declare sample/patient identity, single-sample tumor VCFs, patient HLA and verified
assembly-compatible reference resources with `neo_config()` and `neo_project()`;
`p <- neo_run(p)` returns a saved project. `neo_results(p)` reads its results.
One sample column per VCF is supported; do not silently rewrite names or invent GT.
Multiple samples use separate declarations, including samples of the same patient.

For a real analysis, explicitly configure desired lengths, for example:

```r
cfg <- neo_config(peptides = list(class_i_lengths = 8:12,
                                  class_ii_lengths = 15L))
```

This only declares configuration. Supply actual references/inputs before executing
an analysis. Reference bundles are normally several GB; retain their full required
scope and integrity. Class II candidate generation requires supported complete
presenters, not an invented alpha partner.

Supported tumor VCF FORMAT evidence can supply DNA VAF without BAM or a manual
Event-ID table. DNA BAM is optional fallback/verification; RNA BAM can supply
variant-specific allele support. Expression is imported from external estimates;
Neoantimon2 does not quantify expression from BAM. ASCAT is executed externally;
supported outputs are imported with explicit identity/build/coordinate checks.
Missing observations remain unavailable, and reported fractions do not become
invented read counts. See the [continuous user guide](docs/user/getting-started.md)
and [source evidence and numerical priorities](docs/user/source-evidence-and-priorities.md).

<a id="add-external-predictions-if-needed"></a>
## External predictor workflow

1. Read a completed project and export queries with `neo_export_bundle(p)`.
2. Execute a legitimately obtained predictor separately on supported queries.
3. Convert its unmodified output using a reviewed [Converter Registry adapter](https://github.com/hase62/Neoantimon2-Converters).
4. Import with `neo_analyze_external()` and an explicit policy; inspect
   `neo_external_results()` and `neo_external_explain()`.

| Channel | Modeled quantity | Interpretation |
| --- | --- | --- |
| BA | Peptide–MHC binding | Affinity or provider-defined rank; retain original unit/direction |
| EL | Eluted ligand/presentation | Provider-specific presentation score or rank |
| IM | Immunogenicity | Model-specific output, not experimental demonstration of a T-cell response |

Provider, predictor/model versions, run identity, task, metric, unit, direction,
status and all matched Origins are retained. External evidence does not overwrite
Native results. Raw BigMHC-IM scores are not calibrated experimental positive
probabilities. BA, EL and IM scales are not interchangeable and are not combined
into a universal score. The Registry supplies converters, not predictor engines.
See the [external contract and offline example](docs/user/external-predictions.md).

## Experimental Yield Estimation

`neo_yield_evidence()` describes experimental-reference evidence for each exact
peptide–HLA pair, separately from Native A–E Grades:

- **A:** exact peptide and exact HLA tested; retain positive/negative/conflicting records.
- **B1:** exact peptide, different HLA. **B2:** exact HLA, equal-length peptide
  differing at one residue. Both subtypes may apply; neither proves equivalent TCR recognition.
- **C:** sufficient exact-HLA reference records under configurable provisional thresholds.
- **D:** limited reference evidence, not low immunogenicity. An unavailable database
  gives an explicit unavailable status rather than an invented category.

`neo_experimental_yield()` independently assigns four score bins and profile
rates. Evidence A–D/B1/B2 **never modifies the probabilities**. Frozen,
predictor/endpoint-specific empirical profiles and explicitly assumed exploratory
scenarios are distinct modes. The bundled IMPROVE-derived scenario is **not** a
validated BigMHC-IM/BA/EL calibration or patient-specific probability.

Default planning counts each unique peptide sequence once per sample. Its
best-scoring HLA is a labeled primary proxy; all alternative HLA/candidate/Origin
associations remain available. Expected positive count `sum(p_i)` differs from
probability of at least one positive `1 - prod(1 - p_i)`, which assumes independent
outcomes. Targets 50/80/90/95% return unavailable plus a reason when unreachable;
no artificial candidates are added.

Run the explicitly synthetic numerical illustration:

```r
source(system.file("examples", "experimental-yield.R", package = "Neoantimon2"), local = TRUE)
plan <- experimental_yield_example()
plan$summary
plan$targets
plot(plan)
```

This standalone numerical example is not evidence for the genomic demo above.
For saved Native or imported External results, use the
[Experimental Yield guide](docs/user/experimental-yield.md) and an explicit reference/profile.

<a id="understand-your-results"></a>
## Complete result reporting

`neo_report_bundle()` is a read-only export, separate from the existing
`neo_export_bundle()` interoperability/hash contract. For the completed Quick Start:

```r
report <- neo_report_bundle(p, dir = tempfile("neo-report-"))
report$tables$event_summary
report$tables$event_best_by_channel
report$tables$pmhc_candidates
```

Its Event representatives are unavailable in this unpredicted demo, truthfully.
For a completed research project, optionally supply explicit existing external
analyses and a yield plan **derived from that same project's candidates**:

```r
# These variables refer to your existing completed project/analysis/plan:
research_report <- neo_report_bundle(p_saved, dir = "new-research-report",
  analyses = list(existing_external_analysis_id), yield_plan = plan_for_this_project)
```

| View | Meaning |
| --- | --- |
| Event | Original mutation, associations and separately selected Native/exact external channel representatives |
| Peptide | Unique sequence per sample; optional synthesized-peptide planning and primary-HLA proxy |
| pMHC | Every candidate and every matched selected-analysis prediction, including alternative HLAs |
| Origin | Every molecular path, transcript, gene, source Event and available evidence association |

A best candidate never replaces alternatives. Native selection uses existing
rank semantics; external selection stays within a compatible model/run/task/metric,
unit, sample and MHC class. Complete candidate–Origin, Origin–Event and
candidate–peptide linkage tables retain joint multi-Event source sets and separate
germline context; compound peptides are not claimed to arise independently from
each member Event. List fields use JSON rather than being discarded from TSV.

```text
new-research-report/
  manifest.json
  events.tsv, event_summary.tsv, event_best_by_channel.tsv
  peptide_summary.tsv, pmhc_candidates.tsv, pmhc_predictions.tsv
  origins.tsv, candidate_origin_links.tsv, origin_event_links.tsv
  candidate_peptide_links.tsv, molecular_event_evidence.tsv, ...
  experimental_yield_*.tsv    # only when a plan is supplied
```

See [complete reporting](docs/user/result-reporting.md) for selection policies,
full table schemas, precision boundaries, bounded output and current limitations.
Longitudinal aggregate handles require reporting individual underlying projects.

<a id="scientific-limitations"></a>
## Validation and scientific boundaries

[Versioned validation evidence](docs/development/validation-evidence.md) distinguishes
historical evidence from the final integration checks:

- Independent synthetic molecular oracles cover supported sequence reconstruction;
  qualified/unsupported cases remain qualified, not rewritten as unconditional PASS.
- HCC1395 cell-line saved outputs account for 41,072 source variants, 244,191
  Origins, 110,030 pMHC candidates and 949,015 links. The original driver was
  SIGTERM-terminated (controller 143); later saved-result verification is separate.
- Genuine MHCflurry BA and BigMHC EL/IM pilot executions and Converter imports
  were verified, including exact matching of five HCC1395 Class I alleles.
- A genuine bounded pVACseq comparison audited matched contexts and explained
  eligibility/window/WT differences; it does not establish either tool's clinical superiority.
- Yield/reporting tests verify unchanged mathematics, complete alternatives,
  representative/link integrity and read-only reopening. The preserved reporting
  pilot has 15 peptides, 75 pMHCs, 180 links and 300 separate prediction records.

Software correctness, genuine model execution and interoperability are different
from **experimental immunogenicity validation**. None establishes clinical utility.
Modern adjusted BA/EL/IM priorities remain exploratory and candidate-set-dependent.
Published PI/PR mathematical reproduction has its documented missing-input and CCF
limitations. Human reference, annotation, HLA and predictor support must match the
actual analysis. Preserve inputs, versions, hashes and all completed outputs.

<a id="major-capabilities"></a>
## Documentation and references

- [Documentation index](docs/README.md), [practical guide](docs/user/getting-started.md),
  [Public API reference](docs/user/api-reference.md).
- [Inputs](docs/user/inputs.md), [VCF limits](docs/user/vcf-compatibility.md),
  [setup](docs/user/setup.md), [cohort examples](docs/user/multi-sample-examples.md).
- [Grades/explanations](docs/user/results-and-explain.md),
  [source evidence/priorities](docs/user/source-evidence-and-priorities.md),
  [external import](docs/user/external-predictions.md), [advanced workflows](docs/user/workflows.md).
- [Experimental Yield](docs/user/experimental-yield.md), [result reporting](docs/user/result-reporting.md),
  [reproducibility/source layout](docs/user/reproducibility.md), [scope limitations](docs/user/methods-and-limitations.md).
- [AI-assisted analysis playbook](docs/ai/ANALYSIS_PLAYBOOK.md),
  [original Neoantimon paper](https://doi.org/10.1093/bioinformatics/btaa616),
  [software citation](CITATION.cff), [Converter Registry](https://github.com/hase62/Neoantimon2-Converters).

## License

[PolyForm Noncommercial 1.0.0](LICENSE); retain [NOTICE](NOTICE). Commercial use
requires separate permission/licensing from Takanori Hasegawa
(<hasegawa.t.d968@m.isct.ac.jp>); an inquiry does not grant permission.
Third-party tools, resources and converters retain their own terms.

### Optional indirect-HLA autorun (development branch)

`neo_indirect_autorun()` runs the existing literature-derived SOT donor-to-recipient
workflow from explicit HLA typing and a pinned protein reference. It plans donor
and recipient-self 15mers, invokes a user-authorized NetMHCIIpan **3.2 BA**
installation, validates primary 9mer cores, and exports the existing strict
core–presenter count and complete provenance. Dry-run and completed-run
hash-checked resume are available. No vendor software or models are downloaded.
See [setup, one-command use and limitations](docs/user/indirect-autorun.md).
This API is local development work after `cmpb-code-v1`; it is not part of that
frozen tag or immutable `v0.3.0`. It is not official PIRCHE v4/Frost or a clinical
risk predictor. The bundled end-to-end example is explicitly synthetic/mock.


After supplying complete typed HLA/reference inputs and a legitimate local 3.2
installation, one call performs prediction, raw-output conversion and reporting;
users do **not** prepare or edit prediction result files manually:

```r
run <- neo_indirect_autorun(
  allo = list(donor = donor_typing, recipient = recipient_typing),
  presenters = recipient_dr_presenters,
  reference = list(path = "hla-reference.json", sha256 = reference_sha256),
  executable = netmhciipan32_executable,
  model = list(path = model_root, identity = model_release_id),
  authorization = list(status = "authorized_local_installation",
                       terms = accepted_local_license_description),
  output_dir = "indirect-run")
run$result$contexts[[1]]$complete_count
```

The named typing, presenter, executable/model and license variables above are
explicit user-supplied setup, not inferred patient data. See the
[complete input and offline mock example](docs/user/indirect-autorun.md),
[3.2 acquisition and Mac architecture checks](docs/user/netmhciipan32-setup.md)
and [scope and future candidates](docs/user/indirect-autorun-scope.md).
The current autorun scope is human **SOT donor-to-recipient**, explicit
**DRA*01:01–DRB1**, and the literature-derived **NetMHCIIpan 3.2 BA** baseline.
Official PIRCHE v4 uses Frost; numerical agreement and clinical performance
have not been established. A successful mock run is software validation only.
