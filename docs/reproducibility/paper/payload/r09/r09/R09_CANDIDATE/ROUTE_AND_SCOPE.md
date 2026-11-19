# R09 route and scope — frozen before either native report/import execution

Date: 2026-11-12. This is a post-prediction information/reporting comparison, not a molecular reconstruction, model accuracy, clinical, or unknown-phase comparison. No superiority is presupposed.

## Shared canonical input and expected relations

One synthetic sample `NV_SAMPLE`. Peptide `MKACDFGH`, window1–8, mutation E2K against WT `MEACDFGH`. Two exact Class I HLAs: HLA-A*02:01 and HLA-B*07:02. Three resolved molecular source contexts: chr1:4 G>A / AUD1.1; the SAME genomic Event / AUD2.1; independent chrB:4 G>A / BUD1.1. Each has peptide window1–8. Thus two source Events, three contexts, one sequence, two peptide–HLA identities and six context–HLA associations. No uncertain-phase context enters the primary comparison.

CONTROL_NOT_INFERENCE MT affinity inputs: 100 nM for HLA-A*02:01; 200 nM for HLA-B*07:02, repeated across source contexts because the peptide/HLA identities are identical. These are stipulated test numbers, not measured affinities or trained-model output. WT affinity, all percentiles, RNA expression and RNA measurements are MISSING/NA, never numeric zero. Source DNA support is synthetic tumor8/20, normal0/20. No expression is invented. Core_native prediction remains NOT_RUN/NOT_GRADED.

Expected Q1: both HLA queries recover exactly the three source contexts and two Events. Q2: global unique sequence count1; two distinct HLA scores100/200 with nM/lower-is-better identity. Q3: excluding chr1 Event leaves BUD1 only; excluding chrB leaves AUD1/AUD2. Both retain the sequence and both HLA associations. Q4: excluding both leaves no support. Q5: same results after new-process rereading. All representative choices/exclusion policies are predeclared, not selected for favorable results.

## Boundaries and actual operations

- Neoantimon2 uses unchanged repaired Core archive SHA `5adabad73f3b484338cae8cc594d574f6022eb9040f717f3c276686a1f2a6c1b`. If a matching preserved resolved project cannot be directly reopened at its saved paths, one new tiny resolved project will be constructed from the supplied synthetic reference/VCF; this upstream preparation is NOT compared against pVACseq reconstruction. New CONTROL_NOT_INFERENCE scores are imported through `neo_analyze_external`; results/explanations are obtained through public APIs. Existing native records remain unchanged.
- pVACtools7.1.5 is reused without edits. A byte-identical isolated copy of the installed Python package and entrypoint is used with its already installed dependencies; no downloads/inference. The initial pVAC input is a schema-conforming post-prediction `R09.CONTROL_INPUT.all_epitopes.tsv` constructed from the SAME canonical rows. It is explicitly INPUT, **not a pVACseq-produced raw report**. Genuine native outputs are produced by the public `binding_filter` CLI and `generate_aggregated_report` CLI. This task does not execute `pvacseq run` or its upstream prediction/output parser.
- `MHCflurry ... IC50 Score` column names are schema carriers required by existing code; they do NOT mean MHCflurry was run. An explicit CONTROL_NOT_INFERENCE sidecar and input columns identify the synthetic origin. Missing WT/percentiles/RNA stay NA. Required transcript length/biotype/CDS completeness come from synthetic fixtures; required TSL/MANE/Canonical labels are controlled schema annotations, not genuine VEP output.
- Full native binding-filtered output is examined, including custom canonical context/sample columns passed through by the filter. Aggregate TSV AND metrics JSON are examined together. The input all_epitopes table is inspected separately and never relabelled as observed output. No pVACtools observed output is hand-written.
- Runs: D0 uses unmodified aggregation defaults (binding threshold500nM, aggregate inclusion5000nM, count limit15, top-score metric median, secondary metrics ic50+combined_percentile). D1 explicitly uses `-m2 ic50`, inclusion5000, count limit1. There is only one unique peptide per Event, so neither limit should truncate this pool. D1 is predeclared to evaluate an affinity-only common input without pretending that missing percentiles are supplied. Both output settings/statuses are retained, even if one fails.

## Keys and format conversion

Canonical genomic position is one-based4. For these SNVs, the actual installed pVACseq converter/report convention uses Start3/Stop4 (zero-based half-open); normalization maps Start+1 back to4. REF/ALT, versioned transcript, sequence, exact HLA and window1–8 must match. Source contexts are not counted as distinct synthesized peptides. pVAC transcript-set grouping is expanded using its native transcript list and per-HLA score arrays; standard fields are inspected before any custom-column/sidecar join. Sample is a single-sample file-level scope and is also explicitly carried in full input/output columns.

Origin IDs are Neoantimon2-specific identities. pVACtools is NOT required to generate them; both can be externally mapped to the common natural context key. Run identity, full phase state, native Grade, molecular Origin identity, compound contexts and unknown-phase semantics that are not a common tested input are NOT_COMPARABLE, not failures.

Q3/Q4 are EXTERNAL_QUERY over actual retained relations in BOTH tools, not native candidate exclusion/selection. No reduced object is injected into either tool. Q1/Q2 may combine NATIVE retrieval/export with an explicitly external normalization query. Q5 uses a new R process for public Neoantimon2 reopening and a new Python process for reading pVAC outputs; pVAC file read is EXTERNAL_QUERY, not a fabricated project-reopen API.

## Nonclaims

No model prediction, no VEP annotation claim, no pVACseq molecular reconstruction comparison, no all-candidate enumeration comparison, no clinical benefit or universal provenance advantage. Equivalent answers must be reported. Failure at this boundary is not tool-wide non-support. Aggregation is per Event; its row count is not a global unique peptide count.
