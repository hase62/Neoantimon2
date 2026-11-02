# Fast offline numerical example. ALL scores and references below are SYNTHETIC.
# No predictor, download, neo_run(), or existing project write is performed.
library(Neoantimon2)
experimental_yield_example <- function() {
  aa <- strsplit("ACDEFGHIKLMNPQRSTVWY","")[[1]]
  x <- data.frame(candidate_id=paste0("synthetic-",seq_len(20)),sample_id="synthetic-sample",
    peptide_sequence=paste0("ACDEFGHI",aa),mhc_key="HLA-A*02:01",
    raw_score=seq(0,1,length.out=20),prediction_status="evaluated",source="External",
    provider="synthetic-IM",predictor_version="synthetic-1",model_version="synthetic-1",
    run_id="synthetic-run",task="external:immunogenicity",metric="score",
    metric_kind="raw_score",metric_semantics="raw",unit="unitless",score_direction="higher_is_better")
  x$origin_ids <- I(lapply(seq_len(nrow(x)),function(i)list(paste0("synthetic-origin-",i))))
  profile <- list(schema_version="experimental-yield-profile-0.1.0",id="user-assumptions",
    version="1",mode="scenario",endpoint="example-T-cell-endpoint",assay="example-assay",
    population="Synthetic planning illustration",limitations="No measured experimental probabilities",
    provenance=list(source="User-selected assumed rates, not an empirical dataset"),
    binning=list(basis="fixed_thresholds",thresholds=c(.3,.6,.9),score_direction="higher_is_better"),
    rates=c(.01,.02,.03,.06))
  plan <- neo_experimental_yield(x,profile,unit="peptide")
  stopifnot(plan$summary$N_available==20,all(is.na(plan$units$evidence_category)),
    all(diff(plan$curve$probability_at_least_one)>=0))
  print(plan$summary);print(plan$targets)
  plan
}
if (sys.nframe()==0L) experimental_yield_example()
