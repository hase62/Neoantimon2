# Synthetic graph, independent of production reconstruction and score calculation.
report_test_fixture <- function() {
  event <- function(id, role="somatic", sample="s1") list(event_id=id,sample_id=sample,event_role=role,
    event_type="SNV",contig="chr1",pos=10L,ref="A",alt="T")
  origin <- function(id, events, pep="ACDEFGHIK", tx=id, sample="s1") list(origin_id=id,sample_id=sample,
    peptide_sequence=pep,peptide_length=nchar(pep),source_event_ids=as.list(events),
    germline_context_ids=list("g1"),gene_symbol="GENE",transcript_id=tx,transcript_version="1",
    genes_json=list(list(gene_id="ENSG1",gene_symbol="GENE")),nested=list(null=NULL,empty=list(),truth=TRUE))
  oo <- list(origin("o1",c("e1","e2")),origin("o2","e2"),origin("o3","e3"),
    origin("o4","e1","ACDEFGHIKLMNPQR"),origin("orphan","e4","CDEFGHIKL"))
  candidate <- function(id,oids,hla="HLA-A*02:01",pep="ACDEFGHIK",cl="I")
    list(candidate_id=id,sample_id="s1",peptide_sequence=pep,peptide_length=nchar(pep),
      mhc_key=hla,hla_class=cl,analysis_context="somatic",origin_ids=as.list(oids))
  cc <- list(candidate("c1",c("o1","o2")),candidate("c2",c("o1","o2"),"HLA-B*07:02"),
    candidate("c3","o3"),candidate("c4","o4","HLA-DRB1*01:01","ACDEFGHIKLMNPQR","II"))
  links <- unlist(lapply(cc,function(c)lapply(c$origin_ids,function(oid)list(candidate_id=c$candidate_id,
    origin_id=oid,sample_id=c$sample_id))),recursive=FALSE)
  native <- lapply(seq_along(cc),function(i)list(candidate_id=cc[[i]]$candidate_id,sample_id="s1",
    final_grade=if(i==3)"NOT_GRADED" else "B",primary_grade="B",rank_within_sample_class=as.numeric(c(2,1,3,1)[i]),
    rank_status="RANKED",policy_version="native-policy",priority_origin_id=c("o2","o2","o3","o4")[i],
    primary_prediction_status=if(i==3)"not_run" else "evaluated"))
  external <- function(id,direction="higher_is_better",run="r1",model="model1",unit="unitless",scores=c(.8,.8,.1,.4)) {
    ee<-lapply(seq_along(cc),function(i)list(candidate_id=cc[[i]]$candidate_id,origin_ids=cc[[i]]$origin_ids,
      provider="synthetic-model",predictor_version="1",model_version=model,run_id=run,
      task="external:immunogenicity",metric="score",metric_kind="raw_score",metric_semantics="raw",unit=unit,
      score_direction=direction,scope="peptide_hla",mhc_key=cc[[i]]$mhc_key,status="evaluated",
      numeric_score=scores[i],score=as.character(scores[i]),evidence_id=paste0(id,i),reason=""))
    list(manifest=list(analysis_id=id),payload=list(evidence=ee))
  }
  list(graph=list(events=list(event("e1"),event("e2"),event("e3"),event("e4"),event("empty"),event("g1","germline")),
    origins=oo,candidates=cc,links=links,annotations=list()),native=native,
    external=list(external("high"),external("low","lower_is_better","r2","model2","nM",c(10,5,20,3))))
}
report_test_build <- function(f=report_test_fixture(),...) Neoantimon2:::report_construct(
  f$graph,f$native,external=f$external,...)
