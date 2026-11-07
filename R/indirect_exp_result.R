ie_table <- function(rows,path) {
  if(!length(rows)) {dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);writeLines(character(),path);invisible(path)} else ia_table(rows,path)
}
ie_key <- function(peptide,presenter) paste(peptide,presenter,sep="|")
ie_result <- function(origins,presenters,rows,source_complete,self_complete,profile,thresholds) {
  donor<-Filter(function(o)o$role=="donor",origins);self<-Filter(function(o)o$role=="recipient",origins)
  peps<-sort(unique(vapply(donor,`[[`,"","peptide_sequence")),method="radix")
  lookup<-list();for(r in rows)lookup[[ie_key(r$peptide_sequence,r$presenter)]]<-r
  candidates<-links<-units<-summaries<-list()
  unit_key<-function(r)paste(r$cores[[1]]$core_sequence,r$cores[[1]]$orientation,r$presenter,sep="|")
  for(p in presenters) {
    cutoff<-thresholds[[p$group]]
    self_peps<-unique(vapply(self,`[[`,"","peptide_sequence"))
    hits<-unique(vapply(Filter(function(r)r$presenter==p$presenter&&r$peptide_sequence%in%self_peps&&r$affinity_nM<cutoff,rows),unit_key,""))
    supported<-p$status=="evaluated";keys<-character();evaluated<-0L;excluded<-0L;unknown<-0L
    for(pep in peps) {
      oo<-Filter(function(o)o$peptide_sequence==pep,donor)
      id<-biological_id("candidate_key",list(species="NCBITaxon:9606",hla_class="II",ligand_sequence=pep,mhc_key=p$presenter),peptide_contract())
      r<-lookup[[ie_key(pep,p$presenter)]];state<-if(is.null(r))p$status else "evaluated"
      decision<-if(is.null(r))"not_evaluated"else if(r$affinity_nM>=cutoff)"excluded_nonbinder"else if(unit_key(r)%in%hits)"self_match"else if(self_complete&&supported)"supported_nonself_within_scope"else "self_search_incomplete"
      if(!is.null(r))evaluated<-evaluated+1L
      if(decision=="self_match")excluded<-excluded+1L
      if(decision%in%c("not_evaluated","self_search_incomplete"))unknown<-unknown+1L
      if(!is.null(r)&&r$affinity_nM<cutoff&&!unit_key(r)%in%hits)keys<-c(keys,unit_key(r))
      count_id<-if(!is.null(r))indirect_id("expanded-core-presenter-v1",list(profile=profile,provider_version="4.3k",threshold_nM=cutoff,core=r$cores[[1]]$core_sequence,orientation=r$cores[[1]]$orientation,presenter=p$presenter))else NULL
      cr<-list(candidate_id=id,peptide_sequence=pep,peptide_length=15L,peptide_key=oo[[1]]$peptide_key,presenter=p$presenter,group=p$group,prediction_status=state,
        decision=decision,affinity_nM=if(is.null(r))NULL else r$affinity_nM,rank_BA=if(is.null(r))NULL else r$percentile_rank,
        core=if(is.null(r))NULL else r$cores[[1]],counted_unit_id=count_id,threshold_nM=cutoff,origin_ids=as.list(vapply(oo,`[[`,"","origin_id")))
      candidates[[length(candidates)+1L]]<-cr
      for(o in oo)links[[length(links)+1L]]<-list(candidate_id=id,origin_id=o$origin_id,source_occurrence_id=o$source_occurrence_id,peptide_key=o$peptide_key,presenter=p$presenter,locus=o$locus,role="donor")
      if(decision%in%c("supported_nonself_within_scope","self_search_incomplete"))units[[length(units)+1L]]<-list(counted_unit_id=count_id,candidate_id=id,presenter=p$presenter,group=p$group,core=r$cores[[1]]$core_sequence,orientation=r$cores[[1]]$orientation,origin_ids=cr$origin_ids,status=decision)
    }
    complete<-source_complete&&self_complete&&supported
    summaries[[length(summaries)+1L]]<-list(presenter=p$presenter,group=p$group,status=if(complete)"complete"else"incomplete",prediction_status=p$status,reason=p$reason,
      complete_count=if(complete)length(unique(keys))else NULL,conditional_count=length(unique(keys)),conditional_semantics="provisional_nonself_relative_to_observed_self; not a complete score",
      evaluated_donor_pmhcs=evaluated,self_excluded_pmhcs=excluded,unknown_pmhcs=unknown,threshold_nM=cutoff,threshold_semantics="explicit BA assumption; not cross-model calibration")
  }
  by_locus<-list()
  for(locus in sort(unique(vapply(donor,`[[`,"","locus"))))for(p in presenters){
    ids<-vapply(Filter(function(o)o$locus==locus,donor),`[[`,"","origin_id")
    uu<-Filter(function(u)u$presenter==p$presenter&&any(unlist(u$origin_ids)%in%ids),units)
    summary<-summaries[[match(p$presenter,vapply(summaries,`[[`,"","presenter"))]]
    n<-length(unique(vapply(uu,`[[`,"","counted_unit_id")))
    by_locus[[length(by_locus)+1L]]<-list(locus=locus,presenter=p$presenter,group=p$group,complete_count=if(summary$status=="complete")n else NULL,conditional_count=n,attribution="overlapping union; do not sum loci")
  }
  group_summary<-lapply(sort(unique(vapply(presenters,`[[`,"","group"))),function(g){
    s<-Filter(function(x)x$group==g,summaries);u<-Filter(function(x)x$group==g,units);n<-length(unique(vapply(u,`[[`,"","counted_unit_id")))
    list(group=g,complete_count=if(all(vapply(s,function(x)x$status=="complete",TRUE)))n else NULL,conditional_count=n,
      unit="unique(core,orientation,exact_presenter)",clinical_risk_score=NULL)
  })
  list(candidates=candidates,candidate_origin_links=links,units=units,presenter_summary=summaries,locus_summary=by_locus,group_summary=group_summary)
}
