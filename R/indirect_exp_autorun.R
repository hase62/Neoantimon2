#' Run expanded research-only indirect HLA class II analysis
#'
#' Explicit source loci, exact presenters and scoped self comparisons. Neither
#' official PIRCHE reproduction nor a clinical risk score. Reuses the existing
#' 15mer enumeration and strict 4.3k BA parser; legacy autorun is unchanged.
#' @param donor,recipient Named locus lists. Character vectors declare typed
#'   allele copies; records can declare status typed/absent/untyped/uncertain,
#'   alleles, copy_number and evidence. DRB3/4/5 are never inferred from DRB1.
#' @param profile dr-expanded-v1 or separate dqdp-exploratory-v1.
#' @param installation List(path, authorization); path is the installed 4.3k
#'   root. authorization is list(status, terms), as in neo_indirect_autorun.
#' @param reference List(release, cache, optional download/commit) for official
#'   IPD, or a pinned local list(path, sha256) hla-reference-1.0.0 manifest.
#' @param output_dir New output directory; originals are never overwritten.
#' @param presenters Explicit list of pair/pairing/expression_confirmed records
#'   with optional evidence. Required for DQ/DP. DR can pair an explicitly typed
#'   single DRA*01:01 with explicitly typed DRB1/3/4/5.
#' @param source_loci,self_loci Explicit donor and recipient precursor scopes.
#'   Default is 13 loci; use the historical five explicitly for a reduced scope.
#' @param thresholds Named BA nM assumptions DR/DQ/DP (default 1000 each).
#' @param inversion Vendor default or all (explicit -inv_all).
#' @param dry_run Prepare reference, typing, windows and queries without a predictor.
#' @param resume Reopen saved results after exact request and hash verification.
#' @param batch_size Maximum peptides per local CPU invocation, at most 1024.
#' @param timeout Per-invocation timeout in seconds.
#' @param min_free_gib Required free-space reserve, at least 25 GiB.
#' @return Saved run, typed result and plan. Incomplete counts are NULL, with
#'   reasons and inspectable provisional counts. DQ/DP are never a pooled risk.
#' @export
neo_indirect_class2 <- function(donor, recipient, profile=c("dr-expanded-v1","dqdp-exploratory-v1"),
 installation=NULL, reference, output_dir, presenters=NULL, source_loci=ie_loci(), self_loci=source_loci,
 thresholds=list(DR=1000,DQ=1000,DP=1000), inversion=c("vendor_default","all"), dry_run=FALSE,resume=FALSE,
 batch_size=256L,timeout=600L,min_free_gib=25) {
  profile<-match.arg(profile);inversion<-match.arg(inversion);setup_flag(dry_run,"dry_run");setup_flag(resume,"resume")
  for(x in list(source_loci,self_loci))if(!is.character(x)||!length(x)||anyNA(x)||anyDuplicated(x)||any(!x%in%ie_loci()))stop("Explicit unique supported source/self loci required")
  source_loci<-sort(source_loci);self_loci<-sort(self_loci)
  allo_keys(thresholds,c("DR","DQ","DP"),label="separate BA assumptions")
  for(v in thresholds)if(!is.numeric(v)||length(v)!=1L||!is.finite(v)||v<=0)stop("Invalid BA threshold assumption")
  if(!is.numeric(min_free_gib)||length(min_free_gib)!=1L||!is.finite(min_free_gib)||min_free_gib<25)stop("Free-space reserve must be at least 25 GiB")
  for(v in list(batch_size,timeout))if(!is.numeric(v)||length(v)!=1L||!is.finite(v)||v<1||v!=floor(v))stop("Positive integer execution limits required")
  if(batch_size>1024)stop("Bounded CPU batches must not exceed 1024")
  dt<-ie_typing(donor,source_loci,"donor");rt<-ie_typing(recipient,sort(unique(c(self_loci,"DRA","DRB1","DRB3","DRB4","DRB5","DQA1","DQB1","DPA1","DPB1"))),"recipient")
  pp<-ie_presenters(presenters,rt,profile)
  presenter_scope_complete<-!is.null(presenters)||all(vapply(Filter(function(t)t$locus%in%c("DRA","DRB1","DRB3","DRB4","DRB5"),rt),function(t)t$status%in%c("typed","absent"),TRUE))
  presenter_scope<-if(!is.null(presenters))"explicit_declared_pairs_only; not all expressed molecules"else"all_explicitly_characterized_DRB1_3_4_5_with_single_DRA"
  alleles<-unique(c(unlist(lapply(c(dt,rt),`[[`,"alleles")),unlist(strsplit(vapply(pp,`[[`,"","presenter"),"~",fixed=TRUE))))
  root<-setup_dest(output_dir)
  if(!resume){
    existing<-function(x){while(!dir.exists(x))x<-dirname(x);x}
    ie_guard(existing(dirname(root)),min_free_gib)
    if(is.null(reference$path))ie_guard(existing(reference$cache),min_free_gib)
  }
  if(dir.exists(root)&&!resume)stop("Output directory exists; never overwritten")
  ref_source<-if(!is.null(reference$path)){
    if(!file.exists(reference$path)||ia_hash(reference$path)!=reference$sha256)stop("Reference SHA mismatch")
    list(path=normalizePath(reference$path),sha256=reference$sha256)
  }else{
    spec<-reference;if(resume)spec$download<-FALSE;ie_ipd_cache(spec)$source
  }
  tool<-model<-NULL;mock<-FALSE
  if(!dry_run){
    allo_keys(installation,c("path","authorization"),label="local provider installation")
    auth<-installation$authorization;allo_keys(auth,c("status","terms"),label="local authorization")
    if(!auth$status%in%c("authorized_local_installation","mock_only"))stop("Explicit authorized local installation required")
    scalar_text(auth$terms,"authorization terms");mock<-auth$status=="mock_only"
    exe<-setup_executable(file.path(installation$path,"netMHCIIpan"),"NetMHCIIpan")
    if(is.null(exe))stop("Authorized local 4.3k executable unavailable")
    tool<-list(path=exe,sha256=ia_hash(exe),authorization=auth)
    model<-ia43_model(list(path=installation$path,identity="explicit local 4.3k installation; physical tree pin"))
  }
  request<-list(version="expanded-indirect-class2-request-1.0.0",code_sha256=ie_code_hash(),donor=dt,recipient=rt,
    source_loci=as.list(source_loci),self_loci=as.list(self_loci),profile=profile,presenters=pp,reference_source=ref_source,
    tool=tool,model=if(is.null(model))NULL else model[c("path","rdir","tree_sha256")],thresholds=thresholds,
    inversion=inversion,dry_run=dry_run,batch_size=batch_size,timeout=timeout,min_free_gib=min_free_gib)
  id<-infra_hash(request)
  if(dir.exists(root)){
    if(!file.exists(file.path(root,"run.json")))stop("Existing interrupted directory retained")
    m<-jsonlite::read_json(file.path(root,"run.json"))
    if(!identical(m$request_id,id))stop("Resume request/tool/code/reference mismatch")
    if(!m$status%in%c("COMPLETED","INCOMPLETE","DRY_RUN"))stop("Failed run retained; choose a fresh directory")
    ia_verify_files(root,m);return(ia_return(root))
  }
  if(resume)stop("No saved run to reopen")
  dir.create(root,recursive=TRUE,showWarnings=FALSE);ie_guard(root,min_free_gib)
  lock<-file.path(root,".autorun-lock");if(!dir.create(lock))stop("Output locked")
  writeLines(as.character(Sys.getpid()),file.path(lock,"owner"));on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  manifest<-list(version="expanded-indirect-class2-run-1.0.0",request_id=id,request=request,status="RUNNING",real_predictor_executed=FALSE,
    official_pirche_reproduction=FALSE,clinical_risk_score=NULL,execution_kind=if(dry_run)"NOT_EXECUTED"else if(mock)"MOCK"else"AUTHORIZED_LOCAL_4.3k")
  ia_json(manifest,file.path(root,"run.json"))
  tryCatch({
    ia_json(list(donor=donor,recipient=recipient,presenters=presenters,profile=profile,reference=reference),file.path(root,"inputs.json"))
    ref<-ie_reference(reference,alleles);ia_json(ref$manifest,file.path(root,"reference.json"));ia_json(ref$provenance,file.path(root,"reference-provenance.json"))
    g<-ie_windows(dt,ref$manifest,source_loci,"donor");s<-ie_windows(rt,ref$manifest,self_loci,"recipient");origins<-c(g$origins,s$origins)
    for(i in seq_along(pp)){
      chains<-strsplit(pp[[i]]$presenter,"~",fixed=TRUE)[[1]]
      good<-all(vapply(chains,function(a)any(vapply(ref$manifest$records,function(r)r$allele==a&&r$sequence_status=="complete"&&isTRUE(r$supported),TRUE)),TRUE))
      if(!good&&pp[[i]]$status=="eligible"){pp[[i]]$status<-"reference_unavailable";pp[[i]]$reason<-"presenter_chain_reference_incomplete"}
    }
    links<-list();for(p in pp)for(o in origins)links[[length(links)+1L]]<-list(query_id=indirect_id("indirect-autorun-query-v1",list(peptide_sequence=o$peptide_sequence,presenter=p$presenter)),
      peptide_sequence=o$peptide_sequence,presenter=p$presenter,origin_id=o$origin_id,source_occurrence_id=o$source_occurrence_id,role=o$role,locus=o$locus)
    queries<-allo_records(lapply(links,function(x)x[c("query_id","peptide_sequence","presenter")]))
    plan<-list(profile=profile,presenter_scope=presenter_scope,presenter_scope_complete=presenter_scope_complete,source_loci=as.list(source_loci),self_loci=as.list(self_loci),presenters=pp,queries=queries,source_complete=g$complete,self_complete=s$complete)
    ia_json(plan,file.path(root,"plan.json"));ie_table(dt,file.path(root,"donor_typing.tsv"));ie_table(rt,file.path(root,"recipient_typing.tsv"));ie_table(origins,file.path(root,"origins.tsv"));ie_table(links,file.path(root,"query_origin_links.tsv"));ie_table(c(g$ledger,s$ledger),file.path(root,"coverage.tsv"));ie_table(queries,file.path(root,"queries.tsv"))
    rows<-commands<-list();batch<-0L
    if(!dry_run){
      probe<-ie_probe(tool,model,root,timeout,mock);commands<-probe$commands
      ia_json(model,file.path(root,"execution/model-tree.json"))
      for(i in seq_along(pp)){
        p<-pp[[i]]
        if(p$status!="eligible")next
        if(!p$vendor%in%probe$tokens){pp[[i]]$status<-"unsupported";pp[[i]]$reason<-"exact_model_allele_absent; no fuzzy substitution";next}
        peptides<-sort(unique(vapply(queries[vapply(queries,function(q)q$presenter==p$presenter,TRUE)],`[[`,"","peptide_sequence")),method="radix")
        batches<-if(length(peptides))split(peptides,ceiling(seq_along(peptides)/batch_size))else list()
        for(peps in batches){
          ie_guard(root,min_free_gib);batch<-batch+1L;stem<-sprintf("raw/batch-%05d",batch);input<-file.path(root,paste0(stem,".pep"));dir.create(dirname(input),showWarnings=FALSE);writeLines(peps,input)
          args<-c("-rdir",model$rdir,"-inptype","1","-f",input,"-a",p$vendor,"-BA",if(inversion=="all")"-inv_all"else character())
          ran<-ia_exec(tool$path,args,root,stem,timeout);commands<-c(commands,list(ran$record))
          parsed<-ia43_parse(ran$output,peps,p$presenter,vendor=p$vendor)
          for(j in seq_along(parsed)){parsed[[j]]$raw_output$artifact<-paste0(stem,".stdout");parsed[[j]]$evidence_id<-indirect_id("expanded-prediction-evidence-v1",list(row=parsed[[j]],executable=tool$sha256,model=model$tree_sha256))}
          rows<-c(rows,parsed)
        }
        pp[[i]]$status<-"evaluated"
      }
      if(ia_hash(tool$path)!=tool$sha256||ia43_model(list(path=model$path,identity=model$identity))$tree_sha256!=model$tree_sha256)stop("Provider changed during execution")
      manifest$real_predictor_executed<-!mock&&length(rows)>0L
      manifest$provider_version<-probe$version
    }
    plan$presenters<-pp;ia_json(plan,file.path(root,"plan.json"));ia_json(commands,file.path(root,"execution/commands.json"))
    ia_json(rows,file.path(root,"predictions.json"));ie_table(rows,file.path(root,"predictions.tsv"));if(length(rows))ie_table(ia43_metrics_table(rows),file.path(root,"model_metrics.tsv"))
    derived<-ie_result(origins,pp,rows,g$complete,s$complete,profile,thresholds)
    # Explicitly report configured genes even when they have no eligible windows.
    represented<-unique(vapply(g$origins,`[[`,"","locus"))
    for(locus in setdiff(source_loci,represented))for(p in pp){
      ledger<-Filter(function(l)l$locus==locus,g$ledger)
      known_empty<-length(ledger)>0L&&all(vapply(ledger,function(l)l$status%in%c("complete","absent"),TRUE))
      derived$locus_summary[[length(derived$locus_summary)+1L]]<-list(locus=locus,presenter=p$presenter,group=p$group,
        status=if(known_empty)"known_no_eligible_windows"else"incomplete_source",complete_count=if(known_empty)0L else NULL,
        conditional_count=0L,attribution="no windows; inspect coverage for absence versus unknown")
    }
    if(!presenter_scope_complete)for(i in seq_along(derived$group_summary)){
      derived$group_summary[[i]]["complete_count"]<-list(NULL)
      derived$group_summary[[i]]$reason<-"untyped_or_uncertain_DR_presenter_locus; inspect recipient typing"
    }
    complete<-g$complete&&s$complete&&presenter_scope_complete&&all(vapply(pp,function(p)p$status=="evaluated",TRUE))
    result<-c(list(schema_version=if(profile=="dr-expanded-v1")"indirect-class2-dr-expanded-1.0.0"else"indirect-class2-dqdp-exploratory-1.0.0",
      profile=profile,provider_version=if(dry_run)NULL else"4.3k",status=if(dry_run)"DRY_RUN"else if(complete)"COMPLETED"else"INCOMPLETE",
      source_complete=g$complete,self_complete=s$complete,presenter_scope=presenter_scope,presenter_scope_complete=presenter_scope_complete,official_pirche_reproduction=FALSE,clinical_risk_score=NULL,
      thresholds=thresholds,counting_policy="unique(core,orientation,exact_presenter); source-locus attributions overlap",
      self_scope=as.list(self_loci),source_scope=as.list(source_loci),donor_origin_count=length(g$origins),self_origin_count=length(s$origins),unique_query_count=length(queries),
      limitations=list("SOT donor-to-recipient only","no full human proteome self search","BA thresholds are assumptions, not calibrated clinical probabilities","DQ/DP is exploratory; no pooled clinical score","no official PIRCHE/Frost/PIRCHE-I/B-cell equivalence")),derived)
    ia_json(result,file.path(root,"result.json"));ia_json(list(result=result,coverage=c(g$ledger,s$ledger),origins=origins,query_origin_links=links,reference_provenance=ref$provenance),file.path(root,"explanation.json"))
    for(n in names(derived))ie_table(derived[[n]],file.path(root,paste0(n,".tsv")))
    manifest$status<-result$status;manifest$source_complete<-g$complete;manifest$self_complete<-s$complete;manifest$presenter_scope_complete<-presenter_scope_complete
    manifest$execution_settings<-list(BA=TRUE,context_encoding=FALSE,inversion=inversion,thresholds=thresholds,threads=1L,temporary_policy="vendor default",calibration="none; explicit BA assumptions")
    if(!is.null(reference$path)){if(ia_hash(reference$path)!=reference$sha256)stop("Reference changed during execution")}else{spec<-reference;spec$download<-FALSE;if(!identical(canonical_bytes(ie_ipd_cache(spec)$source),canonical_bytes(ref_source)))stop("IPD reference changed during execution")}
    ia_finish(root,manifest);ia_return(root)
  },error=function(e){manifest$status<-"FAILED";manifest$error<-conditionMessage(e);ia_finish(root,manifest);stop(e)})
}
