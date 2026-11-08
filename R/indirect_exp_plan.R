# Current IPD names include four-digit protein numbers; legacy global parsers stay frozen.
ie_pair <- function(raw) {
  scalar_text(raw,"explicit class II pair");chains<-strsplit(raw,"~",fixed=TRUE)[[1L]]
  if(length(chains)!=2L)stop("Explicit alpha~beta pair required")
  chains<-vapply(chains,ie_allele,"");genes<-vapply(chains,ie_gene,"")
  a<-which(genes%in%c("DRA","DQA1","DPA1"));b<-which(genes%in%c("DRB1","DRB3","DRB4","DRB5","DQB1","DPB1"))
  if(length(a)!=1L||length(b)!=1L||substr(genes[[a]],1,2)!=substr(genes[[b]],1,2))stop("Unresolved or incompatible Class II chain pairing")
  list(status="RESOLVED",hla_class="II",hla_alpha=chains[[a]],hla_beta=chains[[b]],mhc_key=paste(chains[[a]],chains[[b]],sep="~"))
}
ie_vendor <- function(pair) {
  p <- ie_pair(pair)
  if(p$status!="RESOLVED"||p$hla_class!="II")stop("Explicit Class II alpha~beta pair required")
  if(!grepl("\\*[0-9]{2,4}:[0-9]{2,4}$",p$hla_alpha)||!grepl("\\*[0-9]{2,4}:[0-9]{2,4}$",p$hla_beta))stop("Model presenters require exact two-field protein names; no truncation")
  if(startsWith(p$hla_beta,"HLA-DRB")) {
    if(p$hla_alpha!="HLA-DRA*01:01")stop("DR model requires explicit DRA*01:01")
    sub("^HLA-","",gsub(":","",sub("\\*","_",p$hla_beta)))
  } else paste0("HLA-",gsub("[*:]","",sub("^HLA-","",p$hla_alpha)),"-",gsub("[*:]","",sub("^HLA-","",p$hla_beta)))
}
ie_presenters <- function(declarations, typing, profile) {
  auto_dr<-is.null(declarations)
  if(is.null(declarations)) {
    if(profile!="dr-expanded-v1")stop("DQ/DP requires explicit confirmed alpha/beta pairing declarations; no Cartesian combinations")
    alpha <- Filter(function(t)t$locus=="DRA"&&t$status=="typed",typing)
    if(length(alpha)!=1L||!identical(unique(unlist(alpha[[1]]$alleles)),"HLA-DRA*01:01"))stop("DR auto-pairing requires a single explicitly typed DRA*01:01; otherwise declare pairs")
    beta <- Filter(function(t)t$locus%in%c("DRB1","DRB3","DRB4","DRB5")&&t$status=="typed",typing)
    declarations<-lapply(sort(unique(unlist(lapply(beta,`[[`,"alleles"))),method="radix"),function(b)list(pair=paste0("HLA-DRA*01:01~",b),pairing="explicit_DR_single_alpha",expression_confirmed=TRUE))
  }
  if(!is.list(declarations)||!is.null(names(declarations))||!length(declarations))stop("At least one explicit presenter declaration required")
  out<-lapply(declarations,function(x){
    allo_keys(x,c("pair","pairing","expression_confirmed"),"evidence","presenter declaration")
    p<-ie_pair(x$pair)
    if(p$status!="RESOLVED"||p$hla_class!="II")stop("Unresolved alpha/beta pair")
    group<-substr(ie_gene(p$hla_beta),1,2)
    if((profile=="dr-expanded-v1" && group!="DR")||(profile=="dqdp-exploratory-v1"&&!group%in%c("DQ","DP")))stop("Presenter outside selected profile")
    if(!x$pairing%in%c("confirmed_cis","confirmed_trans","explicit_expressed","explicit_DR_single_alpha"))stop("Pairing must be explicitly established, not inferred")
    if(group!="DR" && x$pairing=="explicit_DR_single_alpha")stop("DQ/DP cannot use inferred DR pairing")
    allo_boolean(x$expression_confirmed,"presenter expression")
    chains<-c(p$hla_alpha,p$hla_beta)
    typed<-all(vapply(chains,function(a)any(vapply(typing,function(t)t$status=="typed"&&a%in%unlist(t$alleles),TRUE)),TRUE))
    status<-if(!typed)"unresolved_typing" else if(!x$expression_confirmed)"expression_unconfirmed" else "eligible"
    list(presenter=p$mhc_key,group=group,vendor=ie_vendor(p$mhc_key),pairing=x$pairing,
      # Compatibility flag admits the presenter; it is not an experimental observation.
      expression_confirmed=x$expression_confirmed,
      expression_basis=if(auto_dr)"typing_based_assumption"else"user_declaration",
      experimental_expression_status="not_verified_by_workflow",
      status=status,reason=if(status=="eligible")NULL else status,
      evidence=if(is.null(x$evidence))"explicit hypothetical input, not patient observation" else allo_text(x$evidence,"pair evidence"))
  })
  ids<-vapply(out,`[[`,"","presenter");if(anyDuplicated(ids))stop("Duplicate presenter declaration")
  out[order(ids,method="radix")]
}
ie_windows <- function(typing,reference,loci,role) {
  ledger<-origins<-list();policy<-list(version="expanded-full-precursor-15mer-v1",loci=as.list(sort(loci)),stride=1L,length=15L,padding=FALSE,completion="none")
  gp<-indirect_id("expanded-generation-policy-v1",policy)
  for(t in Filter(function(t)t$locus%in%loci,typing)) {
    if(t$status!="typed") {
      ledger[[length(ledger)+1L]]<-list(role=role,locus=t$locus,status=t$status,reason=if(t$status=="absent")"explicit_absence_no_inference"else"typing_unresolved",origin_count=0L);next
    }
    for(a in sort(unique(unlist(t$alleles)),method="radix")) {
      rr<-Filter(function(r)r$allele==a,reference$records)
      r<-if(length(rr))rr[[1]]else NULL
      status<-if(is.null(r))"unavailable"else r$sequence_status
      if(grepl("[NLSCAQ]$",a))status<-"expression_uncertain"
      if(!is.null(r)&&!isTRUE(r$supported))status<-"unsupported"
      n<-0L
      if(status=="complete") {
        r$unknown_positions<-list()
        seqid<-allo_id("sequence",allo_sequence_identity_payload("NCBITaxon:9606",r))
        source<-indirect_id("expanded-hla-source-v1",list(role=role,locus=t$locus,allele=a,sequence_id=seqid,typing=t))
        g<-indirect_sequence_windows(r)
        if(!g$status%in%c("evaluated","not_applicable"))status<-"partial"
        for(o in g$windows)if(o$status=="evaluated") {
          pk<-biological_id("peptide_key",list(ligand_sequence=o$peptide_sequence),peptide_contract())
          id<-indirect_id("indirect-peptide-occurrence-id-v1",indirect_occurrence_payload(source,o$start,o$end,pk,gp))
          origins[[length(origins)+1L]]<-list(origin_id=id,source_occurrence_id=source,role=role,locus=t$locus,allele=a,
            sequence_id=seqid,peptide_key=pk,peptide_sequence=o$peptide_sequence,start=o$start,end=o$end,
            sequence_coordinates="precursor_1based_inclusive",typing=t,reference_sequence_identifier=r$sequence_identifier)
          n<-n+1L
        }
      }
      ledger[[length(ledger)+1L]]<-list(role=role,locus=t$locus,allele=a,status=status,reason=if(status=="complete")NULL else paste0("reference_or_expression_",status),origin_count=n)
    }
  }
  list(policy=policy,origins=origins,ledger=ledger,complete=all(vapply(ledger,function(l)l$status%in%c("complete","absent"),TRUE)))
}
ie_probe <- function(tool,model,root,timeout,mock) {
  h<-ia_exec(tool$path,"-h",root,"execution/help",timeout)
  if(!mock&&any(grepl("NEOANTIMON2_TEST_MOCK",h$output,fixed=TRUE)))stop("Mock requires mock_only authorization")
  text<-paste(h$output,collapse="\n")
  for(f in c("BA","context","inv_all"))if(!grepl(paste0("\\[-",f,"\\][[:space:]]+0([[:space:]]|$)"),text))stop("Cannot establish provider BA/context/inversion defaults")
  l<-ia_exec(tool$path,c("-rdir",model$rdir,"-list"),root,"execution/alleles",timeout)
  ia43_version(c(h$output,l$output))
  list(version="4.3k",tokens=unique(unlist(strsplit(paste(l$output,collapse=" "),"[[:space:],]+"))),commands=list(h$record,l$record))
}
ie_code_hash <- function() {
  e<-environment(ie_code_hash);ns<-sort(c(grep("^ie_",ls(e,all.names=TRUE),value=TRUE),"neo_indirect_class2"))
  infra_hash(list(expanded=lapply(ns,function(n){f<-get(n,e);list(name=n,formals=as.list(deparse(formals(f))),body=as.list(deparse(body(f))))}),reused=ia_code_hash()))
}
ie_guard <- function(root,min_free_gib) {
  df<-Sys.which("df");if(!nzchar(df))stop("Disk guard requires POSIX df")
  x<-system2(df,c("-Pk",shQuote(root)),stdout=TRUE)
  fields<-strsplit(trimws(tail(x,1L)),"[[:space:]]+")[[1L]]
  free<-suppressWarnings(as.numeric(fields[[4L]]))*1024
  if(!is.finite(free)||free<min_free_gib*1024^3)stop("Disk guard: insufficient free space; outputs retained")
  invisible(free)
}
