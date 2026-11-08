# Expanded workflow only. No change to the frozen allo/3.2 reference resolver.
ie_loci <- function() c("A","B","C","DRB1","DRB3","DRB4","DRB5","DQA1","DQB1","DPA1","DPB1","MICA","MICB")
ie_allele <- function(x, locus = NULL) {
  scalar_text(x,"explicit allele")
  x <- toupper(trimws(x)); x <- sub("^HLA-","",x)
  if (!grepl("^[A-Z0-9]+\\*[0-9]{2,4}(:[0-9]{2,4}){0,3}[NLSCAQ]?$",x)) stop("Invalid exact allele notation; no fuzzy resolution")
  gene <- sub("\\*.*","",x)
  if (!gene %in% c(ie_loci(),"DRA") || (!is.null(locus) && gene != locus)) stop("Allele/locus conflict")
  paste0("HLA-",x)
}
ie_gene <- function(x) sub("\\*.*","",sub("^HLA-","",x))
ie_typing <- function(x, loci, role) {
  if (!is.list(x) || (length(x)>0L && is.null(names(x))) || anyDuplicated(names(x)) || any(!names(x) %in% c(ie_loci(),"DRA"))) stop("Typing must be a named locus list")
  lapply(loci,function(locus) {
    t <- x[[locus]]
    if (is.null(t)) t <- list(status="untyped",alleles=character())
    if (is.character(t)) t <- list(status="typed",alleles=t,copy_number=length(t))
    allo_keys(t,"status",c("alleles","copy_number","evidence"),"expanded typing")
    if (!t$status %in% c("typed","absent","untyped","uncertain")) stop("Invalid typing status")
    if (is.null(t$alleles)) t$alleles <- character()
    if (!is.character(t$alleles) || anyNA(t$alleles)) stop("Alleles must be explicit character values")
    a <- unname(vapply(t$alleles,ie_allele,"",locus=locus))
    if (t$status=="typed" && !length(a)) stop("Typed locus requires an allele")
    if (t$status %in% c("absent","untyped") && length(a)) stop("Absent/untyped locus cannot supply alleles")
    n <- t$copy_number
    if (!is.null(n) && (!is.numeric(n)||length(n)!=1L||is.na(n)||n!=floor(n)||n<0||
        (t$status=="absent" && n!=0)|| (t$status=="typed" && n<length(a)))) stop("Invalid explicitly declared copy number")
    list(role=role,locus=locus,status=t$status,alleles=as.list(a),copy_number=if(t$status=="absent")0L else n,
      evidence=if(is.null(t$evidence))"user declaration; no DRB1-based inference" else allo_text(t$evidence,"typing evidence"))
  })
}
ie_ipd_cache <- function(spec) {
  allo_keys(spec,c("release","cache"),c("download","commit"),"IPD reference specification")
  release <- spec$release
  if (!is.character(release)||length(release)!=1L||!grepl("^3\\.[0-9]{2}(\\.[0-9])?$",release)) stop("Explicit IPD release required")
  if (length(strsplit(release,".",fixed=TRUE)[[1]])==2L) release <- paste0(release,".0")
  root <- file.path(spec$cache,release); pin <- file.path(root,"source.json")
  files <- c("LICENCE.md","release_version.txt","Allele_status.txt","hla_prot.fasta")
  if (file.exists(pin)) {
    source <- jsonlite::read_json(pin)
    if (!identical(source$release,release)||!grepl("^[0-9a-f]{40}$",source$commit)||!setequal(names(source$files),files)) stop("Invalid IPD cache pin")
    if (!is.null(spec$commit) && spec$commit!=source$commit) stop("IPD cache commit mismatch")
    for (f in files) if (!file.exists(file.path(root,f)) || ia_hash(file.path(root,f))!=source$files[[f]]$sha256) stop("IPD cache integrity mismatch: ",f)
    if (!any(grepl(paste0("IPD-IMGT/HLA ",release),readLines(file.path(root,"release_version.txt"),warn=FALSE),fixed=TRUE))) stop("IPD release content mismatch")
    return(list(root=normalizePath(root),source=source))
  }
  if (identical(spec$download,FALSE)) stop("Pinned IPD cache unavailable; downloads disabled")
  if (dir.exists(root) && length(list.files(root,all.files=TRUE,no..=TRUE))) stop("Incomplete IPD cache retained; choose a fresh cache")
  dir.create(root,recursive=TRUE,showWarnings=FALSE)
  commit <- spec$commit
  if (is.null(commit)) {
    if (release=="3.66.0") commit <- "b2055ea33ef464781b9e8ecab20b275b6704f878" else {
      api <- paste0("https://api.github.com/repos/ANHIG/IMGTHLA/commits/",gsub(".","",release,fixed=TRUE))
      commit <- jsonlite::fromJSON(paste(readLines(url(api),warn=FALSE),collapse="\n"))$sha
    }
  }
  if (!is.character(commit)||length(commit)!=1L||!grepl("^[0-9a-f]{40}$",commit)) stop("Invalid immutable IPD commit")
  source <- list(release=release,commit=commit,files=list())
  for (f in files) {
    u <- paste0("https://raw.githubusercontent.com/ANHIG/IMGTHLA/",commit,"/",f)
    target <- file.path(root,f)
    if (utils::download.file(u,target,mode="wb",quiet=TRUE,method="libcurl")!=0L) stop("IPD download failed; incomplete cache preserved")
    source$files[[f]] <- list(url=u,sha256=ia_hash(target),size=unname(file.info(target)$size))
  }
  if (!any(grepl(paste0("IPD-IMGT/HLA ",release),readLines(file.path(root,"release_version.txt"),warn=FALSE),fixed=TRUE))) stop("Downloaded release content disagrees")
  ia_json(source,pin);list(root=normalizePath(root),source=source)
}
ie_ipd_records <- function(cache, alleles) {
  lines <- readLines(file.path(cache$root,"hla_prot.fasta"),warn=FALSE)
  headers <- which(startsWith(lines,">")); ends <- c(headers[-1L]-1L,length(lines))
  names <- vapply(strsplit(lines[headers],"[[:space:]]+"),`[[`,"",2L)
  seqs <- vapply(seq_along(headers),function(i)paste(lines[seq.int(headers[[i]]+1L,ends[[i]])],collapse=""),"")
  if (anyDuplicated(names)) stop("Duplicate IPD allele FASTA identity")
  status <- utils::read.csv(file.path(cache$root,"Allele_status.txt"),comment.char="#",stringsAsFactors=FALSE,check.names=FALSE)
  if (anyDuplicated(status$Allele)) stop("Duplicate IPD completeness record")
  full <- function(i) {
    st <- status[match(names[[i]],status$Allele),,drop=FALSE]
    n <- nchar(seqs[[i]])
    isTRUE(nrow(st)==1L && !is.na(st$Partial) && st$Partial=="Full" && st$Start==1L &&
      st$End==3L*(n+1L) && grepl("^[ACDEFGHIKLMNPQRSTVWY]+$",seqs[[i]]))
  }
  lapply(sort(unique(alleles),method="radix"),function(a) {
    raw <- sub("^HLA-","",a);hit <- which(names==raw);operation <- "exact_allele"
    # Two-field HLA means a protein name, never a guessed full genotype.
    # A fully specified partial allele must NOT fall back to another allele.
    protein_name <- (grepl("^[A-Z0-9]+\\*[0-9]{2,4}:[0-9]{2,4}$",raw) && !ie_gene(a)%in%c("MICA","MICB")) || (ie_gene(a)%in%c("MICA","MICB") && grepl("^[A-Z]+\\*[0-9]{3,4}$",raw))
    if (!length(hit) && protein_name) {
      candidates <- which(startsWith(names,paste0(raw,":")) & !grepl("[NLSCAQ]$",names))
      hit <- candidates[vapply(candidates,full,TRUE)]
      operation <- "explicit_two_field_protein_reference; identical_complete_synonymous_records; not genotype imputation"
      if (length(unique(seqs[hit]))>1L) hit <- integer()
    }
    seq <- if(length(hit))seqs[[hit[[1L]]]] else NULL
    complete <- length(hit)>0L && all(vapply(hit,full,TRUE))
    n <- if(is.null(seq))0L else nchar(seq)
    rec <- list(allele=a,sequence_identifier=paste0("IPD:",cache$source$release,":",raw),protein_sequence=seq,
      sequence_status=if(is.null(seq))"unavailable" else if(complete)"complete" else "partial",supported=TRUE,
      protein_definition="IPD_literal_precursor",isoform="submitted_protein",coordinate_map=list(sequence_coordinates="precursor",precursor_start=if(complete)1L else NULL,mature_start=NULL),
      domains=if(complete)list(list(name="full_precursor",start=1L,end=n))else list(),
      covered_regions=if(complete)list(list(name="full_precursor",start=1L,end=n))else list(),
      resolution_provenance=paste0(operation,"; ",jsonlite::toJSON(list(source_alleles=as.list(names[hit]),commit=cache$source$commit,files=cache$source$files),auto_unbox=TRUE)))
    # Partial FASTA endpoints are not assumed to be complete precursor coordinates.
    rec
  })
}
ie_reference <- function(spec, alleles) {
  if (!is.null(spec$path)) {
    allo_keys(spec,c("path","sha256"),label="reference binding")
    if (!file.exists(spec$path)||ia_hash(spec$path)!=spec$sha256) stop("Reference SHA mismatch")
    m <- jsonlite::read_json(spec$path); if(m$version!="hla-reference-1.0.0"||m$sequence_policy!="precursor"||m$completion_policy!="none")stop("Unsupported reference contract")
    provenance <- list(source_path=normalizePath(spec$path),sha256=spec$sha256)
  } else {
    c <- ie_ipd_cache(spec);records <- ie_ipd_records(c,alleles)
    m <- list(version="hla-reference-1.0.0",status="evaluated",reference_source="IPD-IMGT/HLA",reference_release=c$source$release,
      nomenclature_authority="WHO/IPD-IMGT/HLA",nomenclature_release=c$source$release,species="NCBITaxon:9606",
      provider=list(name="ANHIG/IMGTHLA",version=c$source$commit),parser=list(name="literal_amino_acids",version="1.0.0"),coordinate_convention="one_based_inclusive",
      isoform_policy="explicit_mapped_isoforms",completion_policy="none",sequence_policy="precursor",records=records,
      provenance=list(source=c$source,rights="IPD Attribution-NoDerivs; modified data redistribution requires permission; local generated manifest only"))
    provenance <- c$source
  }
  if (!is.list(m$records)) stop("Missing reference records")
  if(m$status!="evaluated"||m$species!="NCBITaxon:9606"||m$coordinate_convention!="one_based_inclusive"||m$parser$name!="literal_amino_acids"||m$parser$version!="1.0.0")stop("Invalid reference status/species/parser/coordinates")
  seen <- character()
  for(i in seq_along(m$records)) {
    r <- m$records[[i]];r$allele <- ie_allele(r$allele)
    if(r$allele%in%seen)stop("Duplicate exact reference mapping");seen<-c(seen,r$allele)
    if(!r$sequence_status%in%c("complete","partial","unavailable"))stop("Invalid reference status")
    if(!is.null(r$protein_sequence) && !grepl("^[ACDEFGHIKLMNPQRSTVWYX]+$",r$protein_sequence))stop("Invalid literal reference sequence")
    if(r$sequence_status=="unavailable" && !is.null(r$protein_sequence))stop("Unavailable reference has residues")
    if(r$sequence_status=="complete") {
      if(is.null(r$protein_sequence)||grepl("X",r$protein_sequence))stop("Invalid complete reference")
      covered<-unique(unlist(lapply(r$covered_regions,function(z)seq.int(z$start,z$end))))
      if(!setequal(covered,seq_len(nchar(r$protein_sequence))))stop("Incomplete full precursor coverage")
    }
    m$records[[i]]<-r
  }
  list(manifest=m,provenance=provenance)
}
