# Indexed BAM queries use Rsamtools/HTSlib; allele classification reconstructs
# a local haplotype from CIGAR, never a position-only indel pileup.
input_bam_options <- function(options=NULL) {
  defaults <- list(min_mapq=20L,min_baseq=20L,include_duplicates=FALSE,
    mode="fallback",sample_name=NULL)
  if (is.null(options)) return(defaults)
  strict_keys(options,names(defaults),"BAM options")
  for (k in names(options)) defaults[k] <- options[k]
  for (k in c("min_mapq","min_baseq")) {
    defaults[[k]] <- somatic_integer(defaults[[k]],k)
    if (is.null(defaults[[k]]) || defaults[[k]]>255L) stop("Invalid BAM quality threshold")
  }
  if (!is.logical(defaults$include_duplicates) || length(defaults$include_duplicates)!=1L ||
      is.na(defaults$include_duplicates)) stop("Invalid duplicate policy")
  if (!defaults$mode %in% c("fallback","independent")) stop("BAM mode must be fallback or independent")
  if (!is.null(defaults$sample_name)) scalar_text(defaults$sample_name,"BAM sample_name")
  defaults
}

input_bam_haplotype <- function(pos,cigar,sequence,qualities,start,end,min_baseq) {
  tokens <- regmatches(cigar,gregexpr("[0-9]+[MIDNSHP=X]",cigar))[[1L]]
  if (!length(tokens) || paste(tokens,collapse="")!=cigar) return(NULL)
  bases <- strsplit(sequence,"",fixed=TRUE)[[1L]]
  if (nchar(qualities)!=length(bases) || qualities=="*") return(NULL)
  qual <- utf8ToInt(qualities)-33L
  covered <- rep(FALSE,end-start+1L); hap <- rep("",end-start+1L)
  rp <- pos; qp <- 1L
  for (token in tokens) {
    op <- substring(token,nchar(token)); len <- as.integer(substr(token,1L,nchar(token)-1L))
    if (op %in% c("M","=","X")) {
      lo <- max(start,rp); hi <- min(end,rp+len-1L)
      if (lo<=hi) {
        ri <- seq.int(lo,hi); qi <- qp+ri-rp
        if (any(qi>length(bases)) || any(qual[qi]<min_baseq)) return(NULL)
        covered[ri-start+1L] <- TRUE; hap[ri-start+1L] <- bases[qi]
      }
      rp <- rp+len; qp <- qp+len
    } else if (op=="I") {
      anchor <- rp-1L
      if (anchor>=start && anchor<=end) {
        qi <- seq.int(qp,qp+len-1L)
        if (any(qi>length(bases)) || any(qual[qi]<min_baseq)) return(NULL)
        at <- anchor-start+1L; hap[at] <- paste0(hap[at],paste(bases[qi],collapse=""))
      }
      qp <- qp+len
    } else if (op=="D") {
      lo <- max(start,rp);hi <- min(end,rp+len-1L)
      if (lo<=hi) covered[seq.int(lo,hi)-start+1L] <- TRUE
      rp <- rp+len
    } else if (op=="N") {
      # An intron is missing sequence, never a REF base or a deletion allele.
      rp <- rp+len
    } else if (op=="S") qp <- qp+len
  }
  if (!all(covered)) NULL else paste(hap,collapse="")
}

# Target-specific indel edits are normalized with the same reference/VCF
# normalizer as Events. Unrelated aligned substitutions do not negate an indel.
# Full window coverage/quality remains required; N is never a deletion.
input_bam_indel_call <- function(pos,cigar,sequence,qualities,start,end,
    core_end,event,reference,min_baseq) {
  if (is.null(input_bam_haplotype(pos,cigar,sequence,qualities,start,end,min_baseq)))
    return("unusable")
  tokens <- regmatches(cigar,gregexpr("[0-9]+[MIDNSHP=X]",cigar))[[1L]]
  rp <- pos;qp <- 1L;edits <- list()
  for (token in tokens) {
    op <- substring(token,nchar(token));n <- as.integer(substr(token,1L,nchar(token)-1L))
    if (op %in% c("I","D")) {
      anchor <- rp-1L
      if (op=="D" && anchor<start && anchor+n>=event$pos)
        return("ambiguous") # A spanning deletion is not target REF evidence.
      if (anchor>=start && anchor<=end) {
        if (anchor<1L) return("ambiguous")
        ref <- event_reference_fetch(reference,event$contig,anchor,
          if(op=="D") anchor+n else anchor)
        alt <- if(op=="I") paste0(ref,substr(sequence,qp,qp+n-1L)) else substr(ref,1L,1L)
        edit <- tryCatch(event_normalize(reference,event$contig,anchor,ref,alt),error=function(e)NULL)
        if (is.null(edit)) return("ambiguous")
        original_end <- anchor + if(op=="D") n else 0L
        relevant <- (anchor<=core_end && original_end>=event$pos) ||
          (edit$pos<=core_end && edit$end>=event$pos)
        if (relevant) {
          observed_anchor <- input_bam_haplotype(pos,cigar,sequence,qualities,anchor,anchor,min_baseq)
          if (is.null(observed_anchor) || substr(observed_anchor,1L,1L)!=substr(ref,1L,1L))
            return("ambiguous")
          edits[[length(edits)+1L]] <- edit
        }
      }
    }
    if(op %in% c("M","=","X","D","N")) rp <- rp+n
    if(op %in% c("M","=","X","I","S")) qp <- qp+n
  }
  if (!length(edits)) return("ref")
  if (length(edits)!=1L) return("ambiguous")
  edit <- edits[[1L]]
  if (edit$pos==event$pos && edit$ref==event$ref && edit$alt==event$alt) "alt" else "other"
}

input_bam_observations <- function(source,events,project) {
  for (pkg in c("Rsamtools","GenomicRanges","IRanges"))
    if (!requireNamespace(pkg,quietly=TRUE)) stop("BAM evidence requires optional package ",pkg)
  if (!identical(source$reference_build,project$config$scientific$genome_assembly))
    stop("BAM reference_build must match selected reference; do not relabel an incompatible BAM")
  options <- input_bam_options(source$options)
  header <- Rsamtools::scanBamHeader(source$resolved_path)[[1L]]
  groups <- header$text[names(header$text)=="@RG"]
  rg <- lapply(groups,function(g) {
    x <- strsplit(g,":",fixed=TRUE); stats::setNames(vapply(x,function(a)paste(a[-1L],collapse=":"),""),vapply(x,`[[`,"",1L))
  })
  expected <- if (is.null(options$sample_name)) source$sample_id else options$sample_name
  if (!length(rg) || any(vapply(rg,function(x)is.na(x["SM"]) || x["SM"]!=expected || is.na(x["ID"]),TRUE)))
    stop("BAM read-group SM conflicts with sample binding or is missing; declare sample_name only for a verified identity mapping")
  ids <- vapply(rg,function(x)unname(x["ID"]),"")
  if (anyDuplicated(ids)) stop("Duplicate BAM read-group IDs")
  bundle <- reference_bundle(project);reference <- event_reference_open(bundle)
  on.exit(event_reference_close(reference))
  target <- header$targets
  # Non-target ALT/decoy/unplaced sequences do not invalidate a queried locus.
  # Validate every actually queried contig exactly, not unrelated header targets.
  flag <- Rsamtools::scanBamFlag(isUnmappedQuery=FALSE,isSecondaryAlignment=FALSE,
    isSupplementaryAlignment=FALSE,isNotPassingQualityControls=FALSE,
    isDuplicate=if (options$include_duplicates) NA else FALSE)
  out <- list()
  for (event in events) {
    if (event$sample_id!=source$sample_id || event$event_role!="somatic") next
    if (source$kind=="dna_bam" && options$mode=="fallback") {
      available <- input_vcf_assay(event)
      if (!is.null(available$vaf)) next
    }
    assay <- if (source$kind=="rna_bam") "rna" else if (source$kind=="normal_dna_bam") "normal" else "tumor"
    if (is.null(event$ref) || is.null(event$alt) || !event$event_type %in% c("SNV","INS","DEL") ||
        !grepl("^[ACGT]+$",event$alt)) {
      counts <- list(ref_reads=NULL,alt_reads=NULL,depth=NULL,usable_depth=NULL,
        vaf=NULL,status="unsupported",vaf_status="unavailable")
    } else {
      event_validate_ref(reference,event$contig,event$pos,event$ref)
      if (!event$contig %in% names(target)) stop("BAM lacks required variant contig: ",event$contig)
      ri <- match(event$contig,reference$fai$contig)
      if (is.na(ri) || as.numeric(target[[event$contig]])!=reference$fai$length[[ri]])
        stop("Queried BAM contig length incompatible with selected reference: ",event$contig)
      flank <- if (event$event_type=="SNV") 0L else 10L
      start <- max(1L,event$pos-flank)
      right <- event$pos+nchar(event$ref)-1L
      if (event$event_type %in% c("INS","DEL")) {
        motif <- substring(if(event$event_type=="INS") event$alt else event$ref,2L)
        last <- min(target[[event$contig]],event$pos+2000L)
        context <- if(last>event$pos) event_reference_fetch(reference,event$contig,event$pos+1L,last) else ""
        chars <- strsplit(context,"",fixed=TRUE)[[1L]]
        periodic <- rep(strsplit(motif,"",fixed=TRUE)[[1L]],length.out=length(chars))
        different <- which(chars!=periodic)
        repeat_length <- if(length(different)) different[[1L]]-1L else length(chars)
        if (!length(different) && last<target[[event$contig]])
          stop("Indel repeat exceeds 2000-base supported haplotype context; provide independently validated allele evidence")
        right <- max(right,event$pos+repeat_length+1L)
      }
      end <- min(target[[event$contig]],right+flank)
      ref <- event_reference_fetch(reference,event$contig,start,end)
      alt <- paste0(if(event$pos>start) substr(ref,1L,event$pos-start) else "",
        event$alt,substring(ref,event$pos-start+nchar(event$ref)+1L))
      region <- GenomicRanges::GRanges(event$contig,IRanges::IRanges(start,end))
      param <- Rsamtools::ScanBamParam(which=region,flag=flag,mapqFilter=options$min_mapq,
        what=c("qname","pos","cigar","seq","qual"),tag="RG")
      bam <- Rsamtools::scanBam(source$resolved_path,index=source$index_path,param=param)[[1L]]
      if (length(bam$qname) && (is.null(bam$tag$RG) || anyNA(bam$tag$RG) || any(!bam$tag$RG %in% ids)))
        stop("BAM alignment has missing/unknown read group")
      classifications <- lapply(seq_along(bam$qname),function(i) {
        if (event$event_type %in% c("INS","DEL")) return(input_bam_indel_call(
          bam$pos[[i]],bam$cigar[[i]],as.character(bam$seq[[i]]),as.character(bam$qual[[i]]),
          start,end,right,event,reference,options$min_baseq))
        h <- input_bam_haplotype(bam$pos[[i]],bam$cigar[[i]],as.character(bam$seq[[i]]),
          as.character(bam$qual[[i]]),start,end,options$min_baseq)
        if (is.null(h)) "unusable" else if (h==ref) "ref" else if (h==alt) "alt" else "other"
      })
      keys <- paste(bam$tag$RG,bam$qname,sep="\037")
      fragments <- split(as.character(unlist(classifications)),keys)
      calls <- vapply(fragments,function(x) {
        x <- unique(x[x!="unusable"])
        if (!length(x)) "unusable" else if (length(x)>1L) "ambiguous" else x[[1L]]
      },"")
      nref <- sum(calls=="ref");nalt <- sum(calls=="alt");den <- nref+nalt
      unresolved <- event$event_type %in% c("INS","DEL") && any(calls=="ambiguous")
      counts <- list(ref_reads=if(unresolved) NULL else as.integer(nref),
        alt_reads=if(unresolved) NULL else as.integer(nalt),
        resolved_ref_fragments=as.integer(nref),resolved_alt_fragments=as.integer(nalt),
        depth=as.integer(length(calls)),usable_depth=as.integer(den),
        other_reads=as.integer(sum(calls=="other")),ambiguous_reads=as.integer(sum(calls=="ambiguous")),
        vaf_denominator=as.integer(den),vaf=if (den>0 && !unresolved) nalt/den else NULL,
        reported_vaf=NULL,count_unit="fragment",vaf_status=if(den>0 && !unresolved) "derived" else "unavailable",
        status=if(unresolved) "ambiguous" else if (den==0) "insufficient_depth" else if(nalt>0) "observed_support" else "observed_no_alt")
    }
    counts$provenance <- c(somatic_source_provenance(source),list(
      method="indexed-CIGAR-normalized-target-edit-fragments-v2",options=options,
      alignment_reader="Rsamtools",alignment_reader_version=as.character(utils::packageVersion("Rsamtools")),
      read_group_ids=as.list(ids),index_sha256=source$index_sha256,
      overlap_policy="one_RG_qname_fragment; discordant mates ambiguous",
      indel_policy="normalized_target_edit; unrelated substitutions ignored; unresolved local edits withhold assay counts/VAF; resolved fragment counts retained for QC",
      reference_check="declared_assembly+queried_contig_lengths+VCF_REF; non-target header contigs ignored"))
    out[[event$event_id]] <- stats::setNames(list(counts),assay)
  }
  out
}
