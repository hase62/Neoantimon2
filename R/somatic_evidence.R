somatic_coalesce_field <- function(event,row,canonical,aliases=canonical,integer=TRUE) {
  old <- somatic_row_value(event,canonical)
  supplied <- Filter(Negate(is.null),lapply(aliases,function(name) somatic_row_value(row,name)))
  if (length(supplied)>1L && length(unique(unlist(supplied)))>1L)
    stop("Conflicting evidence aliases: ",canonical)
  value <- if (length(supplied)) supplied[[1L]] else NULL
  parse <- if (integer) somatic_integer else somatic_number
  old <- parse(old,canonical); value <- parse(value,canonical)
  if (!is.null(old) && !is.null(value) && !isTRUE(all.equal(old,value)))
    stop("Conflicting event and external evidence: ",canonical)
  if (is.null(value)) old else value
}

somatic_assay <- function(event,row,assay,source) {
  prefix <- paste0(assay,"_")
  get <- function(name,aliases=name) somatic_coalesce_field(event,row,paste0(prefix,name),
    paste0(prefix,aliases))
  ref <- get("ref_count",c("ref_count","ref_reads"))
  alt <- get("alt_count",c("alt_count","alt_reads"))
  depth <- get("raw_depth",c("raw_depth","depth"))
  usable <- get("usable_depth")
  other <- get("other_count")
  ambiguous <- get("ambiguous_count")
  declared_denominator <- get("vaf_denominator")
  denominator <- if (!is.null(declared_denominator)) declared_denominator else
    if (!is.null(ref) && !is.null(alt)) ref+alt else NULL
  if (!is.null(denominator) && denominator>.Machine$integer.max)
    stop("VAF denominator exceeds integer range")
  bad <- (!is.null(depth) && !is.null(alt) && alt>depth) ||
    (!is.null(depth) && !is.null(ref) && !is.null(alt) && ref+alt>depth) ||
    (!is.null(usable) && !is.null(depth) && usable>depth) ||
    (!is.null(declared_denominator) && !is.null(alt) && alt>declared_denominator)
  vaf <- if (!bad && !is.null(denominator) && denominator>0L && !is.null(alt))
    alt/denominator else NULL
  reported_vaf <- somatic_row_value(row,paste0(prefix,"vaf"))
  if (is.null(reported_vaf)) reported_vaf <- somatic_row_value(event,paste0(prefix,"vaf"))
  reported_vaf <- somatic_number(reported_vaf,paste0(prefix,"vaf"),0,1)
  if (!is.null(reported_vaf) && !is.null(vaf) &&
      !isTRUE(all.equal(reported_vaf,vaf,tolerance=1e-8)))
    stop("Conflicting event VAF and derived VAF")
  if (!is.null(reported_vaf) && (bad || (!is.null(depth) && depth==0L)))
    stop("Supplied VAF conflicts with invalid or zero depth")
  vaf_status <- if (!is.null(vaf)) "derived" else if (!is.null(reported_vaf))
    "supplied" else "unavailable"
  if (is.null(vaf)) vaf <- reported_vaf
  status <- if (bad) "invalid_counts" else if (!is.null(depth) && depth==0L)
    "insufficient_depth" else if (is.null(alt)) {
      if (is.null(source) && is.null(depth) && is.null(usable) && is.null(ref))
        "unavailable" else "not_evaluated"
    } else if (!is.null(denominator) && denominator==0L) "insufficient_depth" else
      if (alt>0L) "observed_support" else if (!is.null(ref) ||
        (!is.null(depth) && depth>0L)) "observed_no_alt" else "not_evaluated"
  unit <- somatic_row_value(event,paste0(prefix,"count_unit"))
  if (is.null(unit) && !is.null(row)) unit <- somatic_row_value(row,paste0(prefix,"count_unit"))
  if (is.null(unit) && !is.null(source)) unit <- "read"
  event_fields <- paste0(prefix,c("ref_count","alt_count","raw_depth",
    "usable_depth","vaf_denominator","vaf"))
  from_event <- any(vapply(event[event_fields],function(x) !is.null(x),TRUE))
  external <- somatic_source_provenance(source)
  provenance <- if (!is.null(row) && from_event) c(list(source="frozen_event+external",
    event_id=event$event_id,external=external),external[c("artifact","sha256",
      "source_version","sample_id")]) else if (!is.null(row)) external else
    if (from_event) list(source="frozen_event",event_id=event$event_id) else external
  list(ref_reads=ref,alt_reads=alt,depth=depth,usable_depth=usable,
    other_reads=other,ambiguous_reads=ambiguous,vaf_denominator=denominator,
    vaf=vaf,reported_vaf=reported_vaf,vaf_status=vaf_status,count_unit=unit,status=status,
    provenance=provenance)
}

somatic_purity <- function(event,source,facets_source=NULL) {
  rows <- if (is.null(source)) list() else source$rows
  if (length(rows)>1L) stop("Multiple purity rows for sample")
  row <- if (length(rows)) rows[[1L]] else NULL
  if (is.null(row) && !is.null(facets_source)) {
    values <- unique(vapply(facets_source$rows,function(x)
      if (is.null(x$purity)) NA_character_ else x$purity,""))
    values <- values[!is.na(values)]
    if (length(values)>1L) stop("Conflicting FACETS purity metadata")
    ploidies <- unique(vapply(facets_source$rows,function(x)
      if (is.null(x$ploidy)) NA_character_ else x$ploidy,""))
    ploidies <- ploidies[!is.na(ploidies)]
    if (length(ploidies)>1L) stop("Conflicting FACETS ploidy metadata")
    if (length(values) || length(ploidies)) {
      row <- list(sample_id=event$sample_id,purity=if (length(values)) values[[1L]] else NULL,
        ploidy=if (length(ploidies)) ploidies[[1L]] else NULL,method="FACETS")
      source <- facets_source
    }
  }
  if (!is.null(row$sample_id) && row$sample_id!=event$sample_id)
    stop("Purity sample mismatch")
  value <- somatic_coalesce_field(event,row,"purity","purity",integer=FALSE)
  if (!is.null(value) && (value<0 || value>1)) stop("Invalid purity fraction")
  ploidy <- somatic_coalesce_field(event,row,"ploidy","ploidy",integer=FALSE)
  list(value=value,ploidy=ploidy,status=if (is.null(value)) "unavailable" else "available",
    method=if (!is.null(row)) somatic_row_value(row,"method") else NULL,
    provenance=if (is.null(source)) list(source="frozen_event",event_id=event$event_id) else
      somatic_source_provenance(source))
}

somatic_contig <- function(value) {
  if (is.null(value)) return(NULL)
  x <- toupper(sub("^CHR","",as.character(value),ignore.case=TRUE))
  if (x=="23") x <- "X"
  if (x=="24") x <- "Y"
  if (x %in% c("25","M")) x <- "MT"
  x
}

somatic_facets_segments <- function(source) {
  if (is.null(source)) return(list())
  required <- c("chrom","loc.start","loc.end","tcn.em","lcn.em")
  if (length(source$rows) && !all(required %in% names(source$rows[[1L]])))
    stop("FACETS table missing required cncf columns")
  lapply(seq_along(source$rows),function(i) {
    r <- source$rows[[i]]
    declared_build <- if (!is.null(r$reference_build)) r$reference_build else r$genome_assembly
    if (!is.null(declared_build) && declared_build!=source$reference_build)
      stop("FACETS row/reference build mismatch")
    if (!is.null(r$sample_id) && r$sample_id!=source$sample_id)
      stop("FACETS row/sample mismatch")
    start <- somatic_integer(r[["loc.start"]],"FACETS loc.start")
    end <- somatic_integer(r[["loc.end"]],"FACETS loc.end")
    if (is.null(start) || is.null(end) || start<1L || end<start)
      stop("Invalid FACETS segment coordinates")
    total <- somatic_number(r[["tcn.em"]],"FACETS tcn.em",minimum=0)
    minor <- somatic_number(r[["lcn.em"]],"FACETS lcn.em",minimum=0)
    if (!is.null(minor) && !is.null(total) && minor>total)
      stop("FACETS minor CN exceeds total CN")
    normal <- somatic_number(r[["normal_cn"]],"FACETS normal_cn",minimum=0)
    list(segment_id=paste0("cn-segment:",infra_hash(list(source=source$source_id,row=i,
      chrom=r$chrom,start=start,end=end,total=total,minor=minor))),
      sample_id=source$sample_id,contig=somatic_contig(r$chrom),raw_contig=r$chrom,
      start=start,end=end,total_cn=total,minor_cn=minor,
      major_cn=if (!is.null(total) && !is.null(minor)) total-minor else NULL,
      normal_cn=normal,cellular_fraction=somatic_number(r[["cf.em"]],"FACETS cf.em",0,1),
      num_mark=somatic_integer(r[["num.mark"]],"FACETS num.mark"),
      nhet=somatic_integer(r[["nhet"]],"FACETS nhet"),
      quality=somatic_number(r[["quality"]],"FACETS quality"),
      purity=somatic_number(r[["purity"]],"FACETS purity",0,1),
      ploidy=somatic_number(r[["ploidy"]],"FACETS ploidy",minimum=0),
      reference_build=source$reference_build,provenance=somatic_source_provenance(source))
  })
}

somatic_cn_map <- function(event,segments,source) {
  empty <- function(status) list(status=status,segment_id=NULL,total_cn=NULL,minor_cn=NULL,
    major_cn=NULL,normal_cn=NULL,segment=NULL,provenance=somatic_source_provenance(source))
  if (is.null(source)) {
    total <- somatic_number(event$total_cn,"event total_cn",minimum=0)
    minor <- somatic_number(event$minor_cn,"event minor_cn",minimum=0)
    if (is.null(total)) return(empty("unavailable"))
    return(list(status="provided_event",segment_id=NULL,total_cn=total,
      minor_cn=minor,major_cn=somatic_number(event$major_cn,"event major_cn",minimum=0),
      normal_cn=somatic_number(event$normal_cn,"event normal_cn",minimum=0),
      segment=NULL,provenance=list(source="frozen_event",event_id=event$event_id,
        caller=event$cn_source,method=event$cn_method)))
  }
  if (is.null(event$genome_assembly) || is.null(source$reference_build))
    return(empty("build_unknown"))
  if (event$genome_assembly!=source$reference_build) return(empty("build_incompatible"))
  if (is.null(event$contig) || is.null(event$pos)) return(empty("event_coordinate_unavailable"))
  same <- Filter(function(x) identical(x$contig,somatic_contig(event$contig)),segments)
  if (!length(same)) return(empty(if (length(segments)) "contig_mismatch" else "unavailable"))
  event_end <- if (is.null(event$end)) event$pos else event$end
  hits <- Filter(function(x) x$start<=event_end && x$end>=event$pos,same)
  if (!length(hits)) return(empty("no_overlap"))
  if (length(hits)>1L) return(empty("ambiguous_overlap"))
  hit <- hits[[1L]]
  if (hit$start>event$pos || hit$end<event_end) return(empty("partial_overlap"))
  list(status="mapped",segment_id=hit$segment_id,total_cn=hit$total_cn,
    minor_cn=hit$minor_cn,major_cn=hit$major_cn,normal_cn=hit$normal_cn,
    segment=hit,provenance=hit$provenance)
}

somatic_ccf <- function(event,row,source,dna,purity,cn) {
  external_raw <- if (!is.null(row)) somatic_row_value(row,"ccf_raw") else NULL
  external <- if (!is.null(row)) somatic_row_value(row,"ccf") else NULL
  if (is.null(external_raw)) external_raw <- somatic_row_value(event,"ccf_raw")
  if (is.null(external)) external <- somatic_row_value(event,"ccf")
  if (!is.null(external_raw) || !is.null(external)) {
    raw <- somatic_number(if (is.null(external_raw)) external else external_raw,"supplied CCF")
    bounded <- if (raw>=0 && raw<=1) raw else NULL
    if (!is.null(external)) {
      supplied_bounded <- somatic_number(external,"supplied bounded CCF",0,1)
      if (!is.null(bounded) && !isTRUE(all.equal(bounded,supplied_bounded)))
        stop("Conflicting supplied CCF values")
      bounded <- supplied_bounded
    }
    return(list(status=if (raw<0 || raw>1) "out_of_range" else "external",
      identifiability="external",raw_estimate=raw,bounded_value=bounded,
      lower=somatic_number(if (is.null(row)) event$ccf_lower else row$ccf_lower,
        "CCF lower",0,1),
      upper=somatic_number(if (is.null(row)) event$ccf_upper else row$ccf_upper,
        "CCF upper",0,1),
      solutions=list(),multiplicity=NULL,model_version=NULL,
      method=if (is.null(row)) event$ccf_method else row$method,
      reported_source=if (is.null(row)) event$ccf_source else row$source,
      provenance=if (is.null(source)) list(source="frozen_event",event_id=event$event_id) else
        somatic_source_provenance(source)))
  }
  unavailable <- function(reason) list(status="not_evaluable",reason=reason,
    identifiability="not_evaluable",raw_estimate=NULL,bounded_value=NULL,
    lower=NULL,upper=NULL,solutions=list(),multiplicity=NULL,
    model_version=somatic_ccf_model_version(),provenance=NULL)
  if (!identical(event$event_role,"somatic"))
    return(unavailable("non_somatic_event"))
  if (is.null(dna$vaf)) return(unavailable("tumor_vaf_unavailable"))
  if (is.null(purity$value) || purity$value<=0) return(unavailable("purity_unavailable"))
  if (!cn$status %in% c("mapped","provided_event") || is.null(cn$total_cn) || cn$total_cn<=0)
    return(unavailable("local_total_cn_unavailable"))
  normal_cn <- if (!is.null(cn$normal_cn)) cn$normal_cn else
    somatic_number(event$normal_cn,"event normal_cn",minimum=0)
  if (is.null(normal_cn)) return(unavailable("normal_cn_unavailable"))
  multiplicity <- if (!is.null(row)) somatic_row_value(row,"mutation_multiplicity") else NULL
  if (is.null(multiplicity)) multiplicity <- somatic_row_value(event,"mutation_multiplicity")
  if (!is.null(multiplicity)) {
    possible <- somatic_number(multiplicity,"mutation multiplicity",minimum=0)
    if (possible<=0 || possible>cn$total_cn) return(unavailable("invalid_multiplicity"))
  } else {
    if (cn$total_cn!=floor(cn$total_cn) || cn$total_cn>.Machine$integer.max)
      return(unavailable("multiplicity_unresolved"))
    possible <- seq_len(as.integer(cn$total_cn))
  }
  solutions <- lapply(possible,function(m) {
    raw <- dna$vaf*(purity$value*cn$total_cn+(1-purity$value)*normal_cn)/
      (purity$value*m)
    list(multiplicity=as.numeric(m),raw_estimate=raw,
      bounded_value=if (raw>=0 && raw<=1) raw else NULL,
      compatible=raw>=0 && raw<=1)
  })
  compatible <- Filter(function(x) x$compatible,solutions)
  status <- if (!length(compatible)) "out_of_range" else
    if (length(compatible)>1L) "ambiguous" else "conditional"
  selected <- if (length(compatible)==1L) compatible[[1L]] else NULL
  list(status=status,reason=NULL,
    identifiability=if (status=="ambiguous") "ambiguous" else "conditional",
    raw_estimate=if (is.null(selected)) NULL else selected$raw_estimate,
    bounded_value=if (is.null(selected)) NULL else selected$bounded_value,
    lower=NULL,upper=NULL,solutions=solutions,
    multiplicity=if (is.null(selected)) NULL else selected$multiplicity,
    inputs=list(tumor_vaf=dna$vaf,tumor_alt_reads=dna$alt_reads,
      tumor_ref_reads=dna$ref_reads,purity=purity$value,total_cn=cn$total_cn,
      normal_cn=normal_cn),model_version=somatic_ccf_model_version(),
    method="conditional_mixture_formula",
    provenance=list(source="Neoantimon2",model_version=somatic_ccf_model_version()))
}

somatic_loh <- function(event,row,source,cn) {
  if (!is.null(row)) {
    status <- somatic_row_value(row,"status")
    if (is.null(status) || !status %in% c("no_evidence","loh_supported",
      "allele_specific_loss","copy_neutral_loh","unresolved"))
      stop("Invalid supplied LOH status")
    return(list(status=status,scope="generic_locus",method="supplied",
      reported_allele=row$allele,locus=row$locus,
      provenance=somatic_source_provenance(source)))
  }
  if (!cn$status %in% c("mapped","provided_event") || is.null(cn$minor_cn))
    return(list(status="unavailable",scope="generic_locus",method=NULL,
      provenance=cn$provenance))
  status <- if (cn$minor_cn>0) "no_evidence" else
    if (!is.null(cn$normal_cn) && cn$total_cn==cn$normal_cn)
      "copy_neutral_loh" else "loh_supported"
  list(status=status,scope="generic_locus",method="minor_cn_zero",
    provenance=cn$provenance)
}

somatic_expression_records <- function(sources,origins=list()) {
  out <- list()
  for (source in Filter(function(x) x$kind=="expression",sources)) {
    for (row in source$rows) {
      tx <- somatic_row_value(row,"transcript_id")
      gene <- somatic_row_value(row,"gene_id")
      if (is.null(tx) && is.null(gene)) stop("Expression row needs transcript_id or gene_id")
      if (!is.null(row$sample_id) && row$sample_id!=source$sample_id)
        stop("Expression sample mismatch")
      gene_tpm <- somatic_number(row$gene_tpm,"gene_tpm",minimum=0)
      transcript_tpm <- somatic_number(row$transcript_tpm,"transcript_tpm",minimum=0)
      record <- list(sample_id=source$sample_id,transcript_id=tx,gene_id=gene,
        gene_tpm=gene_tpm,transcript_tpm=transcript_tpm,
        status=if (is.null(gene_tpm) && is.null(transcript_tpm)) "unavailable" else "available",
        provenance=somatic_source_provenance(source))
      record$evidence_id <- paste0("expression-evidence:",infra_hash(record))
      out[[length(out)+1L]] <- record
    }
  }
  keys <- vapply(out,function(x) paste(x$sample_id,x$transcript_id,x$gene_id,sep="\037"),"")
  if (anyDuplicated(keys)) stop("Duplicate expression transcript/gene row")
  for (origin in origins) {
    gene_tpm <- somatic_number(origin$gene_tpm,"origin gene_tpm",minimum=0)
    transcript_tpm <- somatic_number(origin$transcript_tpm,"origin transcript_tpm",minimum=0)
    if (is.null(gene_tpm) && is.null(transcript_tpm)) next
    record <- list(sample_id=origin$sample_id,transcript_id=origin$transcript_id,
      gene_id=NULL,gene_tpm=gene_tpm,transcript_tpm=transcript_tpm,
      status="available",provenance=list(source="frozen_peptide_origin",
        origin_id=origin$origin_id))
    record$evidence_id <- paste0("expression-evidence:",infra_hash(record))
    out[[length(out)+1L]] <- record
  }
  if (length(out)) out <- out[!duplicated(vapply(out,`[[`,"","evidence_id"))]
  out
}

somatic_evidence_row <- function(event,domain,metric,value,status,unit,provenance) {
  if (!is.null(value) && !(is.numeric(value) || is.character(value) || is.logical(value)))
    stop("Unsupported somatic evidence value type")
  row <- list(sample_id=event$sample_id,analysis_context=event$analysis_context,
    target_type="event",target_id=event$event_id,evidence_domain=domain,
    metric=metric,value_number=if (is.numeric(value)) as.numeric(value) else NULL,
    value_text=if (is.character(value)) value else if (is.logical(value)) as.character(value) else NULL,
    unit=unit,status=status,source=if (is.null(provenance)) NULL else provenance$source,
    source_version=if (is.null(provenance)) NULL else provenance$source_version,
    artifact=if (is.null(provenance)) NULL else provenance$artifact,
    checksum=if (is.null(provenance)) NULL else provenance$sha256,
    provenance=provenance,schema_version=somatic_evidence_version())
  row$evidence_id <- paste0("somatic-evidence:",infra_hash(row))
  row
}

somatic_event_metric_rows <- function(event,record) {
  out <- list()
  add <- function(domain,metric,value,status,unit,provenance) {
    out[[length(out)+1L]] <<- somatic_evidence_row(event,domain,metric,value,status,unit,provenance)
  }
  for (assay in c("tumor","normal","rna")) {
    x <- record[[assay]]
    for (field in c("ref_reads","alt_reads","depth","usable_depth","vaf_denominator","vaf"))
      add(if (assay=="rna") "rna" else "dna",paste0(assay,"_",field),x[[field]],
        x$status,if (field=="vaf") "fraction" else x$count_unit,x$provenance)
  }
  add("purity","purity",record$purity$value,record$purity$status,"fraction",
    record$purity$provenance)
  for (field in c("total_cn","minor_cn","major_cn","normal_cn"))
    add("copy_number",field,record$cn[[field]],record$cn$status,"copies",record$cn$provenance)
  add("ccf","ccf",record$ccf$bounded_value,record$ccf$status,"fraction",
    record$ccf$provenance)
  add("loh","generic_loh",record$loh$status,record$loh$status,"category",
    record$loh$provenance)
  out
}

somatic_validate_source_targets <- function(sources,events) {
  by_id <- stats::setNames(events,vapply(events,`[[`,"","event_id"))
  for (source in sources) if (source$kind %in% c("dna_reads","rna_reads","ccf","loh")) {
    ids <- vapply(source$rows,function(row) {
      id <- somatic_row_value(row,"event_id")
      if (is.null(id)) stop("Event evidence row requires event_id")
      id
    },"")
    if (anyDuplicated(ids)) stop("Duplicate evidence event_id")
    for (id in ids) {
      event <- by_id[[id]]
      if (is.null(event) || event$sample_id!=source$sample_id)
        stop("Evidence source event/sample identity mismatch")
    }
  }
  invisible(TRUE)
}

somatic_event_record <- function(event,sources,segments_by_sample) {
  source <- function(kind) somatic_source_index(sources,kind,event$sample_id)
  dna_source <- source("dna_reads"); rna_source <- source("rna_reads")
  cn_source <- source("facets_cn"); ccf_source <- source("ccf")
  loh_source <- source("loh")
  tumor <- somatic_assay(event,somatic_match_event_row(dna_source,event),"tumor",dna_source)
  normal <- somatic_assay(event,somatic_match_event_row(dna_source,event),"normal",dna_source)
  rna <- somatic_assay(event,somatic_match_event_row(rna_source,event),"rna",rna_source)
  purity <- somatic_purity(event,source("purity"),cn_source)
  cn <- somatic_cn_map(event,segments_by_sample[[event$sample_id]],cn_source)
  ccf <- somatic_ccf(event,somatic_match_event_row(ccf_source,event),ccf_source,
    tumor,purity,cn)
  loh <- somatic_loh(event,somatic_match_event_row(loh_source,event),loh_source,cn)
  record <- list(event_id=event$event_id,event_key=event$event_key,
    sample_id=event$sample_id,analysis_context=event$analysis_context,
    source_event_role=event$event_role,tumor=tumor,normal=normal,rna=rna,
    purity=purity,cn=cn,ccf=ccf,loh=loh,
    source_expression_json=event$expression_json,
    source_cn_by_breakpoint_json=event$cn_by_breakpoint_json,
    frozen_ccf_annotations=list(status=event$ccf_status,
      identifiability=event$ccf_identifiability,
      candidates=event$ccf_candidates_json,
      lower=event$ccf_lower,upper=event$ccf_upper),
    schema_version=somatic_evidence_version())
  record$evidence_id <- paste0("somatic-event-evidence:",infra_hash(record))
  record
}

somatic_candidate_links <- function(candidates,origins,peptides,event_evidence,expression,sources) {
  origin_index <- stats::setNames(origins,vapply(origins,`[[`,"","origin_id"))
  peptide_index <- stats::setNames(peptides,vapply(peptides,`[[`,"","peptide_id"))
  event_index <- stats::setNames(event_evidence,vapply(event_evidence,`[[`,"","event_id"))
  lapply(candidates,function(candidate) {
    paths <- lapply(unlist(candidate$origin_ids),function(id) {
      origin <- origin_index[[id]]
      if (is.null(origin)) stop("Candidate origin missing from frozen peptide artifact")
      peptide <- peptide_index[[id]]
      if (is.null(peptide)) stop("Origin peptide provenance missing")
      somatic_ids <- sort(unique(as.character(unlist(origin$source_event_ids))),method="radix")
      germline_ids <- sort(unique(as.character(unlist(origin$germline_context_ids))),method="radix")
      all_ids <- sort(unique(c(somatic_ids,germline_ids)),method="radix")
      if (any(!all_ids %in% names(event_index))) stop("Origin source event evidence missing")
      if (any(vapply(event_index[somatic_ids],function(x)x$source_event_role!="somatic",TRUE)) ||
          any(vapply(event_index[germline_ids],function(x)x$source_event_role!="germline",TRUE)))
        stop("Origin source event role mismatch")
      tx <- Filter(function(x) x$sample_id==candidate$sample_id &&
        !is.null(x$transcript_id) &&
        sub("\\.[0-9]+$","",x$transcript_id)==origin$transcript_id,expression)
      has_expression_source <- any(vapply(sources,function(x)
        x$kind=="expression" && x$sample_id==candidate$sample_id,TRUE))
      list(origin_id=id,peptide_id=origin$sample_peptide_id,
        reconstruction_ids=as.list(sort(unique(vapply(peptide$provenance_occurrences,
          `[[`,"","reconstruction_id")),method="radix")),
        translation_id=origin$sequence_context_id,
        somatic_event_ids=as.list(somatic_ids),germline_event_ids=as.list(germline_ids),
        event_evidence_ids=as.list(vapply(event_index[all_ids],`[[`,"","evidence_id")),
        expression_evidence_ids=as.list(vapply(tx,`[[`,"","evidence_id")),
        expression_status=if (any(vapply(tx,function(x)x$status=="available",TRUE)))
          "available" else if (length(tx)) "unavailable" else
          if (has_expression_source) "not_evaluated" else "unavailable")
    })
    list(candidate_id=candidate$candidate_id,sample_id=candidate$sample_id,
      origin_paths=paths)
  })
}

somatic_raw_construct <- function(events,peptide,sources,run_id=NULL) {
  somatic_validate_source_targets(sources,events)
  segments <- lapply(Filter(function(x) x$kind=="facets_cn",sources),somatic_facets_segments)
  by_sample <- stats::setNames(segments,vapply(Filter(function(x) x$kind=="facets_cn",sources),
    `[[`,"","sample_id"))
  expression <- somatic_expression_records(sources,peptide$origins)
  event_evidence <- lapply(events,somatic_event_record,sources=sources,segments_by_sample=by_sample)
  rows <- unlist(Map(somatic_event_metric_rows,events,event_evidence),recursive=FALSE)
  for (x in expression) for (metric in c("gene_tpm","transcript_tpm")) {
    row <- list(sample_id=x$sample_id,analysis_context=NULL,
      target_type=if (is.null(x$transcript_id)) "gene" else "transcript",
      target_id=if (is.null(x$transcript_id)) x$gene_id else x$transcript_id,
      evidence_domain="expression",metric=metric,
      value_number=x[[metric]],value_text=NULL,unit="TPM",status=x$status,
      source=x$provenance$source,source_version=x$provenance$source_version,
      artifact=x$provenance$artifact,checksum=x$provenance$sha256,
      provenance=x$provenance,schema_version=somatic_evidence_version(),
      expression_evidence_id=x$evidence_id)
    row$evidence_id <- paste0("somatic-evidence:",infra_hash(row))
    rows[[length(rows)+1L]] <- row
  }
  if (!is.null(run_id)) {
    for (i in seq_along(rows)) rows[[i]]$run_id <- run_id
    for (i in seq_along(event_evidence)) event_evidence[[i]]$run_id <- run_id
  }
  list(schema_version=somatic_evidence_version(),run_id=run_id,event_evidence=event_evidence,
    evidence_rows=rows,expression_evidence=expression,
    cn_segments=unlist(segments,recursive=FALSE),
    candidate_links=somatic_candidate_links(peptide$candidates,peptide$origins,
      peptide$peptides,
      event_evidence,expression,sources),
    source_provenance=lapply(sources,somatic_source_provenance))
}
