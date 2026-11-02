# Read-only measurement projection. Never used for grading or native ranking.
input_display_candidates <- function(rows, origins, event_evidence) {
  oi <- external_group_index(vapply(origins, `[[`, "", "origin_id"))
  ei <- external_group_index(vapply(event_evidence, `[[`, "", "event_id"))
  lapply(rows, function(row) {
    paths <- lapply(unlist(row$origin_ids),function(id) {
      at <- external_index_get(oi,id)
      if(length(at)!=1L) return(list(origin_id=id,status="origin_unavailable",events=list(),vaf=NULL))
      origin <- origins[[at]]; ids <- as.character(unlist(origin$source_event_ids))
      obs <- lapply(ids,function(eid) {
        at <- external_index_get(ei,eid)
        e <- if(length(at)==1L) event_evidence[[at]] else NULL
        if(!is.null(e) && !identical(e$sample_id,row$sample_id)) stop("Display evidence sample mismatch")
        list(event_id=eid,evidence_id=e$evidence_id,observation=e$tumor)
      })
      value <- if(length(obs)==1L) obs[[1L]]$observation$vaf else NULL
      list(origin_id=id,source_event_ids=as.list(ids),events=obs,vaf=value,
        status=if(length(ids)>1L) "multiple_events" else if(is.null(value)) "unavailable" else "observed")
    })
    selected <- row$priority_origin_id
    at <- if(is.null(selected)) integer() else which(vapply(paths,function(x)identical(x$origin_id,selected),TRUE))
    chosen <- if(length(at)==1L) paths[at] else paths
    values <- lapply(chosen,`[[`,"vaf")
    valid <- length(values)>0L && all(vapply(values,function(x)
      is.numeric(x) && length(x)==1L && is.finite(x) && x>=0 && x<=1,TRUE))
    same <- valid && length(unique(unlist(values)))==1L
    # No max/mean or missing-as-zero. A frozen selected Origin may define the
    # scalar; otherwise every linked Origin must have a single agreeing value.
    row["tumor_vaf"] <- list(if(same) values[[1L]] else NULL)
    row$dna_vaf_display <- list(status=if(same) "observed" else if(valid) "conflicting_origins" else "unavailable_or_ambiguous",
      policy=if(length(at)==1L) "frozen_selected_origin" else "all_origins_require_agreement",
      selected_origin_id=selected,origin_observations=paths)
    row
  })
}

input_display_rows <- function(project, rows, origins) {
  raw <- input_evidence_read(project)
  if(is.null(raw)) return(rows) # Preserve compatible older persisted projects.
  input_display_candidates(rows,origins,raw$event_evidence)
}
