somatic_test_event <- function(id="event:a",pos=120L,contig="chr1",build="GRCh38")
  list(event_id=id,event_key=paste0("key:",id),sample_id="sample-1",
    analysis_context="somatic:patient-1",event_role="somatic",
    genome_assembly=build,contig=contig,pos=pos,end=pos)

somatic_test_peptide <- function(compound=FALSE) {
  origins <- list(list(origin_id="origin:a",sample_peptide_id="peptide:a",
    sample_id="sample-1",transcript_id="TX1",source_event_ids=list("event:a"),
    germline_context_ids=list()),
    list(origin_id="origin:b",sample_peptide_id="peptide:b",
      sample_id="sample-1",transcript_id="TX2",
      source_event_ids=if (compound) list("event:a","event:b") else list("event:a"),
      germline_context_ids=list()))
  peptides <- lapply(origins,function(x) list(peptide_id=x$origin_id,
    provenance_occurrences=list(list(reconstruction_id=paste0("reconstruction:",x$origin_id)))))
  candidates <- lapply(seq_along(origins),function(i) list(
    candidate_id=paste0("candidate:",i),sample_id="sample-1",
    origin_ids=list(origins[[i]]$origin_id)))
  list(origins=origins,peptides=peptides,candidates=candidates)
}

somatic_test_source <- function(kind,rows=list(),build=NULL) {
  x <- list(kind=kind,sample_id="sample-1",path="fixture.tsv",source=kind,
    source_version="fixture-1",reference_build=build,
    resolved_path="fixture.tsv",sha256=paste(rep("a",64),collapse=""),rows=rows)
  x$source_id <- paste0("source:",kind)
  x
}

somatic_test_row <- function(event_id="event:a",tumor_alt="5",normal_alt="1",
                             tumor_ref="5",normal_ref="19",depth="10")
  list(event_id=event_id,tumor_ref_reads=tumor_ref,tumor_alt_reads=tumor_alt,
    tumor_depth=depth,normal_ref_reads=normal_ref,normal_alt_reads=normal_alt,
    normal_depth="20")
