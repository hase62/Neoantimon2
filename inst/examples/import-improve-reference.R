# Explicit local importer. No automatic download; upstream records are NOT bundled.
# Source: IMPROVE_paper c67094276b52358b5202ad9e5597acd87888ea1d
# data.zip member data/01_data/01_Validated_neoepitopes.txt
# Borch et al. 2024, DOI 10.3389/fimmu.2024.1360281, Methods/T-cell screening.
import_improve_reference <- function(path, output) {
  expected <- "b27854c98a34b6bb83af248979cc3ea121534fa2e8d9e97415654dc0c0a95c95"
  checksum <- digest::digest(file=path,algo="sha256",serialize=FALSE)
  if (!identical(checksum,expected)) stop("Unsupported IMPROVE file revision; audit a new adapter before use")
  if (file.exists(output)) stop("Do not overwrite an existing reference")
  x <- utils::read.delim(path,colClasses="character",check.names=FALSE,quote="",comment.char="")
  required <- c("Patient","HLA_allele","Mut_peptide","response","validation","identity")
  stopifnot(all(required %in% names(x)),nrow(x)==17520,all(x$validation=="Sufficient"),
    all(x$response %in% c("yes","no","1","0")),all(nzchar(x$Patient)))
  # Patient and study identity are preserved, never inferred from filenames.
  records <- lapply(seq_len(nrow(x)),function(i)list(record_id=paste0("IMPROVE-source-line:",i+1L),
    source_record_id=x$identity[i],peptide_sequence=x$Mut_peptide[i],mhc_key=x$HLA_allele[i],
    endpoint="T_cell_recognition",assay="DNA_barcoded_pMHC_multimer",
    result=if(x$response[i] %in% c("yes","1")) "positive" else "negative",
    patient_id=x$Patient[i],study_id="10.3389/fimmu.2024.1360281",citation="https://doi.org/10.3389/fimmu.2024.1360281",
    source_response=x$response[i],source_validation=x$validation[i]))
  ref <- list(schema_version="experimental-reference-0.1.0",id="IMPROVE-tested-reference",
    version="c67094276b52358b5202ad9e5597acd87888ea1d:validated-table-v1",
    patient_identity="study_local",provenance=list(citation="10.3389/fimmu.2024.1360281",
      url="https://github.com/SRHgroup/IMPROVE_paper/tree/c67094276b52358b5202ad9e5597acd87888ea1d",
      archive_member="data/01_data/01_Validated_neoepitopes.txt",sha256=checksum,
      assay_source="Publication Methods: T-cell screening; DNA-barcoded pMHC multimers across patient-derived samples",
      limitations="Reference positives/negatives are assay-specific observations in their source patients; no HCC1395 truth and no BigMHC calibration. Repository redistribution terms unresolved; do not redistribute records."),records=records)
  jsonlite::write_json(ref,output,auto_unbox=TRUE,pretty=FALSE,null="null",digits=NA)
  invisible(output)
}
if (sys.nframe()==0L) {
  a <- commandArgs(TRUE);if(length(a)!=2L) stop("Usage: Rscript import-improve-reference.R original-table.txt new-reference.json")
  import_improve_reference(a[1],a[2])
}
