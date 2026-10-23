#!/usr/bin/env Rscript
# Internal component harness. No public API and no reconstruction/predictors.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
repo <- normalizePath(args[[1L]], mustWork = TRUE)
request <- jsonlite::read_json(args[[2L]], simplifyVector = TRUE)
outfile <- args[[3L]]
pkgload::load_all(repo, quiet = TRUE, helpers = FALSE, attach_testthat = FALSE,
  export_all = FALSE, export_imports = FALSE)
warnings <- character()
started <- proc.time()[['elapsed']]
result <- tryCatch(withCallingHandlers({
  config <- neo_config(env = character(), scientific = list(genome_assembly = request$assembly),
    paths = list(data_root = dirname(request$path), run_root = tempdir()))
  project <- neo_project(samples = data.frame(sample_id = request$sample_id, patient_id = request$patient_id),
    inputs = list(list(input_id = 'public-vcf', input_role = request$role,
      sample_id = request$sample_id, path = request$path)), config = config)
  fai <- read.delim(paste0(request$reference, '.fai'), header = FALSE,
    col.names = c('contig', 'length', 'offset', 'line_bases', 'line_width'))
  con <- file(request$reference, 'rb')
  tryCatch({
    reference <- list(con = con, fai = fai, assembly = request$assembly)
    parsed <- Neoantimon2:::event_parse_vcf(project, project$inputs[[1L]], reference)
    rows <- Neoantimon2:::event_merge_rows(parsed)
    source_view <- function(row, source) {
      list(line = source$line, input_id = source$input_id, sample_id = source$sample_id,
        alt_index = source$selected_alt_index, raw_info = source$raw_info,
        raw_format = source$raw_format, projected_info = source$projected_info,
        projected_format = source$projected_format,
        phase = if (is.null(source$raw_format)) list() else
          lapply(Neoantimon2:::phase_source_assertions(row, source), function(a)
            a[c('status','genotype_state','raw_gt','raw_ps','raw_pid','raw_pgt','mechanism','block_label')]))
    }
    list(success = TRUE, raw_event_count = length(parsed), unique_event_count = length(rows),
      rows = lapply(rows, function(row) list(event_id = row$event_id, event_key = row$event_key,
        patient_id = row$patient_id, sample_id = row$sample_id, contig = row$contig,
        pos = row$pos, ref = row$ref, alt = row$alt, event_type = row$event_type,
        reconstruction_status = row$reconstruction_status,
        sources = lapply(row$event_details_json$source_records, function(s) source_view(row,s)))))
  }, finally = close(con))
}, warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning') }),
error = function(e) list(success = FALSE, error = conditionMessage(e)))
result$seconds <- proc.time()[['elapsed']] - started
result$warnings <- as.list(warnings)
result$package_version <- as.character(packageVersion('Neoantimon2'))
result$exports <- as.list(sort(getNamespaceExports('Neoantimon2')))
jsonlite::write_json(result, outfile, auto_unbox = TRUE, null = 'null', pretty = FALSE, digits = NA)
