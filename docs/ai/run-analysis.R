# Repository-only AI workflow helper; NOT a Neoantimon2 package API or YAML importer.
# Default: read-only preflight. --execute requires the user's analysis authorization.
agent_scalar <- function(x, label) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x) || x != trimws(x))
    stop('Provide one explicit ', label)
  x
}
agent_absolute <- function(x, label, exists = TRUE) {
  x <- path.expand(agent_scalar(x, label))
  if (!grepl('^/|^[A-Za-z]:[/\\\\]', x)) stop(label, ' must be absolute')
  normalizePath(x, mustWork = exists)
}
agent_hla <- function(x) {
  if (!is.character(x) || !length(x) || anyNA(x) || anyDuplicated(x)) stop('Missing/duplicate HLA declarations')
  allele <- 'HLA-[ABC]\\*[0-9]{2}:[0-9]{2}'
  pair <- paste0('(?:HLA-DRA\\*[0-9]{2}:[0-9]{2}~HLA-DRB[1345]\\*[0-9]{2}:[0-9]{2}',
    '|HLA-DPA1\\*[0-9]{2}:[0-9]{2}~HLA-DPB1\\*[0-9]{2}:[0-9]{2}',
    '|HLA-DQA1\\*[0-9]{2}:[0-9]{2}~HLA-DQB1\\*[0-9]{2}:[0-9]{2})')
  if (any(!grepl(paste0('^(?:', allele, '|', pair, ')$'), x, perl = TRUE)))
    stop('Ambiguous HLA: supply canonical two-field Class I alleles or explicit Class II alpha~beta pairs')
  x # Syntax only: not HLA typing or proof of biological allele existence.
}
agent_vcf_header <- function(path, sample_id) {
  con <- if (grepl('\\.gz$', path)) gzfile(path, 'rt') else file(path, 'rt')
  on.exit(close(con))
  for (i in seq_len(100000L)) {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) break
    if (startsWith(line, '#CHROM\t')) {
      h <- strsplit(line, '\t', fixed = TRUE)[[1L]]
      if (length(h) != 10L) stop('Require exactly one VCF sample column; multi-sample VCF unsupported')
      if (h[[10L]] != sample_id) stop('VCF sample identity differs from explicit sample_id')
      return(invisible(TRUE))
    }
    if (!startsWith(line, '#')) stop('Missing VCF #CHROM header')
  }
  stop('Missing/oversized VCF header; no inferred sample mapping')
}
agent_declare <- function(manifest) {
  stopifnot(as.character(utils::packageVersion('Neoantimon2')) == '0.3.0')
  allowed <- c('workflow_declaration', 'mode', 'reference_assembly', 'resource_config',
    'run_directory', 'samples', 'hla', 'prediction_providers', 'optional_evidence', 'external_predictions')
  if (!is.list(manifest) || any(!names(manifest) %in% allowed) ||
      !identical(manifest$workflow_declaration, 'neoantimon2-agent-v1') || !identical(manifest$mode, 'somatic'))
    stop('Require the documented somatic AI declaration; other modes follow the shared playbook')
  assembly <- agent_scalar(manifest$reference_assembly, 'reference assembly')
  config_file <- agent_absolute(manifest$resource_config, 'resource_config')
  output <- agent_absolute(manifest$run_directory, 'new run_directory', exists = FALSE)
  if (file.exists(output)) stop('Output session already exists; retain it and choose a new directory')
  for (key in c('optional_evidence', 'external_predictions'))
    if (length(manifest[[key]])) stop(key, ' requires the explicit documented evidence/external workflow; helper will not ignore it')
  samples <- manifest$samples
  if (!is.list(samples) || !length(samples) || !is.null(names(samples))) stop('samples must be a record array')
  for (s in samples) {
    if (!setequal(names(s), c('sample_id', 'patient_id', 'vcf'))) stop('Require sample_id, patient_id, vcf per sample')
    agent_scalar(s$sample_id, 'sample_id'); agent_scalar(s$patient_id, 'patient_id')
  }
  ids <- vapply(samples, `[[`, '', 'sample_id'); patients <- vapply(samples, `[[`, '', 'patient_id')
  if (anyDuplicated(ids)) stop('Duplicate sample_id')
  if (!is.list(manifest$hla) || anyDuplicated(names(manifest$hla)) || !setequal(names(manifest$hla), unique(patients)))
    stop('Missing/unknown patient-level HLA association')
  hla <- lapply(manifest$hla, agent_hla)
  paths <- vapply(samples, function(s) agent_absolute(s$vcf, 'VCF path'), '')
  for (i in seq_along(paths)) agent_vcf_header(paths[[i]], ids[[i]])
  cfg <- Neoantimon2::neo_config(file = config_file, env = character(),
    paths = list(run_root = file.path(output, 'runs')))
  if (!identical(cfg$scientific$genome_assembly, assembly)) stop('Manifest/config reference assembly mismatch')
  requested <- unlist(manifest$prediction_providers, use.names = FALSE)
  configured <- unique(vapply(cfg$predictors, `[[`, '', 'provider'))
  if (!setequal(requested, configured)) stop('Requested providers and resource configuration differ; review explicit predictor settings')
  # Check roles and file content without rehashing multi-GB references during preflight.
  resources <- Neoantimon2::neo_check_resources(config = cfg, checksum = FALSE)
  if (!isTRUE(resources$ok)) stop('Reference check failed; inspect neo_check_resources(config=cfg); do not download automatically')
  inputs <- lapply(seq_along(samples), function(i) list(input_id = paste0(ids[[i]], '-somatic'),
    input_role = 'somatic_vcf', sample_id = ids[[i]], path = paths[[i]]))
  project <- Neoantimon2::neo_project(samples = data.frame(sample_id = ids, patient_id = patients),
    hla = hla, inputs = inputs, config = cfg)
  list(project = project, config = cfg, output = output, paths = c(config_file, paths), resources = resources)
}
agent_execute <- function(manifest_path, execute = FALSE) {
  manifest_path <- agent_absolute(manifest_path, 'manifest path')
  m <- yaml::read_yaml(manifest_path, eval.expr = FALSE)
  d <- agent_declare(m)
  cat('Preflight passed:', length(d$project$samples), 'samples. New run directory:', d$output, '\n')
  if (!execute) return(invisible(d))
  # Obtain authorization before calling with execute=TRUE. No downloads/installers.
  paths <- unique(c(manifest_path, d$paths, vapply(d$project$resources, `[[`, '', 'resolved_path')))
  hashes <- vapply(paths, function(p) digest::digest(file = p, algo = 'sha256', serialize = FALSE), '')
  if (!dir.create(d$output, recursive = TRUE)) stop('Cannot create fresh output directory')
  saveRDS(d$config, file.path(d$output, 'config.rds'))
  writeLines(capture.output(utils::sessionInfo()), file.path(d$output, 'sessionInfo.txt'))
  helper_path <- normalizePath(file.path('docs', 'ai', 'run-analysis.R'), mustWork = FALSE)
  writeLines(paste('Rscript --vanilla', shQuote(helper_path), shQuote(manifest_path), '--execute'),
    file.path(d$output, 'invocation.txt'))
  p <- Neoantimon2::neo_run(d$project)
  events <- Neoantimon2::neo_events(p); origins <- Neoantimon2::neo_origins(p); results <- Neoantimon2::neo_results(p)
  saveRDS(list(events = events, origins = origins, results = results), file.path(d$output, 'results.rds'))
  # Sensitive per-candidate information stays local in the serialized output.
  if (nrow(results)) saveRDS(Neoantimon2::neo_explain(p, candidate_id = results$candidate_id[[1L]]),
    file.path(d$output, 'explanation.rds'))
  reopened <- Neoantimon2::neo_project(path = p$path)
  stopifnot(identical(Neoantimon2::neo_results(reopened), results),
    identical(hashes, vapply(paths, function(p) digest::digest(file = p, algo = 'sha256', serialize = FALSE), '')))
  report <- list(package_version = as.character(utils::packageVersion('Neoantimon2')),
    package_source_commit = utils::packageDescription('Neoantimon2')$RemoteSha,
    resources = lapply(d$project$resources, function(x) x[c('resource_id', 'version', 'integrity')]),
    predictors = d$config$predictors, input_sha256 = as.list(hashes), project_path = p$path,
    events = nrow(events), origins = nrow(origins), candidates = nrow(results),
    prediction_status = as.list(table(results$prediction_status)), grading_status = as.list(table(results$grading_status)),
    inputs_unchanged = TRUE, reopened_identical = TRUE,
    note = 'NOT_RUN is unavailable prediction, NOT_GRADED is not grade E; research outputs, not clinical validation')
  jsonlite::write_json(report, file.path(d$output, 'report.json'), pretty = TRUE, auto_unbox = TRUE, null = 'null')
  cat('Completed; results, explanations and provenance retained locally.\n')
  invisible(p)
}
if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  if (!length(args) || length(args) > 2L || (length(args) == 2L && args[[2L]] != '--execute'))
    stop('Usage: Rscript docs/ai/run-analysis.R MANIFEST.yaml [--execute]')
  agent_execute(args[[1L]], execute = length(args) == 2L)
}
