# T003 infrastructure contract; independent of the frozen biological ID schema.
infra_version <- function() "1.0.0"
scalar_text <- function(x, label) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x))
    stop(label, " must be a nonempty string", call. = FALSE)
  invisible(x)
}
strict_keys <- function(x, allowed, label) {
  if (!is.list(x) || (length(x) && (is.null(names(x)) || anyDuplicated(names(x)) ||
      any(!nzchar(names(x)))))) stop(label, " must be a unique named list")
  if (length(setdiff(names(x), allowed))) stop("Unknown key in ", label)
}
merge_config <- function(base, x) {
  strict_keys(x, c(names(base), "allo", "indirect_hla", "allo_native_evidence"), "config")
  for (k in names(x)) {
    if (k %in% c("paths", "scientific", "runtime")) {
      strict_keys(x[[k]], names(base[[k]]), k)
      for (n in names(x[[k]])) base[[k]][n] <- x[[k]][n]
    } else base[k] <- x[k]
  }
  base
}

neo_config <- function(file = NULL, ..., env = Sys.getenv()) {
  user_data <- tools::R_user_dir("Neoantimon2", "data")
  cfg <- list(paths = list(data_root = file.path(user_data, "inputs"),
    run_root = file.path(user_data, "runs"), repository = NULL, input_root = NULL,
    tools_root = NULL), scientific = list(genome_assembly = NULL,
    species = "NCBITaxon:9606", policy_id = NULL),
    runtime = list(threads = 1, verbosity = "info"), checksum = "metadata",
    resources = list(), peptides = NULL, predictors = NULL, somatic_evidence = NULL,
    somatic_priority = NULL, hla_loh = NULL, vaccine_set = NULL)
  base <- NULL
  if (!is.null(file)) {
    file <- resolve_path(file)
    base <- dirname(file)
    x <- if (grepl("\\.ya?ml$", file, ignore.case = TRUE))
      yaml::read_yaml(file, eval.expr = FALSE) else jsonlite::read_json(file)
    reject_duplicate_keys(x)
    cfg <- merge_config(cfg, x)
    # File-relative root declarations are anchored before env/argument overlays.
    for (n in names(cfg$paths)) if (!is.null(cfg$paths[[n]]))
      cfg$paths[[n]] <- resolve_path(cfg$paths[[n]], base)
  }
  mapping <- c(DATA_ROOT = "data_root", RUN_ROOT = "run_root",
               REPOSITORY = "repository", INPUT_ROOT = "input_root", TOOLS_ROOT = "tools_root")
  for (n in names(mapping)) {
    v <- unname(env[paste0("NEOANTIMON2_", n)])
    if (!is.null(v) && !is.na(v) && nzchar(v)) cfg$paths[[mapping[[n]]]] <- v
  }
  for (n in c("THREADS", "VERBOSITY", "CHECKSUM")) {
    v <- unname(env[paste0("NEOANTIMON2_", n)])
    if (!is.null(v) && !is.na(v) && nzchar(v)) {
      if (n == "CHECKSUM") cfg$checksum <- v
      else cfg$runtime[[tolower(n)]] <- if (n == "THREADS") suppressWarnings(as.numeric(v)) else v
    }
  }
  cfg <- merge_config(cfg, list(...))
  for (n in names(cfg$paths)) if (!is.null(cfg$paths[[n]])) cfg$paths[[n]] <- resolve_path(cfg$paths[[n]])
  if (is.null(cfg$paths$input_root)) cfg$paths$input_root <- cfg$paths$data_root
  if (is.null(cfg$paths$tools_root)) cfg$paths$tools_root <- file.path(cfg$paths$data_root, "tools")
  for (n in c("data_root", "run_root")) scalar_text(cfg$paths[[n]], n)
  for (n in names(cfg$scientific)) if (!is.null(cfg$scientific[[n]])) scalar_text(cfg$scientific[[n]], n)
  scalar_text(cfg$scientific$species, "species")
  if (!grepl("^NCBITaxon:[1-9][0-9]*$", cfg$scientific$species)) stop("Invalid species")
  t <- cfg$runtime$threads
  if (!is.numeric(t) || length(t) != 1L || is.na(t) || !is.finite(t) || t < 1 || t != floor(t)) stop("Invalid threads")
  cfg$runtime$threads <- as.numeric(t)
  scalar_text(cfg$runtime$verbosity, "verbosity")
  if (!cfg$runtime$verbosity %in% c("quiet", "info", "debug")) stop("Invalid verbosity")
  scalar_text(cfg$checksum, "checksum")
  if (!cfg$checksum %in% c("none", "metadata", "sha256")) stop("Invalid checksum mode")
  cfg$resources <- validate_resources(cfg$resources)
  if (!is.null(cfg$peptides)) cfg$peptides <- peptide_policy(cfg$peptides)
  if (!is.null(cfg$predictors)) cfg$predictors <- predictor_policy(cfg$predictors)
  if (!is.null(cfg$somatic_evidence)) cfg$somatic_evidence <- somatic_evidence_policy(cfg$somatic_evidence)
  priority_input <- cfg$somatic_priority
  if (is.null(priority_input)) priority_input <- list()
  if (is.null(priority_input$rna_min_alt_reads) && !is.null(cfg$somatic_evidence))
    priority_input$rna_min_alt_reads <- cfg$somatic_evidence$eligibility$rna_min_alt_reads
  cfg$somatic_priority <- somatic_priority_policy(priority_input)
  if (!is.null(cfg$hla_loh)) cfg$hla_loh <- hla_loh_config(cfg$hla_loh)
  if (!is.null(cfg$vaccine_set)) cfg$vaccine_set <- vaccine_policy_bind(cfg$vaccine_set, base)
  else cfg$vaccine_set <- NULL
  if (!is.null(cfg[["allo"]])) cfg$allo <- allo_config(cfg[["allo"]], base)
  else cfg$allo <- NULL
  if (!is.null(cfg$indirect_hla)) {
    if (is.null(cfg[["allo"]])) stop("Indirect HLA baseline requires opt-in allo foundation")
    cfg$indirect_hla <- indirect_config(cfg$indirect_hla, base)
  } else cfg$indirect_hla <- NULL
  if (!is.null(cfg$allo_native_evidence)) {
    if (is.null(cfg$indirect_hla)) stop("Native allo evidence requires frozen indirect HLA baseline")
    cfg$allo_native_evidence <- native_allo_config(cfg$allo_native_evidence, base)
  }
  structure(cfg, class = c("neo_config", "list"))
}

print.neo_config <- function(x, ...) {
  cat("<neo_config> contract ", infra_version(), "; ", length(x$resources), " resources\n", sep = "")
  invisible(x)
}

validate_resources <- function(resources) {
  strict_keys(resources, names(resources), "resources")
  for (n in names(resources)) {
    if (!grepl("^(reference|tools|models|legacy|public_datasets|external_outputs)\\.[A-Za-z0-9_.-]+$", n))
      stop("Invalid logical resource name")
    resources[[n]] <- resource_spec(resources[[n]])
  }
  resources
}

# jsonlite retains duplicate object names; reject them before any overlay or
# $ lookup can hide them. yaml rejects duplicate map keys during parsing.
reject_duplicate_keys <- function(x) {
  if (is.list(x)) {
    if (!is.null(names(x)) && anyDuplicated(names(x))) stop("Duplicate config key")
    for (value in x) reject_duplicate_keys(value)
  }
  invisible(TRUE)
}
