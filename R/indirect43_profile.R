# Explicit new analysis policy. The 3.2 profile and molecular generation IDs stay frozen.
ia43_profile <- function() "neo-netmhciipan43k-indirect-hla-ii-dr-ba1000-orientation-v1"
ia43_import_version <- function() "indirect-binding-netmhciipan43k-ba-1.0.0"
ia43_is_provider <- function(p) identical(p$semantics$version, "4.3k")
ia43_is_analysis <- function(x) identical(x$profile, ia43_profile())
ia43_schema_names <- function() c(
  "indirect-binding-import-1.0.0" = ia43_import_version(),
  "indirect-hla-baseline-1.0.0" = "indirect-hla-netmhciipan43k-ba-1.0.0",
  "indirect-pmhc-evidence-1.0.0" = "indirect-pmhc-netmhciipan43k-ba-1.0.0",
  "indirect-hla-result-1.0.0" = "indirect-hla-result-netmhciipan43k-ba-1.0.0",
  "indirect-hla-explanation-1.0.0" = "indirect-hla-explanation-netmhciipan43k-ba-1.0.0")
indirect_provider_threshold <- function(provider) {
  x <- indirect_threshold()
  if (ia43_is_provider(provider)) {
    x$version <- "netmhciipan43k-ba1000-assumption-1.0.0"; x$profile <- ia43_profile()
    x$citation <- "Explicit exploratory 1000 nM BA threshold; not calibrated across versions or equivalent to classical PIRCHE-II v3"
  }
  x
}
indirect_core_key <- function(provider, core, presenter, species) {
  if (ia43_is_provider(provider)) return(indirect_id("indirect43-core-presenter-id-v1", list(
    profile = ia43_profile(), species = species, core_sequence = core$core_sequence,
    orientation = core$orientation, presenter = presenter)))
  biological_id("epitope_key", list(species = species, hla_class = "II",
    binding_core = core$core_sequence, mhc_key = presenter), peptide_contract())
}
ia43_arguments <- function(x) {
  allo_keys(x, c("BA", "context_encoding", "inversion"), label = "4.3k semantic arguments")
  if (!identical(x$BA, TRUE) || !identical(x$context_encoding, FALSE) ||
      !x$inversion %in% c("vendor_default", "all")) stop("4.3k BA profile requires explicit BA=true, context_encoding=false and inversion policy")
  x
}
ia43_number <- function(x, maximum = Inf) {
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      !grepl("^([0-9]+(\\.[0-9]*)?|\\.[0-9]+)([eE][+-]?[0-9]+)?$", x)) stop("Invalid/missing 4.3k numeric metric")
  n <- as.numeric(x)
  if (!is.finite(n) || n < 0 || n > maximum) stop("Invalid 4.3k metric range")
  n
}
ia43_core_start <- function(peptide, core, offset, inverted) {
  if (!grepl("^[ACDEFGHIKLMNPQRSTVWY]{15}$", peptide) ||
      !grepl("^[ACDEFGHIKLMNPQRSTVWY]{9}$", core) || !offset %in% 0:6 ||
      !inverted %in% c(0L, 1L)) stop("Invalid 4.3k core/offset/Inverted")
  oriented <- if (inverted == 1L) paste(rev(strsplit(peptide, "", fixed = TRUE)[[1L]]), collapse = "") else peptide
  if (substring(oriented, offset + 1L, offset + 9L) != core)
    stop("4.3k core does not map in reported orientation; Of is relative to the oriented ligand")
  if (inverted == 1L) 15L - 9L - as.integer(offset) else as.integer(offset)
}
ia43_row_metrics <- function(r) {
  c <- r$raw_output$columns
  list(core_rel = ia43_number(c$Core_Rel, 1), inverted = as.integer(c$Inverted),
    score_el = ia43_number(c$Score_EL, 1), rank_el = ia43_number(c[["%Rank_EL"]], 100),
    score_ba = ia43_number(c$Score_BA, 1))
}
ia43_validate_row <- function(r) {
  if (r$status != "evaluated") return(invisible(TRUE))
  c <- r$raw_output$columns
  needed <- c("Peptide", "MHC", "Of", "Core", "Core_Rel", "Inverted", "Score_EL", "%Rank_EL", "Score_BA", "%Rank_BA", "Affinity(nM)")
  if (!is.list(c) || !all(needed %in% names(c)) ||
      !identical(r$raw_output$provider_version, "4.3k")) stop("4.3k evaluated row requires original BA/EL/core fields and exact version")
  if (!identical(c$Peptide, r$peptide_sequence) || !identical(c$MHC, ia_vendor_allele(r$presenter)) ||
      !grepl("^[0-6]$", c$Of) || !c$Inverted %in% c("0", "1")) stop("4.3k row identity/offset/orientation mismatch")
  metrics <- ia43_row_metrics(r)
  if (!identical(as.numeric(r$affinity_nM), ia43_number(c[["Affinity(nM)"]])) ||
      is.null(r$percentile_rank) || !identical(as.numeric(r$percentile_rank), ia43_number(c[["%Rank_BA"]], 100)))
    stop("4.3k BA measurements disagree with raw fields; EL is not BA")
  if (length(r$cores) != 1L) stop("4.3k requires one reported primary core")
  core <- r$cores[[1L]]
  expected <- if (metrics$inverted == 1L) "reverse" else "forward"
  if (!isTRUE(core$primary) || core$status != "evaluated" || core$orientation != expected ||
      core$core_sequence != c$Core || core$offset != as.integer(c$Of)) stop("4.3k primary core conflicts with raw output")
  start <- ia43_core_start(r$peptide_sequence, core$core_sequence, core$offset, metrics$inverted)
  if (!identical(as.integer(r$raw_output$forward_ligand_core_start0), start)) stop("4.3k physical core coordinate mismatch")
  invisible(TRUE)
}
ia43_compatibility <- function() list(schema_version = "netmhciipan43k-compatibility-1.0.0",
  profile = ia43_profile(), status = "not_compared", official_reproduction = FALSE,
  classical_3_2_numerical_equivalence = FALSE, pirche_v4_frost_equivalence = FALSE,
  clinical_performance_established = FALSE)
ia43_public <- function(out, x) {
  if (!ia43_is_analysis(x)) return(out)
  versions <- ia43_schema_names()
  convert <- function(z) {
    if (!is.list(z)) return(z)
    for (k in seq_along(z)) {
      key <- if (is.null(names(z))) "" else names(z)[[k]]
      if (key == "schema_version" && is.character(z[[k]]) && z[[k]] %in% names(versions)) z[[k]] <- unname(versions[[z[[k]]]])
      else if (key == "specification_id" && identical(z[[k]], "neo-allo-indirect-spec-0.1")) z[[k]] <- "neo-netmhciipan43k-indirect-spec-0.1"
      else z[k] <- list(convert(z[[k]]))
    }
    z
  }
  out <- convert(out)
  orientation <- list()
  for (c in x$contexts) for (p in c$pmhcs) for (a in p$core_assignments)
    if (!is.null(a$epitope_key)) orientation[[a$epitope_key]] <- a$semantics$orientation
  for (i in seq_along(out$contexts)) {
    ctx <- out$contexts[[i]]
    if (!is.null(ctx$graph)) {
      for (type in c("epitope", "aggregate_unit")) for (j in seq_along(ctx$graph$nodes[[type]])) {
        node <- ctx$graph$nodes[[type]][[j]]
        key <- if (type == "epitope") node$id else node$data$counted_unit_key
        node$data$orientation <- orientation[[key]]; ctx$graph$nodes[[type]][[j]] <- node
      }
      for (j in seq_along(ctx$included_units)) ctx$included_units[[j]]$orientation <- orientation[[ctx$included_units[[j]]$counted_unit_key]]
    }
    out$contexts[[i]] <- ctx
  }
  out
}
