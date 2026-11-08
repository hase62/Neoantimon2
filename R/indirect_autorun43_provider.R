# Additive 4.3k execution/parser. The legacy ia_version/ia_parse/ia_import stay intact.
ia43_version <- function(lines, requested = "4.3k") {
  if (requested == "3.2") return(ia_version(lines))
  hits <- grep("NetMHCIIpan[[:space:]]+version[[:space:]]+[0-9]", lines, value = TRUE)
  values <- unique(sub(".*NetMHCIIpan[[:space:]]+version[[:space:]]+([^[:space:]]+).*", "\\1", hits))
  if (length(values) != 1L || !(identical(values, "4.3k") ||
      (requested == "auto" && grepl("^3\\.2($|\\.)", values))))
    stop("Require exact NetMHCIIpan 4.3k; select predictor_version='3.2' for the legacy adapter. Unpatched/other 4.x releases are not substituted")
  values[[1L]]
}
ia43_model <- function(x) {
  allo_keys(x, c("path", "identity"), c("sha256", "runtime_subdir"), "4.3k model root")
  scalar_text(x$identity, "model identity"); scalar_text(x$path, "model path")
  root <- normalizePath(x$path, winslash = "/", mustWork = TRUE)
  if (!dir.exists(file.path(root, "data"))) stop("4.3k model root must contain data/")
  links <- list()
  walk <- function(relative = "") {
    entries <- list.files(file.path(root, relative), all.files = TRUE, no.. = TRUE)
    paths <- if (nzchar(relative)) file.path(relative, entries) else entries
    files <- character()
    for (p in paths) {
      full <- file.path(root, p); target <- Sys.readlink(full)
      if (nzchar(target)) {
        real <- normalizePath(full, winslash = "/", mustWork = TRUE)
        if (!identical(real, root) && !startsWith(real, paste0(root, "/"))) stop("Model symlink escapes pinned installation root")
        links[[length(links) + 1L]] <<- list(path = p, sha256 = digest::digest(target, algo = "sha256", serialize = FALSE), link_target = target)
      } else if (dir.exists(full)) files <- c(files, walk(p)) else files <- c(files, p)
    }
    files
  }
  files <- sort(walk(), method = "radix")
  if (!length(files) || !any(startsWith(files, "data/"))) stop("4.3k physical model data missing")
  records <- c(lapply(files, function(p) list(path = p, sha256 = ia_hash(file.path(root, p)))), links)
  records <- records[order(vapply(records, `[[`, "", "path"), method = "radix")]
  hash <- infra_hash(records)
  if (!is.null(x$sha256) && !identical(x$sha256, hash)) stop("Model tree SHA-256 mismatch")
  platform <- paste(Sys.info()[["sysname"]], Sys.info()[["machine"]], sep = "_")
  automatic <- file.path(root, platform)
  runtime <- if (is.null(x$runtime_subdir)) {
    if (dir.exists(file.path(automatic, "bin"))) normalizePath(automatic, winslash = "/", mustWork = TRUE) else root
  } else {
    scalar_text(x$runtime_subdir, "runtime_subdir")
    normalizePath(file.path(root, x$runtime_subdir), winslash = "/", mustWork = TRUE)
  }
  if (!identical(runtime, root) && !startsWith(runtime, paste0(root, "/"))) stop("runtime_subdir must remain inside pinned installation")
  if (!dir.exists(file.path(runtime, "data"))) stop("Runtime root must expose its matching data/ directory")
  list(path = root, rdir = runtime, identity = x$identity, tree_sha256 = hash, files = records,
    link_policy = "record internal links without traversal; physical targets hashed once")
}
ia43_probe <- function(executable, model, presenters, root, timeout, mock, requested, inversion) {
  h <- ia_exec(executable, "-h", root, "execution/help", timeout)
  if (!mock && any(grepl("NEOANTIMON2_TEST_MOCK", h$output, fixed = TRUE))) stop("Mock executable requires explicit mock_only authorization; not a real provider")
  help <- paste(h$output, collapse = "\n")
  tokens_help <- gsub("[", " ", gsub("]", " ", help, fixed = TRUE), fixed = TRUE)
  flag <- function(x) grepl(paste0("(^|[[:space:]])", x, "([[:space:]=,]|$)"), tokens_help)
  if (!all(vapply(c("-f", "-a", "-inptype", "-rdir"), flag, TRUE))) stop("Provider help lacks required local peptide/rdir CLI options")
  listing <- if (flag("-listMHC")) "-listMHC" else if (flag("-list")) "-list" else stop("Cannot establish exact supported alleles from provider CLI")
  l <- ia_exec(executable, c("-rdir", model$rdir, listing), root, "execution/alleles", timeout)
  version <- ia43_version(c(h$output, l$output), requested)
  if (version == "4.3k") {
    if (!flag("-BA")) stop("4.3k BA output must be explicitly supported/requested")
    for (f in c("context", "inv_all")) if (!grepl(paste0("\\[-", f, "\\][[:space:]]+0([[:space:]]|$)"), help))
      stop("Cannot establish 4.3k default context/inversion settings from actual help")
    if (inversion == "all" && !flag("-inv_all")) stop("Provider cannot explicitly enable inversion for all loci")
  } else if (inversion != "vendor_default") stop("Legacy 3.2 profile cannot use 4.3k inversion settings")
  tokens <- unique(unlist(strsplit(paste(l$output, collapse = " "), "[[:space:],]+")))
  required <- vapply(presenters, ia_vendor_allele, "")
  if (any(!required %in% tokens)) stop("Unsupported declared DR allele(s): ", paste(required[!required %in% tokens], collapse = ", "))
  list(version = version, alleles = as.list(required), commands = list(h$record, l$record),
    temp_flag = flag("-tdir"), context_encoding = FALSE, inversion = inversion,
    version_source = "executed help/list; each prediction banner checked again")
}
ia43_parse <- function(lines, peptide_set, presenter, vendor = ia_vendor_allele(presenter), require_ba = TRUE) {
  version <- ia43_version(lines)
  headers <- NULL; rows <- list(); has_ba <- FALSE
  base_fields <- c("Pos", "MHC", "Peptide", "Of", "Core", "Core_Rel", "Inverted", "Identity", "Score_EL", "%Rank_EL")
  ba_fields <- c("Score_BA", "%Rank_BA", "Affinity(nM)")
  for (i in seq_along(lines)) {
    text <- trimws(lines[[i]])
    if (!nzchar(text) || startsWith(text, "#") || grepl("^-+$", text)) next
    if (grepl("^Number of strong binders: [0-9]+ Number of weak binders: [0-9]+$", text)) next
    fields <- strsplit(text, "[[:space:]]+")[[1L]]
    if (fields[[1L]] == "Pos") {
      if (anyDuplicated(fields) || !all(base_fields %in% fields) ||
          any(!fields %in% c(base_fields, ba_fields, "Exp_Bind", "BindLevel"))) stop("Unexpected/missing 4.3k output columns")
      if (any(ba_fields %in% fields) && !all(ba_fields %in% fields)) stop("Incomplete 4.3k BA fields")
      has_ba <- all(ba_fields %in% fields)
      if (require_ba && !has_ba) stop("4.3k BA analysis requires Affinity(nM), %Rank_BA and Score_BA; EL-only is not BA")
      if (!is.null(headers) && !identical(headers, fields)) stop("Inconsistent 4.3k header")
      headers <- fields; next
    }
    if (is.null(headers) || !grepl("^[0-9]+$", fields[[1L]])) stop("Unexpected/error/truncated 4.3k output line ", i)
    base <- headers[headers != "BindLevel"]
    if (length(fields) < length(base)) stop("Truncated 4.3k row")
    tail <- if (length(fields) > length(base)) fields[seq.int(length(base) + 1L, length(fields))] else character()
    if (length(tail) && !paste(tail, collapse = "") %in% c("<=SB", "<=WB")) stop("Unexpected 4.3k row tail")
    c <- stats::setNames(as.list(fields[seq_along(base)]), base)
    if (!identical(c$MHC, vendor) || !c$Peptide %in% peptide_set) stop("Unexpected 4.3k peptide or allele; no substitution")
    if (!grepl("^[0-6]$", c$Of) || !c$Inverted %in% c("0", "1")) stop("Invalid 4.3k offset/Inverted")
    offset <- as.integer(c$Of); inverted <- as.integer(c$Inverted)
    start <- ia43_core_start(c$Peptide, c$Core, offset, inverted)
    ia43_number(c$Core_Rel, 1); ia43_number(c$Score_EL, 1); ia43_number(c[["%Rank_EL"]], 100)
    affinity <- rank <- NULL
    if (has_ba) { affinity <- ia43_number(c[["Affinity(nM)"]]); rank <- ia43_number(c[["%Rank_BA"]], 100); ia43_number(c$Score_BA, 1) }
    rows[[length(rows) + 1L]] <- list(peptide_sequence = c$Peptide, presenter = presenter,
      status = if (has_ba) "evaluated" else "el_only", affinity_nM = affinity, percentile_rank = rank,
      cores = list(list(core_sequence = c$Core, offset = offset,
        orientation = if (inverted == 1L) "reverse" else "forward", register = paste0("primary:", if (inverted == 1L) "reverse:" else "forward:", offset),
        primary = TRUE, status = "evaluated")),
      raw_output = list(line_number = i, line = lines[[i]], columns = c, provider_version = version,
        prediction_mode = if (has_ba) "EL+BA" else "EL", forward_ligand_core_start0 = start,
        offset_frame = "reported Of indexes ligand reversed first when Inverted=1"), reason_codes = list())
  }
  if (is.null(headers) || !setequal(vapply(rows, `[[`, "", "peptide_sequence"), peptide_set)) stop("Missing 4.3k results, including required self queries")
  for (p in peptide_set) {
    rr <- Filter(function(r) r$peptide_sequence == p, rows)
    signatures <- vapply(rr, function(r) canonical_bytes(list(binding = r[setdiff(names(r), "raw_output")],
      metrics = if (has_ba) ia43_row_metrics(r) else r$raw_output$columns)), "")
    if (length(unique(signatures)) != 1L) stop("Conflicting duplicate 4.3k query/metrics/core")
  }
  rows
}
ia43_import <- function(rows, query_path, result_path, executable, model, probe, commands, mock) {
  x <- ia_import(rows, query_path, result_path, executable, model, probe, commands, mock)
  x$version <- ia43_import_version()
  x$provider$arguments <- list(BA = TRUE, context_encoding = FALSE, inversion = probe$inversion)
  x$provider$calibration <- "Uncalibrated cross-version BA cutoff assumption; EL retained separately and not used as nM"
  if (mock) x$provenance$fixture_id <- "autorun43k-controlled-mock-v1"
  x
}
ia43_metrics_table <- function(rows) lapply(rows, function(r) {
  c <- r$raw_output$columns
  list(peptide_sequence = r$peptide_sequence, presenter = r$presenter, MHC = c$MHC, Core = c$Core,
    Core_Rel = c$Core_Rel, Inverted = c$Inverted, Of = c$Of, orientation = r$cores[[1L]]$orientation,
    forward_ligand_core_start0 = r$raw_output$forward_ligand_core_start0,
    affinity_nM = c[["Affinity(nM)"]], rank_BA = c[["%Rank_BA"]], score_BA = c$Score_BA,
    score_EL = c$Score_EL, rank_EL = c[["%Rank_EL"]], version = r$raw_output$provider_version)
})
