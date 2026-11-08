ia_version <- function(lines) {
  hit <- grep("NetMHCIIpan[[:space:]]+(version[[:space:]]+)?[0-9]+\\.[0-9]+", lines, value = TRUE)
  if (!length(hit)) stop("Cannot establish NetMHCIIpan version from actual help/output")
  values <- unique(sub(".*NetMHCIIpan[[:space:]]+(version[[:space:]]+)?([0-9]+\\.[0-9]+(\\.[0-9]+)?).*", "\\2", hit))
  if (length(values) != 1L || !grepl("^3\\.2($|\\.)", values)) stop("Strict autorun requires NetMHCIIpan 3.2 BA; other versions are not substituted")
  values[[1L]]
}
ia_vendor_allele <- function(presenter) {
  p <- peptide_presenters(list(presenter))[[1L]]
  # NetMHCIIpan 3.2 DR queries expose beta only; alpha is not freely variable.
  if (p$status != "RESOLVED" || p$hla_class != "II" ||
      !identical(p$hla_alpha, "HLA-DRA*01:01") ||
      !grepl("^HLA-DRB1\\*[0-9]{2}:[0-9]{2}$", p$hla_beta))
    stop("Unsupported strict autorun presenter: require explicit DRA*01:01~DRB1 with two-field beta; no truncation or substitution: ", presenter)
  paste0("DRB1_", gsub(":", "", sub("^HLA-DRB1\\*", "", p$hla_beta)))
}
ia_probe <- function(executable, model, presenters, root, timeout, mock) {
  h <- ia_exec(executable, "-h", root, "execution/help", timeout)
  version <- ia_version(h$output)
  if (!mock && any(grepl("NEOANTIMON2_TEST_MOCK", h$output, fixed = TRUE))) stop("Mock executable requires explicit mock_only authorization; not a real provider")
  help <- paste(h$output, collapse = "\n")
  flag <- function(x) grepl(paste0("(^|[[:space:]])", x, "([[:space:]=,]|$)"), help)
  if (!all(vapply(c("-f", "-a", "-inptype", "-rdir"), flag, TRUE)))
    stop("Unsupported 3.2 CLI: help must advertise -f, -a, -inptype and -rdir for explicit data binding")
  listing <- if (flag("-listMHC")) "-listMHC" else if (flag("-list")) "-list" else
    stop("Installed provider help does not expose -list/-listMHC; cannot verify DR support")
  l <- ia_exec(executable, c("-rdir", model$path, listing), root, "execution/alleles", timeout)
  tokens <- unique(unlist(strsplit(paste(l$output, collapse = " "), "[[:space:],]+")))
  required <- vapply(presenters, ia_vendor_allele, "")
  if (any(!required %in% tokens)) stop("Unsupported declared DR allele(s): ", paste(required[!required %in% tokens], collapse = ", "))
  list(version = version, alleles = as.list(required), commands = list(h$record, l$record),
    temp_flag = flag("-tdir"))
}
ia_parse <- function(lines, peptide_set, presenter) {
  version <- ia_version(lines)
  vendor <- ia_vendor_allele(presenter)
  headers <- NULL; layout <- NULL; rows <- list()
  aa <- "^[ACDEFGHIKLMNPQRSTVWY]+$"
  metric <- function(x, maximum = Inf) {
    if (!grepl("^([0-9]+(\\.[0-9]*)?|\\.[0-9]+)([eE][+-]?[0-9]+)?$", x)) stop("Invalid/missing numeric provider metric")
    v <- as.numeric(x)
    if (!is.finite(v) || v < 0 || v > maximum) stop("Invalid provider metric range")
    v
  }
  for (i in seq_along(lines)) {
    text <- trimws(lines[[i]])
    if (!nzchar(text) || startsWith(text, "#") || grepl("^-+$", text)) next
    fields <- strsplit(text, "[[:space:]]+")[[1L]]
    if (fields[[1L]] %in% c("Seq", "Pos")) {
      type <- if (fields[[1L]] == "Seq") "seq" else "pos"
      needed <- if (type == "seq") c("Seq", "Allele", "Peptide", "Identity", "Pos", "Core", "Core_Rel", "1-log50k(aff)", "Affinity(nM)") else
        c("Pos", "MHC", "Peptide", "Of", "Core", "Core_Rel", "Identity", "Score", "Affinity(nM)")
      optional <- if (type == "seq") c("%Rank", "Exp_Bind", "BindingLevel") else c("%Rank", "Exp_Bind", "BindLevel")
      if (anyDuplicated(fields) || !all(needed %in% fields) || any(!fields %in% c(needed, optional))) stop("Unexpected NetMHCIIpan 3.2 output columns")
      if (!is.null(headers) && !identical(headers, fields)) stop("Inconsistent provider header")
      headers <- fields; layout <- type; next
    }
    if (is.null(headers) || !grepl("^[0-9]+$", fields[[1L]])) stop("Unexpected/truncated provider output line ", i)
    base <- headers[!headers %in% c("BindingLevel", "BindLevel")]
    if (length(fields) < length(base)) stop("Truncated provider row")
    tail <- if (length(fields) > length(base)) fields[seq.int(length(base) + 1L, length(fields))] else character()
    if (length(tail) && !paste(tail, collapse = "") %in% c("<=SB", "<=WB")) stop("Unexpected provider row tail")
    r <- stats::setNames(as.list(fields[seq_along(base)]), base)
    allele <- r[[if (layout == "seq") "Allele" else "MHC"]]
    if (!identical(allele, vendor)) stop("Unexpected provider allele; no fuzzy HLA substitution")
    peptide <- r$Peptide
    if (!peptide %in% peptide_set || nchar(peptide) != 15L) stop("Unexpected provider peptide")
    of <- r[[if (layout == "seq") "Pos" else "Of"]]
    if (!grepl("^[0-6]$", of)) stop("Invalid zero-based core offset")
    offset <- as.integer(of); core <- r$Core
    if (nchar(core) != 9L || !grepl(aa, core) || substring(peptide, offset + 1L, offset + 9L) != core)
      stop("Ambiguous/invalid primary forward core mapping; no inferred core")
    metric(r$Core_Rel, 1)
    affinity <- metric(r[["Affinity(nM)"]]); rank <- if ("%Rank" %in% base) metric(r[["%Rank"]], 100) else NULL
    rows[[length(rows) + 1L]] <- list(peptide_sequence = peptide, presenter = presenter, status = "evaluated",
      affinity_nM = affinity, percentile_rank = rank,
      cores = list(list(core_sequence = core, offset = offset, orientation = "forward", register = paste0("primary:", offset), primary = TRUE, status = "evaluated")),
      raw_output = list(line_number = i, line = lines[[i]], columns = r, provider_version = version), reason_codes = list())
  }
  if (is.null(headers) || !setequal(vapply(rows, `[[`, "", "peptide_sequence"), peptide_set))
    stop("Missing provider results, including required self queries; truncated or incomplete output")
  for (peptide in peptide_set) {
    rr <- Filter(function(r) r$peptide_sequence == peptide, rows)
    if (length(unique(vapply(rr, function(r) canonical_bytes(r[setdiff(names(r), "raw_output")]), ""))) != 1L)
      stop("Conflicting duplicate/ambiguous primary core for provider query")
  }
  rows
}
ia_import <- function(rows, query_path, result_path, executable, model, probe, commands, mock) {
  build <- paste0("NetMHCIIpan-", probe$version, ";sha256:", executable$sha256)
  modelid <- paste0(model$identity, ";tree-sha256:", model$tree_sha256)
  artifacts <- list(query_artifact = list(identity = "exact_sorted_15mer_presenter_queries", reference = "queries.tsv", sha256 = ia_hash(query_path)),
    result_artifact = list(identity = "original_stdout_index", reference = "execution/raw-output-index.json", sha256 = ia_hash(result_path)),
    executable_artifact = list(identity = build, reference = paste0("external:sha256:", executable$sha256), sha256 = executable$sha256),
    model_artifact = list(identity = modelid, reference = "execution/model-tree.json", sha256 = ia_hash(file.path(dirname(result_path), "model-tree.json"))))
  provenance <- if (mock) list(kind = "synthetic_fixture", fixture_id = "autorun-controlled-mock-v0.1",
    license = "synthetic software test; no genuine predictor", source = "controlled mock stdout; not observed biological evidence",
    execution = "not_executed", metadata = c(artifacts, list(mock_process_executed = TRUE, real_predictor_executed = FALSE, argv = commands))) else
    c(list(kind = "precomputed_predictor_output"), artifacts, list(
      license = list(status = "authorized_external_precomputed", terms = executable$authorization$terms),
      execution = list(origin = "imported_external_execution", source = "explicit authorized local installation",
        provider_name = "NetMHCIIpan", provider_version = probe$version,
        command = canonical_bytes(commands)), metadata = list(local_authorization_status = executable$authorization$status)))
  list(version = "indirect-binding-import-1.0.0", provider = list(name = "NetMHCIIpan", version = probe$version,
    build = build, model = modelid, task = "mhc2_binding", prediction_mode = "BA", species = "NCBITaxon:9606",
    arguments = list(BA = TRUE), supported_presenters = as.list(names(probe$presenter_map)), supported_lengths = list(15L),
    calibration = "none; raw BA nM, not clinical probability", batching = "sorted unique ligand/presenter queries",
    evidence_kind = if (mock) "synthetic_fixture" else "precomputed_predictor_output"), rows = rows, provenance = provenance)
}
