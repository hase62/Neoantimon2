# Optional orchestration only; the existing molecular/counting contracts are reused.
ia_hash <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
ia_json <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(x, path, auto_unbox = TRUE, null = "null", digits = NA, pretty = TRUE)
  invisible(path)
}
ia_table <- function(rows, path) {
  fields <- sort(unique(unlist(lapply(rows, names))), method = "radix")
  if (!length(fields)) { writeLines(character(), path); return(invisible(path)) }
  cell <- function(x) if (is.null(x)) "" else if (is.list(x) || length(x) != 1L)
    as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", digits = NA)) else as.character(x)
  out <- as.data.frame(stats::setNames(lapply(fields, function(k)
    vapply(rows, function(r) cell(r[[k]]), "")), fields), stringsAsFactors = FALSE)
  utils::write.table(out, path, sep = "\t", row.names = FALSE, quote = TRUE, na = "")
  invisible(path)
}
ia_model <- function(x) {
  allo_keys(x, c("path", "identity"), "sha256", "local NetMHCIIpan model root")
  scalar_text(x$path, "model path"); scalar_text(x$identity, "model identity")
  root <- normalizePath(x$path, winslash = "/", mustWork = TRUE)
  if (!dir.exists(root) || !dir.exists(file.path(root, "data")))
    stop("Model root must be the installed 3.2 -rdir with its data/ directory; no automatic vendor download")
  # Inspect each directory entry before descending: never follow a model symlink.
  walk <- function(relative = "") {
    entries <- list.files(file.path(root, relative), all.files = TRUE, no.. = TRUE)
    paths <- if (nzchar(relative)) file.path(relative, entries) else entries
    if (any(nzchar(Sys.readlink(file.path(root, paths)))))
      stop("Model pinning requires a physical tree without symlinks; supply an explicit physical installation root")
    info <- file.info(file.path(root, paths))
    if (anyNA(info$isdir)) stop("Cannot inspect model directory entries")
    c(paths[!info$isdir], unlist(lapply(paths[info$isdir], walk), use.names = FALSE))
  }
  files <- sort(walk(), method = "radix")
  if (!length(files) || !any(startsWith(files, "data/"))) stop("NetMHCIIpan model data is unavailable")
  records <- lapply(files, function(p) list(path = p, sha256 = ia_hash(file.path(root, p))))
  hash <- infra_hash(records)
  if (!is.null(x$sha256) && !identical(x$sha256, hash)) stop("Model tree SHA-256 mismatch")
  list(path = root, identity = x$identity, tree_sha256 = hash, files = records)
}
ia_code_hash <- function() {
  env <- environment(ia_code_hash)
  names <- sort(c(grep("^ia(_|43_)", ls(env, all.names = TRUE), value = TRUE),
    "neo_indirect_autorun", "indirect_baseline", "indirect_context", "indirect_generate", "indirect_self_catalog", "indirect_provider",
    "indirect_provider_threshold", "indirect_core_key", "indirect_primary", "indirect_compare", "indirect_self_prediction"))
  infra_hash(lapply(names, function(n) {
    f <- get(n, envir = env); list(name = n, formals = as.list(deparse(formals(f))), body = as.list(deparse(body(f))))
  }))
}
ia_files <- function(root) {
  files <- sort(list.files(root, recursive = TRUE, all.files = TRUE, no.. = TRUE), method = "radix")
  files <- setdiff(files, c("run.json", "SHA256SUMS.txt", ".autorun-lock/owner"))
  files <- files[!file.info(file.path(root, files))$isdir]
  lapply(files, function(p) list(path = p, sha256 = ia_hash(file.path(root, p))))
}
ia_verify_files <- function(root, manifest) {
  current <- ia_files(root)
  if (!identical(canonical_bytes(current), canonical_bytes(manifest$files)))
    stop("Saved autorun artifact checksum/file-set mismatch; original files are preserved")
  checksum_path <- file.path(root, "SHA256SUMS.txt")
  records <- c(manifest$files, list(list(path = "run.json", sha256 = ia_hash(file.path(root, "run.json")))))
  expected <- vapply(records, function(r) paste(r$sha256, r$path, sep = "  "), "")
  if (!file.exists(checksum_path) || !identical(readLines(checksum_path, warn = FALSE), unname(expected)))
    stop("Saved autorun checksum manifest mismatch")
  invisible(TRUE)
}
ia_finish <- function(root, manifest) {
  manifest$files <- ia_files(root)
  ia_json(manifest, file.path(root, "run.json"))
  records <- c(manifest$files, list(list(path = "run.json", sha256 = ia_hash(file.path(root, "run.json")))))
  writeLines(vapply(records, function(r) paste(r$sha256, r$path, sep = "  "), ""), file.path(root, "SHA256SUMS.txt"))
  manifest
}
ia_return <- function(root) {
  read <- function(n) if (file.exists(file.path(root, n))) jsonlite::read_json(file.path(root, n)) else NULL
  structure(list(status = read("run.json")$status, output_dir = root, manifest = read("run.json"),
    plan = read("plan.json"), result = read("result.json"), explanation = read("explanation.json")),
    class = "neo_indirect_autorun")
}
ia_exec <- function(executable, args, root, stem, timeout) {
  if (any(grepl("[\r\n]", c(executable, args)))) stop("Newline in executable argument")
  stdout <- file.path(root, paste0(stem, ".stdout")); stderr <- file.path(root, paste0(stem, ".stderr"))
  dir.create(dirname(stdout), recursive = TRUE, showWarnings = FALSE)
  # Base R quotes the command; each argument is separately shell-quoted. No
  # user-supplied command string, sh -c, redirection or command interpolation.
  status <- tryCatch(suppressWarnings(system2(executable, args = vapply(args, shQuote, ""),
    stdout = stdout, stderr = stderr, timeout = timeout,
    env = c("OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1", "MKL_NUM_THREADS=1"))), error = function(e) {
      writeLines(conditionMessage(e), stderr); 127L
    })
  rec <- list(executable = executable, argv = as.list(args), exit_status = as.integer(status),
    stdout = paste0(stem, ".stdout"), stderr = paste0(stem, ".stderr"),
    environment = list(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"))
  ia_json(rec, file.path(root, paste0(stem, ".command.json")))
  if (status != 0L) stop("Local predictor command failed (exit ", status, "): ", stem, "; inspect retained stdout/stderr")
  list(record = rec, output = readLines(stdout, warn = FALSE))
}
