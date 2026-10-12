resolve_path <- function(path, base = NULL) {
  scalar_text(path, "path")
  path <- path.expand(path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    if (is.null(base)) stop("Relative path requires an explicit absolute base")
    path <- file.path(base, path)
  }
  # Resolve existing ancestors as well as full paths (including symlinks).
  if (file.exists(path)) return(normalizePath(path, winslash = "/", mustWork = TRUE))
  link <- Sys.readlink(path)
  if (!is.na(link) && nzchar(link)) stop("Broken symlink")
  parent <- dirname(path)
  if (identical(parent, path)) stop("Cannot resolve path")
  parent <- resolve_path(parent)
  leaf <- basename(path)
  if (leaf == "..") return(dirname(parent))
  if (leaf == ".") return(parent)
  file.path(parent, leaf)
}
within_path <- function(path, root) {
  # Preserve filesystem roots while removing trailing separators.
  trim <- function(x) {
    x <- sub("/+$", "", x)
    if (!nzchar(x)) "/" else x
  }
  path <- trim(path); root <- trim(root)
  identical(path, root) || startsWith(path, if (root == "/") "/" else paste0(root, "/"))
}

resource_spec <- function(x) {
  if (is.character(x)) x <- list(path = x)
  default <- list(path = NULL, required = FALSE, type = "file", resource_id = NULL,
                  version = NULL, checksum = NULL, source = NULL)
  strict_keys(x, names(default), "resource")
  for (n in names(x)) default[n] <- x[n]
  x <- default
  if (!is.logical(x$required) || length(x$required) != 1L || is.na(x$required)) stop("Invalid required flag")
  scalar_text(x$type, "resource type")
  if (!x$type %in% c("file", "directory", "executable", "model")) stop("Invalid resource type")
  for (n in c("path", "resource_id", "version", "source")) if (!is.null(x[[n]])) scalar_text(x[[n]], n)
  if (!is.null(x$source) && grepl("(@|[?]|token|password|secret)", x$source, ignore.case = TRUE))
    stop("Credential-bearing source is not allowed")
  if (!is.null(x$checksum)) {
    strict_keys(x$checksum, c("algorithm", "scope", "value"), "checksum")
    if (!identical(x$checksum$algorithm, "sha256") || !identical(x$checksum$scope, "file") ||
        !is.character(x$checksum$value) || length(x$checksum$value) != 1L ||
        is.na(x$checksum$value) || !grepl("^[a-f0-9]{64}$", x$checksum$value)) stop("Invalid declared checksum")
    if (x$type == "directory") stop("Directory checksum requires a future explicit tree contract")
  }
  x
}

integrity_record <- function(path, mode, supplied = NULL) {
  exists <- !is.null(path) && file.exists(path)
  observed <- list(algorithm = "none", scope = "none", value = NULL)
  if (exists && mode != "none") {
    info <- file.info(path)
    if (mode == "sha256") {
      if (info$isdir) stop("Directory sha256 is unsupported; provide a file inventory")
      observed <- list(algorithm = "sha256", scope = "file",
        value = digest::digest(file = path, algo = "sha256", serialize = FALSE))
      if (!is.null(supplied) && !identical(supplied$value, observed$value)) stop("Declared checksum mismatch")
    } else observed <- list(algorithm = "stat", scope = if (info$isdir) "directory-entry-only" else "file-metadata",
      value = list(size = as.numeric(info$size), mtime = sprintf("%.6f", as.numeric(info$mtime))))
  }
  list(mode = mode, exists = exists, declared = supplied, observed = observed,
       verified = exists && mode == "sha256" && !is.null(supplied))
}

resolve_resources <- function(config) {
  out <- lapply(names(config$resources), function(n) {
    s <- config$resources[[n]]
    base <- if (grepl("^(tools|models)\\.", n)) config$paths$tools_root else config$paths$data_root
    path <- if (is.null(s$path)) NULL else resolve_path(s$path, base)
    exists <- !is.null(path) && file.exists(path)
    if (s$required && !exists) stop("Required resource missing: ", n)
    if (exists && (dir.exists(path) != (s$type == "directory"))) stop("Resource type mismatch: ", n)
    c(list(logical_name = n, declared_path = s$path, resolved_path = path), s[setdiff(names(s), "path")],
      list(integrity = integrity_record(path, config$checksum, s$checksum)))
  })
  out[order(vapply(out, `[[`, "", "logical_name"), method = "radix")]
}

infra_hash <- function(x) {
  # JSON round trips may choose integer vs double; T003 numbers are binary64.
  normalize <- function(v) {
    if (is.list(v)) return(lapply(v, normalize))
    if (is.integer(v)) return(as.numeric(v))
    v
  }
  digest::digest(charToRaw(canonical_bytes(normalize(x))), algo = "sha256", serialize = FALSE)
}
sorted_records <- function(x) {
  if (!length(x)) return(list())
  x[order(vapply(x, canonical_bytes, ""), method = "radix")]
}
