vaccine_inputs <- function(project) {
  binding <- vaccine_policy_bind(project$config$vaccine_set)
  if (length(binding$physical)) {
    p <- binding$physical$path
    if (!file.exists(p) || !identical(digest::digest(file = p, algo = "sha256", serialize = FALSE), binding$physical$sha256))
      stop("Vaccine policy file changed since binding")
  }
  list(binding = binding, projection = vaccine_projection(project))
}
vaccine_fingerprint <- function(inputs) infra_hash(list(schema = vaccine_selection_version(),
  algorithm = vaccine_algorithm(), policy_id = inputs$binding$policy_id,
  # This is an execution checkpoint digest, never a biological/set/policy ID.
  # RDS payload hashing avoids recursively canonicalizing large frozen graphs.
  bound_evidence = digest::digest(inputs, algo = "sha256", serialize = TRUE)))
vaccine_root <- function(project, inputs) file.path(project$path, "results", "vaccine_set", vaccine_fingerprint(inputs))
vaccine_validate <- function(x, inputs) {
  expected <- vaccine_construct(inputs$projection, inputs$binding)
  if (!identical(x, expected)) stop("Vaccine selection derived contract mismatch")
  invisible(TRUE)
}
vaccine_postrun <- function(project) {
  if (!vaccine_enabled(project)) return(invisible(NULL))
  inputs <- vaccine_inputs(project); root <- vaccine_root(project, inputs)
  if (file.exists(file.path(root, "checkpoint.json"))) return(invisible(vaccine_read(project)))
  parent <- dirname(root); dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  if (!within_path(resolve_path(parent), project$path)) stop("Unsafe vaccine output path")
  tmp <- tempfile(".vaccine-", tmpdir = parent); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  x <- vaccine_construct(inputs$projection, inputs$binding)
  vaccine_validate(x, inputs)
  saveRDS(x, file.path(tmp, "selection.rds"), version = 3L)
  atomic_json(longitudinal_plain(x), file.path(tmp, "selection.json"))
  atomic_json(list(schema_version = vaccine_selection_version(), fingerprint = vaccine_fingerprint(inputs),
    outputs = lapply(c("selection.rds", "selection.json"), function(f) list(file = f,
      sha256 = digest::digest(file = file.path(tmp, f), algo = "sha256", serialize = FALSE)))), file.path(tmp, "checkpoint.json"))
  if (dir.exists(root) || !file.rename(tmp, root)) stop("Could not publish vaccine overlay")
  invisible(x)
}
vaccine_read <- function(project) {
  if (!vaccine_enabled(project)) return(NULL)
  inputs <- vaccine_inputs(project); root <- vaccine_root(project, inputs)
  path <- file.path(root, "checkpoint.json")
  if (!file.exists(path)) stop("Vaccine selection sidecar missing; run neo_run")
  cp <- jsonlite::read_json(path)
  if (!identical(cp$schema_version, vaccine_selection_version()) || !identical(cp$fingerprint, vaccine_fingerprint(inputs)) ||
      !identical(sort(vapply(cp$outputs, `[[`, "", "file")), sort(c("selection.rds", "selection.json"))))
    stop("Vaccine checkpoint mismatch")
  for (o in cp$outputs) {
    p <- resolve_path(file.path(root, o$file))
    if (!within_path(p, root) || !file.exists(p) || !identical(o$sha256, digest::digest(file = p, algo = "sha256", serialize = FALSE)))
      stop("Vaccine selection integrity failure")
  }
  x <- readRDS(file.path(root, "selection.rds"))
  vaccine_validate(x, inputs)
  x
}
vaccine_summary <- function(x) list(policy_version = x$policy_version, policy_id = x$policy_id,
  algorithm = x$algorithm, sets = lapply(x$selections, function(s) s[c("vaccine_set_id", "scope", "status", "selected_count",
    "max_peptides", "class_i_count", "class_ii_count", "coverage", "unmet_constraints", "stop_reason")]))
vaccine_results_append <- function(out, project) {
  if (!vaccine_enabled(project)) return(out)
  attr(out, "vaccine_set") <- vaccine_summary(vaccine_read(project))
  out
}
vaccine_explain_append <- function(items, project) {
  if (!vaccine_enabled(project)) return(items)
  x <- vaccine_read(project)
  audits <- unlist(lapply(x$selections, `[[`, "candidate_audit"), recursive = FALSE)
  lapply(items, function(item) {
    matches <- Filter(function(a) a$candidate_id == item$result$candidate_id, audits)
    if (length(matches) != 1L) stop("Vaccine explanation candidate mismatch")
    item$vaccine_selection <- matches[[1L]]
    item
  })
}
vaccine_explain_text <- function(item) {
  x <- item$vaccine_selection
  if (is.null(x)) return("")
  paste0("\nVACCINE SET SELECTION: ", x$status, "; ", vaccine_algorithm(), " (deterministic lexicographic heuristic)",
    "; reasons=", paste(vapply(x$reasons, `[[`, "", "code"), collapse = ","),
    "; inspect structured explanation for marginal vectors, evidence and trace.")
}
