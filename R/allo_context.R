# Explicit frozen identity fields: record/claim metadata must not enter this ID.
allo_context_identity_payload <- function(x) x[c("transplant_donor_subject_id",
  "transplant_recipient_subject_id", "source_subject_id", "self_subject_id",
  "presenter_subject_id", "responding_subject_id", "source_subject_role",
  "self_subject_role", "presenter_subject_role", "responding_subject_role",
  "transplant_type", "scenario_id", "direction", "tissue", "apc_compartment",
  "occasion_id")]

allo_contexts <- function(config) {
  subjects <- vapply(config$subjects, function(s) {
    allo_keys(s, "subject_id", label = "allo subject")
    allo_text(s$subject_id, "subject_id")
  }, "")
  if (length(subjects) < 2L || anyDuplicated(subjects)) stop("Distinct explicit allo subjects required")
  out <- lapply(config$contexts, function(x) {
    required <- c("transplant_donor_subject_id", "transplant_recipient_subject_id",
      "source_subject_id", "self_subject_id", "presenter_subject_id", "responding_subject_id",
      "source_subject_role", "self_subject_role", "presenter_subject_role", "responding_subject_role",
      "transplant_type", "scenario_id", "direction", "tissue", "apc_compartment", "occasion_id", "claim_level")
    allo_keys(x, required, "provenance", "allo context")
    semantic <- stats::setNames(lapply(required, function(k) allo_text(x[[k]], k)), required)
    d <- semantic$transplant_donor_subject_id; r <- semantic$transplant_recipient_subject_id
    if (d == r || !all(unlist(semantic[required[seq_len(6L)]]) %in% subjects)) stop("Invalid allo subject binding")
    scenarios <- list(sot_recipient_response = c("SOT", "donor_to_recipient", d, r),
      hsct_graft_versus_host_donor_apc = c("HSCT", "recipient_to_donor", r, d),
      hsct_host_versus_graft_recipient_apc = c("HSCT", "donor_to_recipient", d, r))
    expected <- scenarios[[semantic$scenario_id]]
    if (is.null(expected)) stop("Unsupported or ambiguous allo orientation; broader APC scenarios deferred")
    if (!identical(c(semantic$transplant_type, semantic$direction, semantic$source_subject_id,
      semantic$self_subject_id), expected) || semantic$presenter_subject_id != expected[[4L]] ||
      semantic$responding_subject_id != expected[[4L]]) stop("Contradictory allo orientation/role assignment")
    for (role in c("source", "self", "presenter", "responding")) {
      id <- semantic[[paste0(role, "_subject_id")]]
      if (semantic[[paste0(role, "_subject_role")]] != if (id == d) "donor" else "recipient")
        stop("Contradictory transplant role")
    }
    if (semantic$claim_level != "hypothetical_molecular_scenario") stop("Phase A establishes only hypothetical molecular scenarios")
    list(schema_version = "allo-context-1.0.0", allo_context_id = allo_id("context", allo_context_identity_payload(semantic)),
      semantics = semantic, frozen_projection = list(mode = "allo_hla", patient_id = r,
        source_subject_id = semantic$source_subject_id, self_subject_id = semantic$self_subject_id,
        presenter_subject_id = semantic$presenter_subject_id, direction = semantic$direction),
      provenance = if (is.null(x$provenance)) list() else allo_provenance(x$provenance),
      evaluability = allo_state("evaluated", scope = "explicit_role_assignment"))
  })
  if (!length(out) || anyDuplicated(vapply(out, function(x) x$semantics$scenario_id, "")))
    stop("Require nonempty contexts with unique scenario_id; orientations are separate analyses")
  allo_records(out)
}
allo_provenance <- function(x) {
  allo_keys(x, character(), c("path", "sha256", "retrieved_at", "source_url", "license",
    "cell_lineage", "cellular_origin_evidence_ref", "sampling_time", "display_label"), "allo provenance")
  for (k in names(x)) if (!is.null(x[[k]])) scalar_text(x[[k]], k)
  x
}
