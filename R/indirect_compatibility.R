# Descriptive, patient-independent transcription of the frozen local matrix.
indirect_compatibility <- function() {
  path <- system.file("extdata", "indirect-hla-frozen-audit-1.0.0.json", package = "Neoantimon2")
  x <- jsonlite::read_json(path)
  indirect_audit_projection(x)
}
indirect_audit_projection <- function(x) {
  x$targets <- indirect_report_order(x$targets, "target_label")
  for (i in seq_along(x$targets)) {
    t <- x$targets[[i]]
    t$dimensions <- indirect_report_order(t$dimensions, "dimension")
    # Non-biological audit record: normalized target/configuration and every
    # dimension, including the independent reference_release axis.
    # It deliberately has no patient/context/provider-artifact input.
    t$audit_record_id <- paste0("audit-record-sha256:", digest::digest(canonical_bytes(
      t[c("target_label", "application_version", "configuration", "matrix_profile", "dimensions")]),
      algo = "sha256", serialize = FALSE))
    x$targets[[i]] <- t
  }
  x
}
