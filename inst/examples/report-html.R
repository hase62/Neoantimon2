# Sourceable helper; no project, predictor, network or output is started by sourcing.
# This is an example utility, not an additional exported Neoantimon2 API.
write_neo_report_html <- function(project, file, analyses = list(), candidate_id = NULL,
    privacy = "share", confirm_private = FALSE, max_candidates = 500L,
    max_bytes = 10 * 1024^2, report_options = list(), label = "Saved scientific results") {
  utils::getFromNamespace("html_report_write", "Neoantimon2")(
    project, file, analyses, candidate_id, privacy, confirm_private,
    max_candidates, max_bytes, report_options, label)
}
