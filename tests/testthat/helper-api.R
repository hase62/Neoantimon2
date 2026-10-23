# Approved package API, grouped by integration milestone. This is test data,
# independent of NAMESPACE, so accidental exports cannot update the oracle.
public_api_by_milestone <- list(
  T003 = c("neo_project", "neo_config"),
  T004 = c("neo_run", "neo_results", "neo_events", "neo_origins",
           "neo_multisample", "neo_explain"),
  T058_T063 = c("neo_export_bundle", "neo_analyze_external",
               "neo_external_results", "neo_external_explain"),
  phase3_setup = c("neo_download_resources", "neo_check_resources",
                   "neo_check_predictors", "neo_setup_predictor"),
  later = c("neo_regrade", "neo_import_prediction", "neo_doctor")
)
approved_public_api <- unlist(public_api_by_milestone, use.names = FALSE)
t004_public_api <- unlist(public_api_by_milestone[c("T003", "T004")],
                          use.names = FALSE)
current_public_api <- unlist(public_api_by_milestone[c("T003", "T004", "T058_T063", "phase3_setup")],
                            use.names = FALSE)
