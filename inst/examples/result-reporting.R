# Executable minimal SYNTHETIC reporting demonstration. Not biological validation.
library(Neoantimon2)
run_example <- function(output_root) {
  cfg <- neo_config(env = character(), paths = list(run_root = file.path(output_root, "runs")))
  p <- neo_run(neo_project(config = cfg, fixture_kind = "synthetic"))
  r <- neo_report_bundle(p, file.path(output_root, "report"))
  reopened <- neo_project(path = p$path)
  again <- neo_report_bundle(reopened)
  stopifnot(identical(r$tables$pmhc_candidates, again$tables$pmhc_candidates))
  invisible(r)
}
# A fresh temporary directory; no existing project or reference is modified.
example_reporting_result <- run_example(tempfile("neo-reporting-example-"))
