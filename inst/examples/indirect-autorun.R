# Installed end-to-end SOFTWARE demonstration only. No real HLA proteins/models.
library(Neoantimon2)
root <- tempfile("neo-indirect-mock-demo-"); dir.create(root)
fixture <- system.file("extdata", "indirect-autorun", package = "Neoantimon2")
stopifnot(nzchar(fixture), nzchar(Sys.which("python3")))
exe <- file.path(root, "mock-netmhciipan")
stopifnot(file.copy(file.path(fixture, "mock-netmhciipan.py"), exe)); Sys.chmod(exe, "0755")
model <- file.path(root, "model"); dir.create(file.path(model, "data"), recursive = TRUE)
rows <- jsonlite::read_json(file.path(fixture, "synthetic-binding.json"))$rows
jsonlite::write_json(rows, file.path(model, "data/binding-map.json"), auto_unbox = TRUE, null = "null", digits = NA)
writeLines("good", file.path(model, "mode.txt"))
a <- jsonlite::read_json(file.path(fixture, "synthetic-foundation.json"))
ref <- file.path(fixture, "synthetic-reference.json")
run <- neo_indirect_autorun(a, "HLA-DRA*01:01~HLA-DRB1*01:01",
  list(path = ref, sha256 = digest::digest(file = ref, algo = "sha256", serialize = FALSE)),
  executable = exe, model = list(path = model, identity = "SYNTHETIC MOCK; not real NetMHCIIpan data"),
  authorization = list(status = "mock_only", terms = "Controlled software fixture; no real model execution"),
  output_dir = file.path(root, "run"))
stopifnot(run$status == "COMPLETED", !run$manifest$real_predictor_executed,
          run$result$contexts[[1L]]$complete_count == 5L)
cat("TESTED_WITH_MOCK: 5 synthetic core-presenter units\nREAL_PREDICTOR_E2E_NOT_TESTED\nOutput:", run$output_dir, "\n")
