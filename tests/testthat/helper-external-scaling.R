sys.source(test_path("..", "fixtures", "external-scaling", "workload.R"), envir = environment())
external_scaling_legacy <- function() {
  env <- new.env(parent = asNamespace("Neoantimon2"))
  sys.source(test_path("..", "fixtures", "external-scaling", "baseline.R"), envir = env)
  env
}
