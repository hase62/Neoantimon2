# T003 expectations use tiny metadata/files only; no biological pipeline.
test_that("configuration has explicit precedence and strict keys/types", {
  d <- tempfile(); dir.create(d)
  f <- file.path(d, "config.json")
  jsonlite::write_json(list(runtime = list(threads = 2)), f, auto_unbox = TRUE)
  env <- c(NEOANTIMON2_THREADS = "3")
  expect_equal(neo_config(f, env = env)$runtime$threads, 3)
  expect_equal(neo_config(f, env = env, runtime = list(threads = 4))$runtime$threads, 4)
  expect_error(neo_config(unknown = 1), "Unknown")
  expect_error(neo_config(runtime = list(threads = "4")), "threads")
  expect_error(neo_config(scientific = list(secret = "hidden")), "Unknown")
  expect_error(neo_config(paths = list(data_root = "relative")), "absolute")
  expect_equal(neo_config(paths = list(data_root = "~/tiny"))$paths$data_root,
               file.path(normalizePath("~"), "tiny"))
})

t003_config <- function(...) neo_config(env = character(), paths = list(
  run_root = file.path(tempdir(), "t003-runs"), data_root = file.path(tempdir(), "t003-data")), ...)
t003_samples <- function() data.frame(sample_id = c("s2", "s1"), patient_id = "p1")

test_that("project identity is stable and runtime independent", {
  p <- neo_project(samples = t003_samples(), config = t003_config())
  expect_s3_class(p, "neo_project")
  expect_null(p$run_id)
  expect_identical(p$analysis_id, neo_project(samples = t003_samples()[2:1, ],
    config = t003_config(runtime = list(threads = 7)))$analysis_id)
  expect_false(identical(p$analysis_id, neo_project(samples = t003_samples(),
    config = t003_config(scientific = list(genome_assembly = "GRCh37")))$analysis_id))
  expect_error(neo_project(mode = "allo_hla"), "unsupported")
})

test_that("runs coexist and resume verifies identity", {
  p <- neo_project(samples = t003_samples(), config = t003_config(), initialize = TRUE)
  q <- neo_project(samples = t003_samples(), config = t003_config(), initialize = TRUE)
  expect_false(identical(p$run_id, q$run_id))
  expect_true(all(dir.exists(file.path(p$path, c("manifest", "logs", "checkpoints",
    "intermediate", "cache", "results")))))
  r <- neo_project(path = p$path)
  expect_identical(r$analysis_id, p$analysis_id)
  expect_error(neo_project(path = p$path, samples = data.frame(sample_id = "other",
    patient_id = "p1")), "identity")
  m <- Neoantimon2:::read_manifest(p$path)
  expect_identical(m$run_manifest$schema_version, "1.0.0")
  expect_identical(m$run_manifest$id_schema_version, "1.0.0")
  expect_silent(Neoantimon2:::validate_table(list(m$run_manifest), "run_manifest"))
  expect_error(Neoantimon2:::update_run_status(p, "COMPLETED"), "transition")
  Neoantimon2:::update_run_status(p, "RUNNING")
  Neoantimon2:::update_run_status(p, "COMPLETED")
  expect_error(Neoantimon2:::update_run_status(p, "RUNNING"), "transition")
})

test_that("atomic writes preserve previous metadata on interruption", {
  f <- tempfile(fileext = ".json")
  Neoantimon2:::atomic_json(list(value = "old"), f)
  expect_error(Neoantimon2:::atomic_json(list(value = "new"), f,
    before_rename = function(...) stop("crash")), "crash")
  expect_equal(jsonlite::read_json(f)$value, "old")
  Neoantimon2:::atomic_json(list(value = "new"), f)
  expect_equal(jsonlite::read_json(f)$value, "new")
})

test_that("resources resolve by category with optional/required checks", {
  d <- tempfile(); dir.create(d)
  f <- file.path(d, "tiny.fa"); writeLines("tiny", f)
  cfg <- neo_config(env = character(), paths = list(data_root = d, tools_root = d),
    resources = list(reference.genome_fasta = list(path = "tiny.fa", required = TRUE),
      tools.netmhcpan = list(path = "missing", type = "executable")))
  p <- neo_project(config = cfg)
  expect_equal(p$resources[[1]]$resolved_path, normalizePath(f))
  expect_false(p$resources[[2]]$integrity$exists)
  cfg$resources[[2]]$required <- TRUE
  expect_error(neo_project(config = cfg), "Required resource missing")
  cfg$resources[[2]]$required <- FALSE
  cfg$resources[[1]]$type <- "directory"
  expect_error(neo_project(config = cfg), "type mismatch")
  y <- file.path(d, "config.yaml")
  writeLines(c("paths:", "  data_root: ./external", "runtime:", "  threads: 2"), y)
  expect_equal(neo_config(y, env = character())$paths$data_root, file.path(normalizePath(d), "external"))
  expect_error(neo_config(resources = list(reference.genome_fasta = list(pth = f))), "Unknown")
})

test_that("checksum modes are explicit, declarations do not imply verification", {
  f <- tempfile(); writeLines("tiny", f)
  supplied <- list(algorithm = "sha256", scope = "file",
    value = digest::digest(file = f, algo = "sha256", serialize = FALSE))
  stat <- Neoantimon2:::integrity_record(f, "metadata", supplied)
  expect_false(stat$verified)
  expect_equal(stat$observed$scope, "file-metadata")
  expect_true(Neoantimon2:::integrity_record(f, "sha256", supplied)$verified)
  expect_equal(Neoantimon2:::integrity_record(f, "none")$observed$algorithm, "none")
  supplied$value <- paste(rep("0", 64), collapse = "")
  expect_error(Neoantimon2:::integrity_record(f, "sha256", supplied), "mismatch")
  expect_error(Neoantimon2:::integrity_record(tempdir(), "sha256"), "Directory")
  expect_equal(Neoantimon2:::integrity_record(tempdir(), "metadata")$observed$scope, "directory-entry-only")
})

test_that("input integrity and resume do not silently accept changed files", {
  f <- tempfile(); writeLines("one", f)
  inputs <- list(list(input_id = "bam1", input_role = "tumor_dna_alignment", sample_id = "s1", path = f))
  p <- neo_project(samples = t003_samples(), inputs = inputs,
    config = t003_config(), initialize = TRUE)
  expect_identical(neo_project(path = p$path)$analysis_id, p$analysis_id)
  writeLines("different size", f)
  expect_error(neo_project(path = p$path), "identity")
  q <- neo_project(samples = t003_samples(), inputs = inputs, config = t003_config())
  expect_false(identical(p$analysis_id, q$analysis_id))
})

test_that("secret environment is not printed or persisted", {
  env <- c(AWS_SECRET_ACCESS_KEY = "do-not-leak", NEOANTIMON2_TOKEN = "do-not-leak")
  cfg <- neo_config(env = env, paths = t003_config()$paths)
  p <- neo_project(samples = t003_samples(), config = cfg, initialize = TRUE)
  text <- paste(readLines(Neoantimon2:::manifest_file(p$path)), collapse = "")
  expect_false(grepl("do-not-leak|AWS_SECRET|TOKEN", text))
  expect_false(any(grepl("do-not-leak", capture.output(print(p)))))
  expect_false(any(grepl("do-not-leak", capture.output(print(cfg)))))
  expect_error(neo_config(resources = list(tools.netmhcpan = list(path = NULL, token = "secret"))), "Unknown")
  expect_error(neo_config(resources = list(tools.netmhcpan = list(source = "https://host/?token=secret"))), "Credential")
})

test_that("checkpoints have safe state and step-scoped fingerprints", {
  p <- neo_project(samples = t003_samples(), config = t003_config(), initialize = TRUE)
  fp <- Neoantimon2:::checkpoint_fingerprints(list(input = "a"), list(reference = "v1"), list())
  other <- Neoantimon2:::checkpoint_fingerprints(list(input = "b"), list(reference = "v1"), list())
  Neoantimon2:::checkpoint_update(p, "reference", fingerprints = fp)
  expect_false(Neoantimon2:::checkpoint_reusable(p, "reference", fp))
  expect_error(Neoantimon2:::checkpoint_update(p, "reference", "COMPLETED", fp), "transition")
  Neoantimon2:::checkpoint_update(p, "reference", "RUNNING", fp)
  expect_false(Neoantimon2:::checkpoint_reusable(neo_project(path = p$path), "reference", fp))
  Neoantimon2:::checkpoint_update(p, "reference", "COMPLETED", fp)
  expect_true(Neoantimon2:::checkpoint_reusable(p, "reference", fp))
  expect_error(Neoantimon2:::checkpoint_reusable(p, "reference", other), "fingerprint mismatch")
  Neoantimon2:::checkpoint_update(p, "reference", "INVALIDATED", fp)
  expect_false(Neoantimon2:::checkpoint_reusable(p, "reference", fp))
  Neoantimon2:::checkpoint_update(p, "reference", "PENDING", other)
  Neoantimon2:::checkpoint_update(p, "reference", "RUNNING", other)
  Neoantimon2:::checkpoint_update(p, "reference", "FAILED", other, error = list(code = "INPUT_CHANGED"))
  expect_equal(Neoantimon2:::read_manifest(p$path)$checkpoints[[1]]$error$code, "INPUT_CHANGED")
  # Grade changes are absent from reference/event step slices.
  a <- list(reference = "v1", grade = "policy1"); b <- a; b$grade <- "policy2"
  expect_identical(Neoantimon2:::checkpoint_fingerprints(list(), a["reference"], list()),
    Neoantimon2:::checkpoint_fingerprints(list(), b["reference"], list()))
  expect_false(identical(Neoantimon2:::checkpoint_fingerprints(list(), a["grade"], list()),
    Neoantimon2:::checkpoint_fingerprints(list(), b["grade"], list())))
})

test_that("checkpoint success requires intact outputs and never follows partial writes", {
  p <- neo_project(config = t003_config(), initialize = TRUE)
  fp <- Neoantimon2:::checkpoint_fingerprints(list(), list(), list())
  Neoantimon2:::checkpoint_update(p, "step", fingerprints = fp)
  Neoantimon2:::checkpoint_update(p, "step", "RUNNING", fp)
  output <- list(list(path = "results/tiny.txt", checksum = list(algorithm = "sha256", scope = "file",
    value = paste(rep("0", 64), collapse = ""))))
  expect_error(Neoantimon2:::checkpoint_update(p, "step", "COMPLETED", fp, output), "missing")
  expect_equal(Neoantimon2:::read_manifest(p$path)$checkpoints[[1]]$status, "RUNNING")
  f <- file.path(p$path, "results/tiny.txt"); writeLines("done", f)
  output[[1]]$checksum$value <- digest::digest(file = f, algo = "sha256", serialize = FALSE)
  Neoantimon2:::checkpoint_update(p, "step", "COMPLETED", fp, output)
  expect_true(Neoantimon2:::checkpoint_reusable(p, "step", fp))
  writeLines("changed", f)
  expect_error(Neoantimon2:::checkpoint_reusable(p, "step", fp), "checksum mismatch")
})

test_that("run locks, manifest corruption and unsafe locations reject", {
  p <- neo_project(config = t003_config(), initialize = TRUE)
  lock <- file.path(p$path, ".metadata-lock"); dir.create(lock)
  expect_error(neo_project(path = p$path), "locked")
  unlink(lock, recursive = TRUE)
  m <- Neoantimon2:::read_manifest(p$path)
  m$run_manifest$schema_version <- "99.0.0"
  Neoantimon2:::atomic_json(m, Neoantimon2:::manifest_file(p$path))
  expect_error(neo_project(path = p$path), "enum")
  expect_error(neo_project(path = tempfile(), config = t003_config()), "inside run_root")
  expect_error(neo_project(config = neo_config(env = character(), paths = list(run_root = tempdir(),
    data_root = tempdir())), initialize = TRUE), "overlaps")
})

test_that("symlinks are canonical and external substitutions are rejected", {
  d <- tempfile(); dir.create(d)
  f <- file.path(d, "input"); writeLines("tiny", f)
  link <- file.path(d, "link")
  expect_true(file.symlink(f, link))
  expect_equal(Neoantimon2:::resolve_path(link), normalizePath(f))
  expect_equal(Neoantimon2:::resolve_path(file.path(d, "new", "..", "leaf")), file.path(normalizePath(d), "leaf"))
  p <- neo_project(samples = t003_samples(), config = t003_config(), initialize = TRUE,
    inputs = list(list(input_id = "one", input_role = "tumor_dna_alignment", sample_id = "s1", path = link)))
  unlink(f)
  expect_error(neo_project(path = p$path), "Broken symlink")
})

test_that("moving a run preserves run-relative output locators", {
  cfg <- t003_config()
  p <- neo_project(config = cfg, initialize = TRUE)
  target <- file.path(cfg$paths$run_root, paste0("moved-", basename(p$path)))
  expect_true(file.rename(p$path, target))
  q <- neo_project(path = target)
  expect_identical(q$analysis_id, p$analysis_id)
  expect_identical(q$run_id, p$run_id)
  expect_equal(q$path, target)
})

test_that("reopening unrelated directories performs no writes", {
  d <- tempfile(); dir.create(d)
  expect_error(neo_project(path = d), "manifest missing")
  expect_length(list.files(d, all.files = TRUE, no.. = TRUE), 0)
})

test_that("runtime/output location and resource ordering do not affect scientific identity", {
  a <- t003_config(resources = list(reference.hla = list(version = "v1"),
    reference.mane = list(version = "v2")))
  b <- a; b$paths$run_root <- file.path(tempdir(), "other-runs")
  b$runtime$verbosity <- "debug"; b$resources <- rev(b$resources)
  expect_identical(neo_project(config = a)$analysis_id, neo_project(config = b)$analysis_id)
  a$resources[[1]]$version <- "changed"
  expect_false(identical(neo_project(config = a)$analysis_id, neo_project(config = b)$analysis_id))
})

test_that("abrupt process exit before rename preserves committed JSON", {
  f <- tempfile(fileext = ".json")
  Neoantimon2:::atomic_json(list(value = "old"), f)
  script <- tempfile(fileext = ".R")
  code <- c("atomic_json <-", deparse(Neoantimon2:::atomic_json),
    paste0("atomic_json(list(value = 'new'), ", encodeString(f, quote = '"'),
      ", before_rename = function(...) quit(save = 'no', status = 77, runLast = FALSE))"))
  writeLines(code, script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)),
    stdout = tempfile(), stderr = tempfile())
  expect_equal(status, 77)
  expect_equal(jsonlite::read_json(f)$value, "old")
})

test_that("a relocated run still rejects retargeted declared-checksum inputs", {
  d <- tempfile(); dir.create(d)
  a <- file.path(d, "a"); b <- file.path(d, "b"); writeLines("same", a); writeLines("same", b)
  link <- file.path(d, "link"); expect_true(file.symlink(a, link))
  checksum <- list(algorithm = "sha256", scope = "file",
    value = digest::digest(file = a, algo = "sha256", serialize = FALSE))
  p <- neo_project(config = t003_config(), samples = t003_samples(), initialize = TRUE,
    inputs = list(list(input_id = "x", input_role = "tumor_dna_alignment", sample_id = "s1", path = link, checksum = checksum)))
  unlink(link); expect_true(file.symlink(b, link))
  expect_error(neo_project(path = p$path), "path/integrity")
})

test_that("a manifest directory symlink cannot redirect metadata updates", {
  p <- neo_project(config = t003_config(), initialize = TRUE)
  outside <- tempfile(); dir.create(outside)
  original <- file.path(p$path, "manifest")
  moved <- file.path(p$path, "saved-manifest")
  expect_true(file.rename(original, moved))
  expect_true(file.symlink(outside, original))
  expect_error(Neoantimon2:::update_run_status(p, "RUNNING"), "Unsafe metadata")
  expect_length(list.files(outside, all.files = TRUE, no.. = TRUE), 0)
})

test_that("HLA genotype declarations bind once per subject", {
  hla <- list(p1 = c("A*02:01", "B*07:02"))
  p <- neo_project(samples = t003_samples(), hla = hla, config = t003_config(), initialize = TRUE)
  expect_named(p$hla, "p1")
  expect_length(p$hla, 1)
  expect_identical(neo_project(path = p$path)$hla, p$hla)
  expect_identical(p$analysis_id, neo_project(samples = t003_samples()[2:1, ],
    hla = hla, config = t003_config())$analysis_id)
  expect_error(neo_project(samples = t003_samples(), hla = list(s1 = "A*02:01"),
    config = t003_config()), "hla subject")
  expect_false(identical(p$analysis_id, neo_project(samples = t003_samples(),
    hla = list(p1 = "A*03:01"), config = t003_config())$analysis_id))
})

test_that("input roles define scientific identity and preserve multiplicity", {
  f <- tempfile(); writeLines("one", f)
  g <- tempfile(); writeLines("two", g)
  inputs <- list(list(input_id = "a", input_role = "somatic_variants",
    sample_id = "s1", path = f, version = "v1"))
  project <- function(x) neo_project(samples = t003_samples(), inputs = x,
    config = t003_config(checksum = "sha256"))
  p <- project(inputs)
  renamed <- inputs; renamed[[1]]$input_id <- "renamed"
  expect_identical(p$analysis_id, project(renamed)$analysis_id)
  changed <- inputs; changed[[1]]$input_role <- "germline_variants"
  expect_false(identical(p$analysis_id, project(changed)$analysis_id))
  changed <- inputs; changed[[1]]$version <- "v2"
  expect_false(identical(p$analysis_id, project(changed)$analysis_id))
  multiple <- c(inputs, renamed); multiple[[2]]$path <- g
  q <- project(multiple)
  expect_length(q$inputs, 2)
  expect_identical(q$analysis_id, project(rev(multiple))$analysis_id)
  expect_false(identical(p$analysis_id, project(c(inputs, renamed))$analysis_id))
  expect_error(project(c(inputs, inputs)), "Duplicate input_id")
  for (role in list(NULL, "", "a/b", "../x", "has space", "Somatic")) {
    invalid <- inputs; invalid[[1]]$input_role <- role
    expect_error(project(invalid), "input_role")
  }
  writeLines("new content", f)
  expect_false(identical(p$analysis_id, project(inputs)$analysis_id))
  run <- Neoantimon2:::initialize_run(q)
  manifest <- Neoantimon2:::read_manifest(run$path)
  expect_setequal(vapply(manifest$input_records, `[[`, "", "input_id"), c("a", "renamed"))
})

test_that("run overlap rejects either containment direction", {
  base <- tempfile(); dir.create(base); base <- normalizePath(base)
  run <- file.path(base, "run")
  for (key in c("data_root", "repository")) {
    for (forbidden in c(file.path(run, "nested"), base, run)) {
      paths <- list(run_root = base, data_root = file.path(base, "data"))
      paths[[key]] <- forbidden
      cfg <- neo_config(env = character(), paths = paths)
      expect_error(Neoantimon2:::check_run_location(run, cfg), "overlaps")
    }
  }
  cfg <- neo_config(env = character(), paths = list(run_root = base,
    data_root = paste0(run, "-data"), repository = paste0(run, "-repo")))
  expect_identical(Neoantimon2:::check_run_location(run, cfg), run)
  expect_true(Neoantimon2:::within_path("/a/b/", "/a/"))
  expect_true(Neoantimon2:::within_path("/a", "/"))
  expect_true(Neoantimon2:::within_path("/", "/"))
  expect_false(Neoantimon2:::within_path("/ab", "/a/"))
})

test_that("reference overrides use the config resource validator", {
  f <- tempfile(); writeLines("tiny", f)
  spec <- list(path = f, required = TRUE, version = "v1")
  cfg <- t003_config(resources = list(reference.genome = list(path = "missing")))
  expect_identical(neo_project(config = cfg, reference = list())$analysis_id,
    neo_project(config = cfg)$analysis_id)
  p <- neo_project(config = cfg, reference = list(genome = spec))
  expect_identical(p$config$resources$reference.genome$path, f)
  expect_identical(unclass(p$config), unclass(do.call(neo_config,
    c(unclass(p$config), list(env = character())))))
  for (key in c("bad/name", "bad name", "../escape", "")) {
    expect_error(neo_project(config = cfg, reference = setNames(list(spec), key)))
    expect_error(neo_config(env = character(),
      resources = setNames(list(spec), paste0("reference.", key))))
  }
  expect_error(neo_project(config = cfg, reference = list(genome = list(type = "invalid"))),
    "resource type")
  expect_error(neo_project(config = cfg, reference = list(genome = list(checksum = list()))),
    "checksum")
  expect_error(neo_project(config = cfg, reference = list(genome = list(type = "directory", path = f))),
    "type mismatch")
})

test_that("biological manifest projects only existing declared resources", {
  f <- tempfile(); writeLines("tiny", f)
  checksum <- list(algorithm = "sha256", scope = "file",
    value = digest::digest(file = f, algo = "sha256", serialize = FALSE))
  spec <- list(path = f, version = "v1", checksum = checksum)
  missing <- spec; missing$path <- paste0(f, "-missing")
  cfg <- t003_config(resources = list(reference.present = spec, reference.absent = missing))
  p <- neo_project(config = cfg, initialize = TRUE)
  m <- Neoantimon2:::read_manifest(p$path)
  expect_length(m$resource_records, 2)
  expect_false(m$resource_records[[1]]$integrity$exists)
  expect_identical(m$resource_records[[1]]$checksum, checksum)
  expect_length(m$run_manifest$resources, 1)
  expect_identical(m$run_manifest$resources[[1]]$resource_id, "reference.present")
  cfg$resources$reference.absent$required <- TRUE
  expect_error(neo_project(config = cfg), "Required resource missing")
})

test_that("checkpoint output namespaces exclude infrastructure and cache", {
  p <- neo_project(config = t003_config(), initialize = TRUE)
  fp <- Neoantimon2:::checkpoint_fingerprints(list(), list(), list())
  checksum <- list(algorithm = "sha256", scope = "file", value = strrep("0", 64))
  for (prefix in c("results", "intermediate", "checkpoints")) {
    output <- list(list(path = paste0(prefix, "/x"), checksum = checksum))
    expect_silent(Neoantimon2:::checkpoint_update(p, prefix, fingerprints = fp,
      output_records = output))
  }
  for (path in c("manifest/run.json", "logs/x", "cache/x", ".metadata-lock",
    "x", "../x", "/results/x", "results/../manifest/run.json", "results/./x",
    "results//x", "results/", "results/\\x")) {
    expect_error(Neoantimon2:::checkpoint_update(p, "bad", fingerprints = fp,
      output_records = list(list(path = path, checksum = checksum))),
      "run-relative|namespace")
  }
  # A valid lexical prefix must not permit symlinks into manifest/.
  link <- file.path(p$path, "results", "redirect")
  expect_true(file.symlink(file.path(p$path, "manifest"), link))
  cp <- list(output_records = list(list(path = "results/redirect/run.json", checksum = checksum)))
  expect_error(Neoantimon2:::verify_checkpoint_outputs(cp, p$path), "unsafe")
})

test_that("the distributed project YAML loads and duplicate keys never disappear", {
  source <- test_path("..", "..", "config", "project.example.yaml")
  fixture <- test_path("..", "fixtures", "config", "project.example.yaml")
  if (file.exists(source)) expect_identical(readLines(source), readLines(fixture))
  # R CMD build excludes the repository config directory; ship its exact fixture.
  example <- normalizePath(if (file.exists(source)) source else fixture, mustWork = TRUE)
  cfg <- neo_config(example, env = character())
  expect_s3_class(cfg, "neo_config")
  expect_length(cfg$resources, 6)
  expect_identical(cfg$resources$reference.genome_fasta$type, "file")
  cases <- list(
    json = c('{"checksum":"none","checksum":"metadata"}',
      '{"runtime":{"threads":1,"threads":2}}',
      '{"resources":{"reference.hla":{"type":"file","type":"model"}}}',
      '{"resources":{"reference.hla":{"checksum":{"algorithm":"sha256","algorithm":"sha256"}}}}'),
    yaml = c("checksum: none\nchecksum: metadata",
      "runtime:\n  threads: 1\n  threads: 2",
      "resources:\n  reference.hla:\n    type: file\n    type: model",
      "resources:\n  reference.hla:\n    checksum:\n      algorithm: sha256\n      algorithm: sha256"))
  for (format in names(cases)) for (text in cases[[format]]) {
    f <- tempfile(fileext = paste0(".", format)); writeLines(text, f)
    # Reject before explicit arguments could hide the duplicate subtree.
    expect_error(neo_config(f, env = character(), runtime = list(threads = 3),
      resources = list()), "[Dd]uplicate|unique named")
  }
})
