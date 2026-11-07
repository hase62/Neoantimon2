test_that("autorun reuses the frozen molecular result and retains every query path", {
  args <- autorun_setup(); x <- autorun_call(args); ctx <- autorun_ctx(x)
  golden <- indirect_test(); expected <- Neoantimon2:::indirect_public_result(golden$result, golden$foundation)$contexts[[1L]]
  expect_identical(x$status, "COMPLETED")
  for (field in c("complete_count", "conditional_count", "included_unit_count", "self_excluded_count", "unevaluable_count"))
    expect_equal(ctx[[field]], expected[[field]])
  expect_identical(x$manifest$execution_kind, "MOCK_EXECUTION"); expect_false(x$manifest$real_predictor_executed)
  expect_length(x$plan$queries, 7L)
  expect_length(x$plan$generation$occurrences, 30L)
  expect_identical(x$plan$self_window_count, 5L)
  expect_equal(sum(vapply(x$plan$query_links, function(r) r$role == "donor", TRUE)), 30L)
  expect_equal(sum(vapply(x$plan$query_links, function(r) r$role == "self", TRUE)), 5L)
  m <- jsonlite::read_json(file.path(x$output_dir, "binding.json"))
  expect_identical(m$provider$evidence_kind, "synthetic_fixture")
  expect_false(m$provenance$metadata$real_predictor_executed)
  expect_true(indirect_schema_valid(m)); expect_match(readLines(file.path(x$output_dir,"README.txt")), "REAL_PREDICTOR", all = FALSE)
  expect_true(file.exists(file.path(x$output_dir, "explanation.json")))
})

test_that("dry-run is offline and completed resume is read-only and hash checked", {
  args <- autorun_setup(); args$executable <- NULL; args$model <- NULL; args$authorization <- NULL
  x <- autorun_call(args, dry_run = TRUE)
  expect_identical(x$status, "DRY_RUN"); expect_null(x$result)
  expect_false(dir.exists(file.path(x$output_dir,"execution")))
  expect_identical(autorun_call(args, dry_run = TRUE, resume = TRUE)$plan, x$plan)
  expect_error(autorun_call(args, dry_run = TRUE), "directory exists")
  writeLines("tampering", file.path(x$output_dir,"queries.tsv"))
  expect_error(autorun_call(args, dry_run = TRUE, resume = TRUE), "checksum")
})

test_that("provider resume never executes again and rejects changed pins", {
  a <- autorun_setup();x <- autorun_call(a)
  before <- file.info(file.path(x$output_dir,"raw/batch-00001.stdout"))$mtime
  y <- autorun_call(a, resume = TRUE)
  expect_identical(y$result,x$result);expect_identical(file.info(file.path(x$output_dir,"raw/batch-00001.stdout"))$mtime,before)
  a$authorization$terms <- "changed";expect_error(autorun_call(a,resume=TRUE),"identity mismatch")
  a <- autorun_setup();a$reference$sha256 <- paste(rep("a",64),collapse="");expect_error(autorun_call(a),"integrity mismatch")
  a <- autorun_setup();a$model$sha256 <- paste(rep("a",64),collapse="");expect_error(autorun_call(a),"Model tree SHA")
})

test_that("all malformed/failed output modes fail closed and keep raw diagnostics", {
  modes <- c("badcolumns","badcore","bad-offset","wrongallele","badnumber","conflict","truncated","missing-self","failure")
  for (mode in modes) {
    a <- autorun_setup(mode)
    expect_error(autorun_call(a), info=mode)
    manifest <- jsonlite::read_json(file.path(a$output_dir,"run.json"))
    expect_identical(manifest$status,"FAILED",info=mode)
    expect_false(file.exists(file.path(a$output_dir,"result.json")),info=mode)
    expect_true(file.exists(file.path(a$output_dir,"raw/batch-00001.stdout")),info=mode)
    expect_error(autorun_call(a,resume=TRUE),"completed runs/plans only",info=mode)
  }
  a <- autorun_setup("timeout");expect_error(autorun_call(a,timeout=1L),"command failed")
})

test_that("genuine 3.2 identity/support is required and mock cannot be called real", {
  a <- autorun_setup();a$executable <- file.path(tempdir(),"no-such-netmhciipan");expect_error(autorun_call(a),"unavailable")
  a <- autorun_setup();a$model$path <- tempdir();expect_error(autorun_call(a),"model root|Model root")
  a <- autorun_setup("unsupported");expect_error(autorun_call(a),"Unsupported declared DR")
  a <- autorun_setup();writeLines(sub("version 3.2","version 4.3",readLines(a$executable),fixed=TRUE),a$executable)
  expect_error(autorun_call(a),"requires NetMHCIIpan 3.2")
  a <- autorun_setup();a$authorization$status <- "authorized_local_installation";expect_error(autorun_call(a),"Mock executable")
  a <- autorun_setup();a$presenters <- "HLA-DQA1*01:01~HLA-DQB1*02:01";expect_error(autorun_call(a),"Unsupported strict autorun presenter")
  a <- autorun_setup();a$presenters <- "HLA-DRA*01:01~HLA-DRB1*01:01:01";expect_error(autorun_call(a),"no truncation")
})

test_that("alternate documented BA layout and optional percentile remain explicit", {
  for (mode in c("pos-layout","no-rank")) {
    x <- autorun_call(autorun_setup(mode));expect_identical(x$status,"COMPLETED")
    m <- jsonlite::read_json(file.path(x$output_dir,"binding.json"))
    expect_true(all(vapply(m$rows,function(r) r$cores[[1L]]$offset==0L,TRUE)))
    if(mode=="no-rank")expect_true(all(vapply(m$rows,function(r)is.null(r$percentile_rank),TRUE)))
  }
})

test_that("cutoff is strictly below 1000 and core offset comes from provider", {
  a <- autorun_setup();p <- file.path(a$model$path,"data/binding-map.json");rows <- jsonlite::read_json(p)
  rows[[2L]]$affinity_nM <- 1000;rows[[3L]]$affinity_nM <- 999.999
  rows[[3L]]$cores[[1L]]$offset <- 2L;rows[[3L]]$cores[[1L]]$core_sequence <- substring(rows[[3L]]$peptide_sequence,3,11)
  jsonlite::write_json(rows,p,auto_unbox=TRUE,null="null",digits=NA)
  x <- autorun_call(a);expect_equal(autorun_ctx(x)$binding_policy_excluded_count,1L)
  evidence <- jsonlite::read_json(file.path(x$output_dir,"binding.json"))$rows
  r <- Filter(function(r)r$peptide_sequence==rows[[3L]]$peptide_sequence,evidence)[[1L]]
  expect_identical(r$cores[[1L]]$offset,2L);expect_identical(r$cores[[1L]]$orientation,"forward")
})

test_that("ordering and duplicate presenter declarations do not alter scientific identities", {
  a <- autorun_setup();x <- autorun_call(a)
  a$output_dir <- paste0(a$output_dir,"-second");a$allo$typing <- rev(a$allo$typing);a$presenters <- rep(a$presenters,2)
  y <- autorun_call(a)
  expect_identical(x$manifest$request_id,y$manifest$request_id)
  expect_identical(x$result,y$result)
  expect_identical(x$plan$query_links,y$plan$query_links)
})

test_that("incomplete self and source typing never gives a definitive zero", {
  a <- autorun_setup();a$allo$self_search[[1L]]$catalog_complete <- FALSE
  x <- autorun_call(a);expect_null(autorun_ctx(x)$complete_count)
  expect_identical(autorun_ctx(x)$completeness_status,"incomplete")
  a <- autorun_setup();a$allo$typing[[1L]]$complete <- FALSE
  x <- autorun_call(a);expect_null(autorun_ctx(x)$complete_count)
  a <- autorun_setup();a$reference$path <- tempfile();expect_error(autorun_call(a),"unavailable")
})

test_that("simple diploid typing runs end to end without invented alleles", {
  a <- autorun_setup();donor <- list(A=c("A*02:01","A*02:01"),B=c("B*07:02","B*07:02"),C=c("C*07:02","C*07:02"),DRB1=c("DRB1*04:01","DRB1*04:01"),DQB1=c("DQB1*03:01","DQB1*03:01"))
  recipient <- list(A=c("A*01:01","A*01:01"),B=c("B*08:01","B*08:01"),C=c("C*07:01","C*07:01"),DRB1=c("DRB1*01:01","DRB1*01:01"),DQB1=c("DQB1*02:01","DQB1*02:01"),DRA=c("DRA*01:01","DRA*01:01"))
  a$allo <- list(donor=donor,recipient=recipient);x <- autorun_call(a);expect_identical(x$status,"COMPLETED")
  a$allo$donor$A <- "A*02:01";a$output_dir <- paste0(a$output_dir,"bad");expect_error(autorun_call(a),"two explicit allele")
})

test_that("Origin IDs remain the frozen existing generation identities", {
  x <- autorun_call(); g <- indirect_test()$result$contexts[[1L]]$generation
  expect_setequal(vapply(x$plan$generation$occurrences,`[[`,"","indirect_peptide_occurrence_id"),
    vapply(g$occurrences,`[[`,"","indirect_peptide_occurrence_id"))
  expect_identical(x$plan$generation$generation_policy_id,g$generation_policy_id)
})

test_that("missing source allele and missing presenter scope stay insufficient", {
  a <- autorun_setup();a$allo$typing[[1L]]$alleles <- list("A*99:99")
  x <- autorun_call(a);expect_null(autorun_ctx(x)$complete_count)
  expect_identical(autorun_ctx(x)$completeness_status,"incomplete")
  a <- indirect_second_presenter();b <- autorun_setup();b$allo <- a$a;b$reference <- list(path=a$a$reference$path,sha256=digest::digest(file=a$a$reference$path,algo="sha256",serialize=FALSE))
  x <- autorun_call(b);expect_null(autorun_ctx(x)$complete_count)
  expect_true("presenter_scope_incomplete" %in% unlist(autorun_ctx(x)$reason_codes))
})

test_that("valid zero affinity is preserved and reverse/unknown core is rejected", {
  a <- autorun_setup();p <- file.path(a$model$path,"data/binding-map.json");rows <- jsonlite::read_json(p)
  rows[[2L]]$affinity_nM <- 0;jsonlite::write_json(rows,p,auto_unbox=TRUE,null="null",digits=NA)
  x <- autorun_call(a);m <- jsonlite::read_json(file.path(x$output_dir,"binding.json"));r <- Filter(function(z)z$peptide_sequence==rows[[2L]]$peptide_sequence,m$rows)[[1L]]
  expect_identical(r$affinity_nM,0L);expect_identical(r$status,"evaluated")
  lines <- c("# NetMHCIIpan version 3.2", "Seq Allele Peptide Identity Pos Core Core_Rel 1-log50k(aff) Affinity(nM)",
    "0 DRB1_0101 ACDEFGHIKLMNPQR Sequence 0 KIHGFEDCA 1 0.5 500")
  expect_error(Neoantimon2:::ia_parse(lines,"ACDEFGHIKLMNPQR",a$presenters[[1L]]),"forward core")
})

test_that("real-branch import formatting passes schema in an explicitly synthetic contract test", {
  # Exercises formatting only; no real provider, authorization or clinical evidence.
  a <- autorun_setup();x <- autorun_call(a);root <- x$output_dir
  r <- jsonlite::read_json(file.path(root,"binding.json"))$rows
  model <- Neoantimon2:::ia_model(a$model)
  tool <- list(sha256=digest::digest(file=a$executable,algo="sha256",serialize=FALSE),
    authorization=list(status="CONTRACT_TEST_ONLY",terms="Synthetic schema-test authorization, not real use"))
  probe <- list(version="3.2",presenter_map=stats::setNames("DRB1_0101",a$presenters[[1L]]))
  m <- Neoantimon2:::ia_import(r,file.path(root,"queries.tsv"),file.path(root,"execution/raw-output-index.json"),tool,model,probe,list(list(test="mock formatting only")),FALSE)
  expect_identical(m$provider$evidence_kind,"precomputed_predictor_output")
  expect_true(indirect_schema_valid(m))
  expect_identical(m$provenance$metadata$local_authorization_status,"CONTRACT_TEST_ONLY")
})


test_that("model hashing refuses symlinked content before traversal", {
  f <- autorun_setup()
  outside <- tempfile("model-outside-"); dir.create(outside)
  writeLines("unrelated", file.path(outside, "unrelated.txt"))
  link <- file.path(f$model$path, "data", "external-link")
  skip_if_not(file.symlink(outside, link), "symlinks unavailable")
  expect_error(ia_model(f$model), "without symlinks")
})
