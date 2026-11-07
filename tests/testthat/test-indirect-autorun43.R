autorun43_setup <- function(mode = "good") {
  a <- autorun_setup(mode)
  source <- system.file("extdata", "indirect-autorun", "mock-netmhciipan43k.py", package="Neoantimon2")
  stopifnot(file.copy(source, a$executable, overwrite=TRUE));Sys.chmod(a$executable,"0755")
  a$predictor_version <- "4.3k";a$model$identity <- "CONTROLLED 4.3k MOCK; no vendor assets"
  a
}
schema43_valid <- function(x, schema, definition = NULL) {
  skip_if_not_installed("jsonvalidate")
  f <- system.file("schema",paste0(schema,".json"),package="Neoantimon2")
  s <- jsonlite::read_json(f)
  if(!is.null(definition))s[["$ref"]] <- paste0("#/$defs/",definition)
  jsonvalidate::json_validate(jsonlite::toJSON(x,auto_unbox=TRUE,null="null",digits=NA),
    jsonlite::toJSON(s,auto_unbox=TRUE,null="null",digits=NA),engine="ajv",verbose=TRUE)
}

test_that("4.3k default runs BA explicitly and retains independent EL/core fields", {
  a <- autorun43_setup();a$predictor_version <- NULL;x <- autorun_call(a)
  expect_identical(x$status,"COMPLETED");expect_identical(x$manifest$provider_version,"4.3k")
  expect_identical(x$manifest$resolved_profile,ia43_profile())
  expect_identical(x$result$schema_version,"indirect-hla-result-netmhciipan43k-ba-1.0.0")
  expect_equal(autorun_ctx(x)$complete_count,5L)
  expect_false(x$manifest$real_predictor_executed)
  expect_identical(x$manifest$execution_settings$BA,TRUE)
  expect_identical(x$manifest$execution_settings$context_encoding,FALSE)
  expect_identical(x$manifest$threshold_policy$value,1000L)
  expect_false(identical(x$manifest$threshold_policy$version,indirect_threshold()$version))
  expect_identical(x$manifest$counting_policy,"unique_core_orientation_presenter_union")
  cmd <- jsonlite::read_json(file.path(x$output_dir,"raw/batch-00001.command.json"))
  expect_true("-BA" %in% unlist(cmd$argv));expect_false("-inv_all" %in% unlist(cmd$argv))
  expect_false("-tdir" %in% unlist(cmd$argv))
  expect_match(x$manifest$execution_settings$temporary_policy,"vendor TMPDIR default")
  m <- jsonlite::read_json(file.path(x$output_dir,"binding.json"))
  expect_identical(m$version,ia43_import_version());expect_true(schema43_valid(m,m$version))
  expect_true(schema43_valid(x$result,x$result$schema_version))
  expect_true(schema43_valid(x$explanation,x$explanation$schema_version))
  cfg <- jsonlite::read_json(file.path(x$output_dir,"indirect-config.json"))
  expect_true(schema43_valid(cfg,cfg$version,"config"))
  derived <- indirect_baseline(allo_foundation(allo_config(a$allo)),indirect_config(cfg,x$output_dir))
  expect_true(schema43_valid(derived,derived$schema_version,"baseline"))
  for(r in m$rows) {
    expect_identical(r$raw_output$columns$Score_EL,"0.123456")
    expect_identical(r$raw_output$columns$Core_Rel,"0.750")
    expect_identical(r$raw_output$columns$Inverted,"0")
    expect_equal(r$percentile_rank,1.5)
  }
  expect_true(file.exists(file.path(x$output_dir,"model-metrics.tsv")))
  expect_identical(autorun_call(a,resume=TRUE)$result,x$result)
  expect_length(x$plan$generation$occurrences,30L);expect_length(x$plan$query_links,35L)
})

test_that("documented inverted EL example uses reversed-ligand Of and cannot become BA", {
  f <- system.file("extdata","indirect-autorun","documented43k-el-example.txt",package="Neoantimon2")
  p <- "NPISEFVKWYKSHKL";presenter <- "HLA-DPA1*02:02~HLA-DPB1*19:01";vendor <- "HLA-DPA10202-DPB11901"
  r <- ia43_parse(readLines(f),p,presenter,vendor,require_ba=FALSE)[[1L]]
  expect_identical(r$status,"el_only");expect_null(r$affinity_nM)
  expect_identical(r$cores[[1L]]$core_sequence,"KYWKVFESI")
  expect_identical(r$cores[[1L]]$orientation,"reverse");expect_equal(r$cores[[1L]]$offset,4L)
  expect_identical(r$raw_output$forward_ligand_core_start0,2L)
  expect_identical(substring(p,3,11),"ISEFVKWYK")
  expect_identical(paste(rev(strsplit("ISEFVKWYK","")[[1]]),collapse=""),"KYWKVFESI")
  expect_identical(r$raw_output$columns$Score_EL,"0.600852")
  expect_identical(r$raw_output$columns[["%Rank_EL"]],"0.84")
  expect_error(ia43_parse(readLines(f),p,presenter,vendor),"EL-only is not BA")
})

test_that("inversion is retained and same-orientation self comparison is explicit", {
  a <- autorun43_setup("inverted");x <- autorun_call(a,inversion="all")
  cmd <- jsonlite::read_json(file.path(x$output_dir,"raw/batch-00001.command.json"))
  expect_true("-inv_all" %in% unlist(cmd$argv))
  m <- jsonlite::read_json(file.path(x$output_dir,"binding.json"))
  expect_true(all(vapply(m$rows,function(r)r$cores[[1]]$orientation=="reverse",TRUE)))
  expect_true(schema43_valid(x$explanation,x$explanation$schema_version))
  provider <- indirect_provider(m)
  # Independent explicit sequence constructions: identical displayed core, distinct orientation.
  pep <- c("ACDEFGHIKLMNPQR","AAAAAAKIHGFEDCA","TTTTTTKIHGFEDCA")
  make <- function(p,inv) {
    line <- paste("1 DRB1_0101",p,"0 ACDEFGHIK 1.0",inv,"Synthetic 0.1 20 NA 0.7 2 500")
    ia43_parse(c("# NetMHCIIpan version 4.3k",
      "Pos MHC Peptide Of Core Core_Rel Inverted Identity Score_EL %Rank_EL Exp_Bind Score_BA %Rank_BA Affinity(nM)",line),p,a$presenters[[1L]])[[1L]]
  }
  m$rows <- list(make(pep[[1]],0),make(pep[[2]],1),make(pep[[3]],1));provider <- indirect_provider(m)
  catalog <- list(complete=TRUE,reason_codes=list(),windows=list(list(peptide_sequence=pep[[2]],sequence_id="test",allele="HLA-A*01:01",start=1L,end=15L)),universe=list(self_universe_id="test:self",typing_evidence=list()),self_subject_id="test")
  search <- indirect_self_prediction(catalog,a$presenters[[1]],provider)
  expect_true(search$complete);expect_identical(search$hits[[1]]$orientation,"reverse")
  cmp <- function(inv)indirect_compare("test:query","core_pmhc","ACDEFGHIK",a$presenters[[1]],catalog,search,provider,orientation=inv)
  expect_identical(cmp("reverse")$status,"self_match")
  expect_identical(cmp("forward")$status,"supported_nonself_within_scope")
  keys <- vapply(provider$rows,function(r)indirect_core_key(provider,r$cores[[1]],r$presenter,"NCBITaxon:9606"),"")
  expect_length(unique(keys),2L)
  expect_true(all(grepl("^indirect43-core-presenter-id-v1:",keys)))
})

test_that("EL is never used as nM and malformed 4.3k results fail closed", {
  for(mode in c("el-only","badcolumns","bad-inverted","badcore","wrongallele","badnumber","missing-self","truncated","conflict","exit-zero-error","failure")) {
    a <- autorun43_setup(mode);expect_error(autorun_call(a),info=mode)
    m <- jsonlite::read_json(file.path(a$output_dir,"run.json"));expect_identical(m$status,"FAILED",info=mode)
    expect_false(file.exists(file.path(a$output_dir,"result.json")),info=mode)
  }
  a <- autorun43_setup();f <- file.path(a$model$path,"data/binding-map.json");rows <- jsonlite::read_json(f)
  for(i in seq_along(rows))rows[[i]]$affinity_nM <- 1000
  jsonlite::write_json(rows,f,auto_unbox=TRUE,null="null",digits=NA)
  x <- autorun_call(a);expect_identical(autorun_ctx(x)$complete_count,0L)
  expect_equal(autorun_ctx(x)$binding_policy_excluded_count,6L)
})

test_that("exact release and explicit legacy separation are enforced", {
  a <- autorun43_setup();writeLines(sub("4.3k","4.3j",readLines(a$executable),fixed=TRUE),a$executable)
  expect_error(autorun_call(a),"exact NetMHCIIpan 4.3k")
  a <- autorun43_setup();a$predictor_version <- "3.2";expect_error(autorun_call(a),"3.2")
  a <- autorun_setup();a$predictor_version <- "4.3k";expect_error(autorun_call(a),"4.3k")
  a <- autorun_setup();a$predictor_version <- "auto";x <- autorun_call(a)
  expect_identical(x$manifest$provider_version,"3.2");expect_equal(autorun_ctx(x)$complete_count,5L)
  a <- autorun43_setup();x <- autorun_call(a)
  m <- jsonlite::read_json(file.path(x$output_dir,"binding.json"));m$version <- "indirect-binding-import-1.0.0"
  expect_error(indirect_provider(m),"3.2")
})

test_that("internal distribution data symlink is pinned without following external links", {
  a <- autorun43_setup();root <- a$model$path;dir.create(file.path(root,"Darwin_arm64"))
  stopifnot(file.symlink("../data",file.path(root,"Darwin_arm64/data")))
  a$model$runtime_subdir <- "Darwin_arm64";p <- ia43_model(a$model)
  expect_identical(p$rdir,normalizePath(file.path(root,"Darwin_arm64")))
  expect_true(any(vapply(p$files,function(r)!is.null(r$link_target),TRUE)))
  outside <- tempfile();dir.create(outside);stopifnot(file.symlink(outside,file.path(root,"escape")))
  expect_error(ia43_model(a$model),"escapes pinned")
  a <- autorun43_setup();a$model$runtime_subdir <- "..";expect_error(ia43_model(a$model),"inside pinned")
})


test_that("official platform layout is chosen automatically without changing weights", {
  a <- autorun43_setup();platform <- paste(Sys.info()[["sysname"]],Sys.info()[["machine"]],sep="_")
  dir.create(file.path(a$model$path,platform,"bin"),recursive=TRUE)
  stopifnot(file.symlink("../data",file.path(a$model$path,platform,"data")))
  model <- ia43_model(a$model)
  expect_identical(model$rdir,normalizePath(file.path(a$model$path,platform)))
})
