test_that("registry separates implemented adapters from installed predictors", {
  registry <- Neoantimon2:::predictor_registry()
  expect_setequal(names(registry),c("netmhcpan","netmhciipan","mhcflurry","fixture"))
  expect_true(all(vapply(registry,`[[`,TRUE,"implemented")))
  expect_false("U" %in% strsplit(registry$netmhcpan$alphabet,"")[[1L]])
  expect_true("U" %in% strsplit(registry$fixture$alphabet,"")[[1L]])
  expect_equal(Neoantimon2:::predictor_version_from_text("# NetMHCpan version 4.1b","netmhcpan"),"4.1b")
  expect_equal(Neoantimon2:::predictor_version_from_text("# NetMHCpan version 4.2estatic","netmhcpan"),"4.2estatic")
  expect_equal(Neoantimon2:::predictor_version_from_text("# NetMHCIIpan version 4.3k","netmhciipan"),"4.3k")
  expect_equal(Neoantimon2:::predictor_version_from_text("mhcflurry 2.0.1","mhcflurry"),"2.0.1")
  expect_equal(Neoantimon2:::predictor_version_from_text("mhcflurry 2.3.0rc20","mhcflurry"),"2.3.0rc20")
  expect_null(Neoantimon2:::predictor_version_from_text("unknown","netmhcpan"))
  fake <- tempfile("netmhcpan-version-")
  for (version in c("4.1b","4.2d","4.2e","4.2estatic","4.2.1","4.3","4.20")) {
    writeLines(c("#!/bin/sh",paste("echo '# NetMHCpan version",version,"'")),fake)
    Sys.chmod(fake,"0755")
    cap <- Neoantimon2:::predictor_capability("netmhcpan",executable=fake)
    expect_identical(cap$predictor_version,version)
    if (version %in% c("4.3","4.20")) {
      expect_false(cap$available)
      expect_identical(cap$availability_status,"unsupported_version")
    } else expect_true(cap$available)
  }
})

test_that("capability failures retain biological candidates and U literally", {
  cap <- predictor_test_capability("netmhcpan")
  a <- predictor_test_candidate("a","ACDEFGHIK")
  u <- predictor_test_candidate("u","ACDEFGHUK")
  long <- predictor_test_candidate("long",paste(rep("A",20L),collapse=""))
  unsupported <- predictor_test_candidate("x","ACDEFGHIK","HLA-E*01:01")
  class_ii <- predictor_test_candidate("ii","ACDEFGHIKLMNPQR",
    "HLA-DQA1*05:01~HLA-DQB1*02:01")
  p <- Neoantimon2:::predictor_prepare(list(u,long,unsupported,a),"netmhcpan","mhc1_binding",cap)
  expect_equal(length(p$queries),1L)
  reasons <- vapply(p$statuses,`[[`,"","reason")
  expect_setequal(reasons,c("unsupported_residue","unsupported_length","unsupported_allele_format"))
  expect_identical(p$statuses[[which(reasons=="unsupported_residue")]]$residue,"U")
  expect_identical(u$ligand_sequence,"ACDEFGHUK")
  expect_identical(Neoantimon2:::predictor_prepare(list(class_ii),"netmhcpan",
    "mhc1_binding",cap)$statuses[[1L]]$reason,"unsupported_presenter")
  u_only <- Neoantimon2:::predictor_prepare(list(u),"netmhcpan","mhc1_binding",cap)
  u_run <- Neoantimon2:::predictor_run(u_only,tempfile("unsupported-u-"))
  expect_length(u_run$evidence,0L)
  expect_identical(u_run$statuses[[1L]]$reason,"unsupported_residue")
  cap$available <- FALSE; cap$availability_status <- "predictor_unavailable"
  q <- Neoantimon2:::predictor_prepare(list(a),"netmhcpan","mhc1_binding",cap)
  expect_identical(q$statuses[[1L]]$reason,"predictor_unavailable")
  expect_length(q$queries,0L)
})

test_that("canonical presenters map only at adapter boundary", {
  i <- predictor_test_candidate("i")
  dq <- predictor_test_candidate("dq","ACDEFGHIKLMNPQR",
    "HLA-DQA1*05:01~HLA-DQB1*02:01")
  dr <- predictor_test_candidate("dr","ACDEFGHIKLMNPQR",
    "HLA-DRA*01:01~HLA-DRB1*04:01")
  dp <- predictor_test_candidate("dp","ACDEFGHIKLMNPQR",
    "HLA-DPA1*01:03~HLA-DPB1*04:01")
  expect_identical(Neoantimon2:::predictor_allele(i,"netmhcpan"),"HLA-A02:01")
  expect_identical(Neoantimon2:::predictor_allele(i,"mhcflurry"),"HLA-A*02:01")
  expect_identical(Neoantimon2:::predictor_allele(dq,"netmhciipan"),"HLA-DQA10501-DQB10201")
  expect_identical(Neoantimon2:::predictor_allele(dr,"netmhciipan"),"DRB1_0401")
  expect_identical(Neoantimon2:::predictor_allele(dp,"netmhciipan"),"HLA-DPA10103-DPB10401")
  expect_identical(dq$mhc_key,"HLA-DQA1*05:01~HLA-DQB1*02:01")
  ci <- Neoantimon2:::predictor_prepare(list(i),"netmhcpan","mhc1_binding",
    predictor_test_capability("netmhcpan"))
  cii <- Neoantimon2:::predictor_prepare(list(dq),"netmhciipan","mhc2_binding",
    predictor_test_capability("netmhciipan"))
  cmd_i <- Neoantimon2:::predictor_command(ci,"HLA-A02:01","input.txt")
  cmd_ii <- Neoantimon2:::predictor_command(cii,"HLA-DQA10501-DQB10201","input.txt")
  expect_true(all(c("-p","-a","-BA") %in% cmd_i$args))
  expect_true(all(c("-f","-inptype","1","-a","-BA") %in% cmd_ii$args))
  expect_match(Neoantimon2:::predictor_input(cii,"HLA-DQA10501-DQB10201"),
    "ACDEFGHIKLMNPQR")
})

test_that("batching deduplicates external queries but preserves provenance", {
  cap <- predictor_test_capability("netmhcpan")
  a <- predictor_test_candidate("a")
  b <- predictor_test_candidate("b")
  c <- predictor_test_candidate("c","ACDEFGHIK","HLA-B*07:02")
  first <- Neoantimon2:::predictor_prepare(list(c,b,a),"netmhcpan","mhc1_binding",cap)
  second <- Neoantimon2:::predictor_prepare(list(a,c,b),"netmhcpan","mhc1_binding",cap)
  expect_equal(first$queries,second$queries)
  expect_equal(length(first$queries),2L)
  expect_setequal(unlist(first$queries[[1L]]$candidate_ids),c("a","b"))
  expect_equal(Neoantimon2:::predictor_run_identity(first),Neoantimon2:::predictor_run_identity(second))
  expect_false(identical(Neoantimon2:::predictor_run_identity(first,list(timeout_seconds=10)),
    Neoantimon2:::predictor_run_identity(first,list(timeout_seconds=20))))
  cap$predictor_version <- "4.1b"
  third <- Neoantimon2:::predictor_prepare(list(a,b,c),"netmhcpan","mhc1_binding",cap)
  expect_false(identical(Neoantimon2:::predictor_run_identity(first),Neoantimon2:::predictor_run_identity(third)))
  cap$predictor_version <- "4.1"
  cap$model_version <- "model-B"
  model_changed <- Neoantimon2:::predictor_prepare(list(a,b,c),"netmhcpan","mhc1_binding",cap)
  expect_false(identical(Neoantimon2:::predictor_run_identity(first),
    Neoantimon2:::predictor_run_identity(model_changed)))
  expect_match(Neoantimon2:::predictor_input(first,"HLA-A02:01"),"ACDEFGHIK")
})

test_that("NetMHCpan 4.2 MHC-header fixture normalizes standard EL and BA", {
  cap <- predictor_test_capability("netmhcpan")
  cap$predictor_version <- "4.2d"
  candidates <- list(predictor_test_candidate("a"),predictor_test_candidate("b"),
    predictor_test_candidate("c","LMNPQRSTV"),
    predictor_test_candidate("d","LMNPQRSTV","HLA-B*07:02"))
  p <- Neoantimon2:::predictor_prepare(rev(candidates),"netmhcpan","mhc1_presentation",cap)
  lines <- readLines(testthat::test_path("..","fixtures","predictor","netmhcpan-4.2.txt"))
  header <- lines[grepl("^Pos ",lines)]
  rows <- lines[grepl("^[0-9]",lines)]
  raw <- list("HLA-A02:01"=paste(c(header,rev(rows[grepl("HLA-A\\*02:01",rows)])),collapse="\n"),
    "HLA-B07:02"=paste(c(header,rows[grepl("HLA-B\\*07:02",rows)]),collapse="\n"))
  result <- Neoantimon2:::predictor_normalize(p,raw,"run-42",
    list("HLA-A02:01"="a42.txt","HLA-B07:02"="b42.txt"))
  expect_length(result$statuses,0L)
  expect_length(result$evidence,20L)
  expect_length(result$core_evidence,0L)
  expect_setequal(unique(vapply(result$evidence,`[[`,"","candidate_id")),c("a","b","c","d"))
  a_rows <- Filter(function(x) x$candidate_id=="a",result$evidence)
  expect_setequal(vapply(a_rows,`[[`,"","metric"),
    c("EL_score","EL_rank","BA_log_score","BA_rank","BA_nM"))
  scores <- stats::setNames(vapply(a_rows,`[[`,0.0,"raw_score"),
    vapply(a_rows,`[[`,"","metric"))
  expect_equal(scores[["EL_score"]],0.901234)
  expect_equal(scores[["EL_rank"]],0.04)
  expect_equal(scores[["BA_log_score"]],0.6123)
  expect_equal(scores[["BA_rank"]],0.9)
  expect_equal(scores[["BA_nM"]],79)
  expect_true(all(vapply(result$evidence,function(x) x$predictor_version=="4.2d",TRUE)))
  expect_equal(Filter(function(x) x$candidate_id=="c" && x$metric=="BA_nM",
    result$evidence)[[1L]]$raw_score,410)
})

test_that("NetMHCpan golden parser retains BA, rank, EL and fanout", {
  cap <- predictor_test_capability("netmhcpan")
  a <- predictor_test_candidate("a"); b <- predictor_test_candidate("b")
  d <- predictor_test_candidate("d","LMNPQRSTV")
  c <- predictor_test_candidate("c","LMNPQRSTV","HLA-B*07:02")
  p <- Neoantimon2:::predictor_prepare(list(c,d,b,a),"netmhcpan","mhc1_binding",cap)
  lines <- readLines(testthat::test_path("..","fixtures","predictor","netmhcpan-4.1.txt"))
  header <- lines[grepl("^Pos ",lines)]
  rows <- lines[grepl("^[0-9]",lines)]
  raw <- list("HLA-A02:01"=paste(c(header,rev(rows[grepl("HLA-A02:01",rows)])),collapse="\n"),
    "HLA-B07:02"=paste(c(header,rows[grepl("HLA-B07:02",rows)]),collapse="\n"))
  result <- Neoantimon2:::predictor_normalize(p,raw,"run-1",list("HLA-A02:01"="a.txt",
    "HLA-B07:02"="b.txt"))
  expect_length(result$evidence,20L) # 5 metrics for 4 biological candidates
  expect_length(result$statuses,0L)
  expect_setequal(unique(vapply(result$evidence,`[[`,"","candidate_id")),c("a","b","c","d"))
  expect_equal(Filter(function(x) x$candidate_id=="d" && x$metric=="BA_nM",result$evidence)[[1L]]$raw_score,410)
  a_rows <- Filter(function(x) x$candidate_id=="a",result$evidence)
  expect_setequal(vapply(a_rows,`[[`,"","metric"),c("EL_score","EL_rank","BA_log_score","BA_rank","BA_nM"))
  expect_identical(Filter(function(x) x$metric=="BA_nM",a_rows)[[1L]]$score_direction,"lower_is_better")
  expect_identical(Filter(function(x) x$metric=="EL_score",a_rows)[[1L]]$task,"mhc1_presentation")
  expect_true(all(vapply(a_rows,function(x) !is.null(x$raw_output_ref),TRUE)))
})

test_that("NetMHCIIpan cores remain predictor-specific and share epitope keys", {
  cap <- predictor_test_capability("netmhciipan")
  a <- predictor_test_candidate("a","ACDEFGHIKLMNPQR","HLA-DQA1*05:01~HLA-DQB1*02:01")
  b <- predictor_test_candidate("b","CDEFGHIKLMNPQRS",a$mhc_key)
  dr <- predictor_test_candidate("dr","ACDEFGHIKLMNPQR","HLA-DRA*01:01~HLA-DRB1*04:01")
  p <- Neoantimon2:::predictor_prepare(list(dr,b,a),"netmhciipan","mhc2_binding",cap)
  lines <- readLines(testthat::test_path("..","fixtures","predictor","netmhciipan-4.3.txt"))
  header <- lines[grepl("^Pos ",lines)]; rows <- lines[grepl("^[0-9]",lines)]
  raw <- list("HLA-DQA10501-DQB10201"=paste(c(header,rows[grepl("HLA-DQA10501",rows)]),collapse="\n"),
    "DRB1_0401"=paste(c(header,rows[grepl("DRB1_0401",rows)]),collapse="\n"))
  result <- Neoantimon2:::predictor_normalize(p,raw,"run-II",list("HLA-DQA10501-DQB10201"="dq.txt",
    "DRB1_0401"="dr.txt"))
  expect_length(result$core_evidence,3L)
  cores <- Filter(function(x) x$candidate_id %in% c("a","b"),result$core_evidence)
  expect_identical(cores[[1L]]$epitope_key,cores[[2L]]$epitope_key)
  expect_identical(cores[[1L]]$binding_core_start,3L)
  expect_true(all(vapply(cores,function(x) x$core_mapping_status=="mapped" &&
    x$epitope_key_status=="derived" && x$binding_orientation=="forward",TRUE)))
  expect_identical(substring(a$ligand_sequence,cores[[1L]]$binding_core_start,
    cores[[1L]]$binding_core_start+nchar(cores[[1L]]$binding_core)-1L),cores[[1L]]$binding_core)
  expect_true(all(vapply(result$core_evidence,function(x) !is.null(x$run_id),TRUE)))
  expect_true(all(vapply(result$evidence,function(x) x$task %in% c("mhc2_binding","mhc2_presentation"),TRUE)))
  changed <- result$core_evidence[[1L]]
  changed$model_version <- "other-model"
  expect_identical(changed$epitope_key,result$core_evidence[[1L]]$epitope_key)
  expect_false(identical(changed$model_version,result$core_evidence[[1L]]$model_version))
})

test_that("NetMHCIIpan inverted core maps to original ligand and defers identity", {
  cap <- predictor_test_capability("netmhciipan")
  ligand <- "PIWNPISEFVKWYKS"
  candidate <- predictor_test_candidate("inverted",ligand,
    "HLA-DPA1*02:02~HLA-DPB1*19:01")
  p <- Neoantimon2:::predictor_prepare(list(candidate),"netmhciipan","mhc2_binding",cap)
  allele <- Neoantimon2:::predictor_allele(candidate,"netmhciipan")
  raw <- predictor_test_fixture("netmhciipan-4.3-inverted.txt")
  result <- Neoantimon2:::predictor_run(p,tempfile("inverted-run-"),
    fixture_raw=stats::setNames(list(raw),allele))
  expect_length(result$statuses,0L)
  expect_length(result$core_evidence,1L)
  core <- result$core_evidence[[1L]]
  expect_identical(core$binding_orientation,"reverse")
  expect_identical(core$binding_core,"KYWKVFESI")
  expect_identical(core$binding_core_start,6L)
  expect_identical(core$core_mapping_status,"mapped")
  expect_identical(core$epitope_key_status,"deferred_orientation")
  expect_null(core$epitope_key)
  expect_identical(paste(rev(strsplit(core$binding_core,"",fixed=TRUE)[[1L]]),collapse=""),
    substring(ligand,core$binding_core_start,
      core$binding_core_start+nchar(core$binding_core)-1L))
  row <- Neoantimon2:::predictor_parse_text(raw,"netmhciipan")$rows[[1L]]
  row$Core <- "KYWKVFESA"
  bad <- Neoantimon2:::predictor_core(row,candidate,"netmhciipan","4.3k",NULL,"run-bad")
  expect_identical(bad$binding_orientation,"reverse")
  expect_identical(bad$core_mapping_status,"unresolved")
  expect_identical(bad$epitope_key_status,"unavailable")
  expect_null(bad$binding_core_start)
  expect_null(bad$epitope_key)
})

test_that("affinity_only is restricted to MHCflurry binding", {
  request <- function(provider,task,options) list(list(provider=provider,task=task,options=options))
  expect_error(Neoantimon2:::predictor_policy(request("netmhcpan","mhc1_binding",
    list(affinity_only=TRUE))),"only for MHCflurry")
  expect_error(Neoantimon2:::predictor_policy(request("netmhcpan","mhc1_binding",
    list(affinity_only=FALSE))),"only for MHCflurry")
  expect_error(Neoantimon2:::predictor_policy(request("mhcflurry","mhc1_presentation",
    list(affinity_only=TRUE))),"binding task")
  expect_silent(Neoantimon2:::predictor_policy(request("mhcflurry","mhc1_binding",
    list(affinity_only=TRUE))))
  p <- Neoantimon2:::predictor_prepare(list(predictor_test_candidate("a")),"mhcflurry",
    "mhc1_presentation",predictor_test_capability("mhcflurry"))
  expect_error(Neoantimon2:::predictor_run(p,tempfile("bad-options-"),
    options=list(affinity_only=TRUE)),"binding task")
  expect_error(Neoantimon2:::predictor_policy(request("mhcflurry","mhc1_binding",
    list(unrecognized=TRUE))),"predictor options")
})

test_that("MHCflurry fixture separates affinity and presentation", {
  cap <- predictor_test_capability("mhcflurry")
  a <- predictor_test_candidate("a")
  p <- Neoantimon2:::predictor_prepare(list(a),"mhcflurry","mhc1_binding",cap)
  raw <- predictor_test_fixture("mhcflurry-2.csv")
  lines <- strsplit(raw,"\n",fixed=TRUE)[[1L]]
  result <- Neoantimon2:::predictor_normalize(p,list("HLA-A*02:01"=paste(c(lines[[1L]],
    lines[grepl("ACDEFGHIK",lines)]),collapse="\n")),"run-F",list("HLA-A*02:01"="f.csv"))
  expect_setequal(vapply(result$evidence,`[[`,"","metric"),
    c("BA_nM","BA_rank","presentation_score","presentation_percentile"))
  expect_setequal(unique(vapply(result$evidence,`[[`,"","task")),
    c("mhc1_binding","mhc1_presentation"))
  expect_length(result$statuses,0L)
  pan_cap <- predictor_test_capability("netmhcpan")
  pan <- Neoantimon2:::predictor_prepare(list(a),"netmhcpan","mhc1_binding",pan_cap)
  pan_lines <- readLines(testthat::test_path("..","fixtures","predictor","netmhcpan-4.1.txt"))
  pan_raw <- paste(c(pan_lines[grepl("^Pos ",pan_lines)],
    pan_lines[grepl("HLA-A02:01 ACDEFGHIK",pan_lines,fixed=TRUE)]),collapse="\n")
  pan_result <- Neoantimon2:::predictor_normalize(pan,list("HLA-A02:01"=pan_raw),
    "run-pan",list("HLA-A02:01"="pan.txt"))
  combined <- c(result$evidence,pan_result$evidence)
  expect_setequal(unique(vapply(combined,`[[`,"","predictor")),c("mhcflurry","netmhcpan"))
  expect_true(all(vapply(combined,function(x) x$candidate_id=="a",TRUE)))
  expect_equal(length(Filter(function(x) x$metric=="BA_nM",combined)),2L)
  expect_false(identical(Filter(function(x) x$metric=="BA_nM",result$evidence)[[1L]]$raw_score,
    Filter(function(x) x$metric=="BA_nM",pan_result$evidence)[[1L]]$raw_score))
})

test_that("malformed and missing rows are explicit without fabricated scores", {
  cap <- predictor_test_capability("netmhcpan")
  a <- predictor_test_candidate("a")
  p <- Neoantimon2:::predictor_prepare(list(a),"netmhcpan","mhc1_binding",cap)
  bad <- Neoantimon2:::predictor_normalize(p,list("HLA-A02:01"="Pos HLA Peptide Score_EL\n0 BAD"),
    "run-bad",list("HLA-A02:01"="raw.txt"))
  expect_length(bad$evidence,0L)
  expect_true("malformed_row" %in% vapply(bad$statuses,`[[`,"","reason"))
  expect_true("missing_output_row" %in% vapply(bad$statuses,`[[`,"","reason"))
  none <- Neoantimon2:::predictor_normalize(p,list(),"run-none",list())
  expect_length(none$evidence,0L)
  expect_identical(none$statuses[[1L]]$reason,"execution_failure")
  numeric_bad <- Neoantimon2:::predictor_normalize(p,list("HLA-A02:01"=
    "Pos HLA Peptide Score_EL %Rank_EL\n0 HLA-A02:01 ACDEFGHIK bad NA"),
    "run-numeric-bad",list("HLA-A02:01"="raw.txt"))
  expect_length(numeric_bad$evidence,0L)
  expect_true("malformed_score" %in% vapply(numeric_bad$statuses,`[[`,"","reason"))
  expect_true("metric_unavailable" %in% vapply(numeric_bad$statuses,`[[`,"","reason"))
})

test_that("fixture adapter preserves raw artifacts, U and predictor-returned core", {
  cap <- predictor_test_capability("fixture")
  seq <- "ACDEFGHUKLMNPQR"
  a <- predictor_test_candidate("sec-I",substring(seq,1L,9L))
  b <- predictor_test_candidate("sec-II",seq,"HLA-DQA1*05:01~HLA-DQB1*02:01")
  prepared <- Neoantimon2:::predictor_prepare(list(b,a),"fixture","mhc2_binding",cap)
  # Requested task is Class II, so the Class I candidate remains an explicit
  # unsupported-task status; its literal-U biological input is untouched.
  expect_length(prepared$queries,1L)
  expect_identical(prepared$statuses[[1L]]$reason,"unsupported_task")
  allele <- b$mhc_key
  raw <- paste("Pos MHC Peptide Of Core Score_EL %Rank_EL Score_BA Affinity(nM) %Rank_BA",
    paste("1",allele,seq,"2","DEFGHUKLM","0.5","1.2","0.4","123","2.0"),sep="\n")
  root <- tempfile("predictor-artifacts-")
  result <- Neoantimon2:::predictor_run(prepared,root,fixture_raw=stats::setNames(list(raw),allele))
  expect_length(result$raw_artifacts,1L)
  expect_true(file.exists(file.path(root,result$run_id,result$raw_artifacts[[1L]]$file)))
  input <- file.path(root,result$run_id,result$commands[[allele]]$input_file)
  expect_true(file.exists(input))
  expect_identical(digest::digest(file=input,algo="sha256",serialize=FALSE),
    result$commands[[allele]]$input_sha256)
  expect_equal(length(result$evidence),5L)
  expect_identical(result$core_evidence[[1L]]$binding_core,"DEFGHUKLM")
  expect_match(result$core_evidence[[1L]]$epitope_key,"^epitope_key:")
  expect_identical(result$core_evidence[[1L]]$binding_core_start,3L)
  expect_identical(Neoantimon2:::predictor_run(prepared,root,fixture_raw=stats::setNames(list(raw),allele))$run_id,
    result$run_id)
  writeLines("tampered",input)
  expect_error(Neoantimon2:::predictor_run(prepared,root,fixture_raw=stats::setNames(list(raw),allele)),
    "input cache integrity")
  expect_null(b$epitope_key)
})

test_that("same ligand may have parallel core evidence from different runs", {
  cap <- predictor_test_capability("netmhciipan")
  b <- predictor_test_candidate("b","ACDEFGHIKLMNPQR","HLA-DQA1*05:01~HLA-DQB1*02:01")
  p <- Neoantimon2:::predictor_prepare(list(b),"netmhciipan","mhc2_binding",cap)
  allele <- Neoantimon2:::predictor_allele(b,"netmhciipan")
  header <- "Pos MHC Peptide Of Core Score_EL %Rank_EL"
  raw1 <- paste(header,paste("1",allele,b$ligand_sequence,"1","CDEFGHIKL","0.2","4.2"),sep="\n")
  raw2 <- paste(header,paste("1",allele,b$ligand_sequence,"2","DEFGHIKLM","0.3","3.1"),sep="\n")
  x <- Neoantimon2:::predictor_normalize(p,stats::setNames(list(raw1),allele),"run-1",
    stats::setNames(list("one.txt"),allele))
  y <- Neoantimon2:::predictor_normalize(p,stats::setNames(list(raw2),allele),"run-2",
    stats::setNames(list("two.txt"),allele))
  expect_false(identical(x$core_evidence[[1L]]$epitope_key,y$core_evidence[[1L]]$epitope_key))
  expect_identical(x$core_evidence[[1L]]$candidate_id,y$core_evidence[[1L]]$candidate_id)
  expect_null(b$epitope_key)
})

test_that("postrun prediction sidecars preserve completed biological runs", {
  project <- peptide_test_run()
  before <- Neoantimon2:::read_manifest(project$path)
  request <- list(list(provider="netmhcpan",task="mhc1_binding",
    executable=file.path(tempdir(),"absent-netmhcpan"),options=list()))
  output <- Neoantimon2:::prediction_postrun(project,request)
  expect_length(output,1L)
  run <- output[[1L]]
  expect_length(run$evidence,0L)
  expect_true(length(run$statuses)>0L)
  expect_true(all(vapply(run$statuses,function(x) x$status %in% c("unsupported","unavailable"),TRUE)))
  expect_true(file.exists(Neoantimon2:::predictor_result_path(project,run$run_id)))
  expect_identical(Neoantimon2:::read_manifest(project$path),before)
  expect_identical(Neoantimon2:::prediction_postrun(project,request)[[1L]]$run_id,run$run_id)
})

test_that("configured prediction stage checkpoints an unavailable provider", {
  p <- predictor_test_project()
  p <- neo_run(p)
  snapshot <- Neoantimon2:::read_runtime_snapshot(p)
  expect_true(length(snapshot$neoantigen_results)>0L)
  expect_true(all(vapply(snapshot$neoantigen_results,function(x) x$prediction_status=="NOT_RUN",TRUE)))
  manifest <- Neoantimon2:::read_manifest(p$path)
  cp <- Filter(function(x) startsWith(x$step_id,"prediction_"),manifest$checkpoints)
  expect_length(cp,1L)
  expect_identical(cp[[1L]]$status,"COMPLETED")
  Neoantimon2:::verify_checkpoint_outputs(cp[[1L]],p$path)
  artifact <- readRDS(file.path(p$path,cp[[1L]]$output_records[[1L]]$path))
  expect_length(artifact$evidence,0L)
  expect_true(any(vapply(artifact$statuses,function(x) x$reason=="unsupported_length",TRUE)))
  expect_identical(neo_run(p)$run_id,p$run_id)
})

test_that("execution failure remains a status and never becomes a score", {
  cap <- predictor_test_capability("netmhcpan")
  cap$executable_path <- Sys.which("false")
  p <- Neoantimon2:::predictor_prepare(list(predictor_test_candidate("failed")),
    "netmhcpan","mhc1_binding",cap)
  run <- Neoantimon2:::predictor_run(p,tempfile("failed-predictor-"))
  expect_length(run$evidence,0L)
  expect_length(run$raw_artifacts,1L)
  expect_identical(run$raw_artifacts[[1L]]$kind,"execution_log")
  expect_identical(run$statuses[[1L]]$reason,"execution_failure")
  expect_identical(run$statuses[[1L]]$candidate_id,"failed")
})
