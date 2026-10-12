test_that("package loads side by side and schema is the complete field authority", {
  expect_true(isNamespaceLoaded("Neoantimon2"))
  # The exact two-export assertion was milestone-specific; T004 integrates
  # additional approved APIs while the schema/ID contracts stay frozen.
  expect_length(setdiff(getNamespaceExports("Neoantimon2"), approved_public_api), 0L)
  expect_setequal(names(contract$tables), c("events","peptide_origins","neoantigen_results","multisample_summary","run_manifest"))
  for (def in c(contract$tables,contract$records,contract$identities)) {
    expect_false(anyDuplicated(vapply(def$fields,`[[`,"","name")) > 0)
    for (f in def$fields) {
      expect_true(all(c("name","type","required","nullable","unit","definition") %in% names(f)))
      expect_true(nzchar(f$definition))
    }
  }
  expect_true(length(Neoantimon2:::schema_field_markdown(contract)) > 100)
})

test_that("independent synthetic fixtures validate against frozen human truth", {
  truth <- fixture("expected")
  expect_identical(truth$fixture_kind,"synthetic")
  for (n in names(truth$cases)) {
    x <- fixture(n); ex <- truth$cases[[n]]
    expect_true(valid(x), info=n)
    expect_length(x$events,ex$events)
    expect_length(x$peptide_origins,ex$origins)
    expect_length(x$neoantigen_results,ex$candidates)
    expect_identical(x$run_manifest$fixture_kind,"synthetic")
  }
})

test_that("empty tables, required presence, nullable values and unknown fields", {
  x <- fixture("simple_somatic")
  for (t in c("events","peptide_origins","neoantigen_results")) x[[t]] <- list()
  expect_true(valid(x))
  x <- fixture("simple_somatic"); x$events[[1]]$tumor_vaf <- NULL
  expect_error(valid(x),"required field tumor_vaf")
  x <- fixture("simple_somatic"); x$events[[1]]["tumor_vaf"] <- list(NULL)
  expect_true(valid(x))
  x$events[[1]]$event_id <- NA_character_
  expect_error(valid(x),"null is not allowed")
  x <- fixture("simple_somatic"); x$events[[1]]$unexpected <- 1
  expect_error(valid(x),"unknown fields")
  expect_error(Neoantimon2:::validate_table(data.frame(),"events",contract),"Missing required")
  x <- fixture("simple_somatic"); x$events[[2]] <- x$events[[1]]
  expect_error(valid(x),"duplicate primary key")
})

test_that("enum, units and finite biological ranges are enforced", {
  for (n in c("tumor_vaf","normal_vaf","rna_vaf","purity","ccf","ccf_lower","ccf_upper")) {
    for (value in c(-0.01,1.01,Inf,NaN)) {
      x <- fixture("simple_somatic"); x$events[[1]][[n]] <- value
      expect_error(valid(x),info=paste(n,value))
    }
  }
  x <- fixture("simple_somatic"); x$events[[1]]$tumor_count_unit <- "reads"
  expect_error(valid(x),"enum")
  x <- fixture("simple_somatic"); x$events[[1]]$tumor_alt_count <- 0.5
  expect_error(valid(x),"integer")
  x$events[[1]]$tumor_alt_count <- -1L
  expect_error(valid(x),"minimum")
  for (rank in c(0,2,100)) {
    x <- fixture("simple_somatic"); x$neoantigen_results[[1]]$prediction_status <- "SUCCESS"
    x$neoantigen_results[[1]]$mhc_rank_pct <- rank
    expect_true(valid(x))
  }
  for (rank in c(-1,100.1)) {
    x <- fixture("simple_somatic"); x$neoantigen_results[[1]]$mhc_rank_pct <- rank
    expect_error(valid(x))
  }
  x <- fixture("simple_somatic"); x$neoantigen_results[[1]]$mhc_affinity_nm <- -1
  expect_error(valid(x),"minimum")
  expect_equal(Neoantimon2:::rank_percent(.02,"fraction_0_1"),2)
  expect_equal(Neoantimon2:::rank_percent(2,"percent_0_100"),2)
  expect_error(Neoantimon2:::rank_percent(.02),"Declare")
  expect_error(Neoantimon2:::rank_percent(2,"fraction_0_1"),"range")
})

test_that("zero, missing, unsupported and low confidence have separate meanings", {
  x <- fixture("zero_vs_missing")
  expect_identical(x$events[[1]]$rna_alt_count,0L)
  expect_null(x$events[[1]]$normal_alt_count)
  expect_true(valid(x))
  x$events[[1]]$rna_assessment_status <- "NOT_RUN"
  expect_error(valid(x),"EVALUATED")
  x <- fixture("zero_vs_missing"); x$neoantigen_results[[1]]$evidence_records_json[[1]]$value <- FALSE
  expect_error(valid(x),"cannot become negative")
  x <- fixture("zero_vs_missing"); x$neoantigen_results[[1]]$evidence_records_json[[2]]$details <- structure(list(), names=character())
  expect_error(valid(x),"scope")
  x <- fixture("simple_somatic"); x$neoantigen_results[[1]]$prediction_status <- "UNSUPPORTED"
  x$neoantigen_results[[1]]$not_evaluated_codes <- list("UNSUPPORTED_INPUT")
  x$neoantigen_results[[1]]$warning_codes <- list("LOW_CONFIDENCE")
  expect_true(valid(x))
  x$neoantigen_results[[1]]$mhc_rank_pct <- 0
  expect_error(valid(x),"Unperformed prediction")
  x <- fixture("simple_somatic"); x$neoantigen_results[[1]]$final_grade <- "INVALID"
  expect_error(valid(x),"Valid candidate")
  for (grade in c("A","B","C","D","E")) {
    x <- fixture("simple_somatic"); x$neoantigen_results[[1]]$final_grade <- grade
    x$neoantigen_results[[1]]$grading_status <- "GRADED"
    x$neoantigen_results[[1]]$prediction_status <- "SUCCESS"
    expect_true(valid(x)) # Contract validates representation, never computes grades.
  }
  for (n in c("final_grade","grading_status","prediction_status")) {
    x <- fixture("simple_somatic"); x$neoantigen_results[[1]][[n]] <- "BOGUS"
    expect_error(valid(x),"enum")
  }
  x <- fixture("simple_somatic"); x$neoantigen_results[[1]]$warning_codes <- list("NOT_RUN")
  expect_error(valid(x),"enum")
})

test_that("identity oracle, Unicode and canonical typing are stable", {
  for (z in fixture("canonical_vectors")) {
    bytes <- Neoantimon2:::canonical_bytes(z$value)
    expect_identical(bytes,z$canonical)
    expect_identical(digest::digest(charToRaw(bytes),algo="sha256",serialize=FALSE),z$sha256)
  }
  canon <- Neoantimon2:::canonical_bytes
  expect_identical(canon("\u00e9"),canon("e\u0301"))
  expect_identical(canon(NA_character_),canon(NULL))
  expect_false(identical(canon(NULL),canon("")))
  expect_false(identical(canon(1L),canon(1)))
  expect_identical(canon(list(b=1L,a="x")),canon(list(a="x",b=1L)))
  expect_error(canon(Inf),"Non-finite")
  expect_error(canon(NaN),"Non-finite")
  x <- fixture("multiple_origins")
  before <- x$neoantigen_results[[1]]$candidate_id
  x$events <- rev(x$events); x$peptide_origins <- rev(x$peptide_origins)
  x$neoantigen_results[[1]]$origin_ids <- rev(x$neoantigen_results[[1]]$origin_ids)
  expect_true(valid(x)); expect_identical(before,x$neoantigen_results[[1]]$candidate_id)
  o <- x$peptide_origins[[1]]; v <- o[c("sample_peptide_id","origin_identity")]
  v$origin_identity$transcript_version <- "2"
  expect_false(identical(id("origin_id",v),o$origin_id))
  expect_identical(before,x$neoantigen_results[[1]]$candidate_id)
  for (t in c("events","peptide_origins","neoantigen_results")) for (i in seq_along(x[[t]])) x[[t]][[i]]$analysis_id <- "rerun"
  x$run_manifest$analysis_id <- "rerun"
  expect_true(valid(x))
})

test_that("candidate and epitope identities preserve shared cores and registers", {
  x <- fixture("class_ii_shared_core"); r <- x$neoantigen_results
  expect_false(identical(r[[1]]$candidate_id,r[[2]]$candidate_id))
  expect_identical(r[[1]]$epitope_key,r[[2]]$epitope_key)
  v <- r[[1]][c("species","hla_class","binding_core","mhc_key")]
  v$mhc_key <- "HLA-DRA*01:01~HLA-DRB1*04:01"
  expect_false(identical(id("epitope_key",v),r[[1]]$epitope_key))
  x <- fixture("multiple_registers"); r <- x$neoantigen_results[[1]]
  expect_length(r$epitope_assignments_json,2)
  expect_length(unique(vapply(r$epitope_assignments_json,`[[`,"","epitope_key")),2)
  expect_null(r$binding_core);expect_true(valid(x))
  x$neoantigen_results[[1]]$epitope_assignments_json[[2]]$binding_core_start <- 99L
  expect_error(valid(x),"coordinates")
  x <- fixture("multiple_registers"); x$neoantigen_results[[1]]$epitope_assignments_json[[1]]$candidate_id <- "other"
  expect_error(valid(x),"candidate mismatch")
  x <- fixture("class_ii_shared_core"); x$neoantigen_results[[1]]$binding_core <- "CDEFGHIKL"
  expect_error(valid(x),"projection mismatch")
  x <- fixture("class_ii_shared_core"); x$neoantigen_results[[1]]$epitope_assignments_json[[1]]$model_checksum <- "other-model"
  expect_error(valid(x),"provenance")
})

test_that("subject roles, foreign keys and multiple-origin provenance cannot drift", {
  x <- fixture("subject_roles"); expect_true(valid(x))
  expect_false("source_subject_id" %in% names(x$events[[1]]))
  x$events[[1]]$source_subject_role <- "recipient"
  expect_error(valid(x),"Role projection")
  x <- fixture("subject_roles"); x$run_manifest$analysis_contexts[[1]]$presenter_subject_role <- "donor"
  expect_error(valid(x),"direction")
  x <- fixture("simple_somatic"); x$peptide_origins[[1]]$source_event_ids <- list("nonexistent")
  expect_error(valid(x),"foreign key")
  x <- fixture("simple_somatic"); x$events[[1]]$analysis_context <- "other"
  expect_error(valid(x),"foreign key")
  x <- fixture("multiple_origins"); x$neoantigen_results[[1]]$origin_ids <- x$neoantigen_results[[1]]$origin_ids[1]
  expect_error(valid(x),"origins must be preserved")
  x <- fixture("multiple_origins"); x$neoantigen_results[[1]]$rna_alt_count <- 9L
  expect_error(valid(x),"Cross-origin")
  x <- fixture("multiple_origins"); x$neoantigen_results[[1]]$gene_tpm <- 20.0
  expect_error(valid(x),"Cross-origin")
})

test_that("schema version tracking rejects implicit migrations", {
  for (t in c("events","peptide_origins","neoantigen_results")) {
    x <- fixture("simple_somatic"); x[[t]][[1]]$schema_version <- "0.2"
    expect_error(valid(x),"enum")
  }
  x <- fixture("simple_somatic"); x$run_manifest$id_schema_version <- "2.0.0"
  expect_error(valid(x),"enum")
  x <- fixture("simple_somatic"); x$run_manifest$output_schema_versions$events <- "2.0.0"
  expect_error(valid(x),"version mismatch")
})

test_that("typed evidence units, targets and unresolvable registers are explicit", {
  x <- fixture("zero_vs_missing")
  ev <- x$neoantigen_results[[1]]$evidence_records_json[[2]]
  ev$evidence_type <- "hla_binding_rank"; ev$value <- 2.0; ev$unit <- "percent_0_100"
  x$neoantigen_results[[1]]$evidence_records_json[[2]] <- ev
  expect_true(valid(x))
  x$neoantigen_results[[1]]$evidence_records_json[[2]]$unit <- "fraction_0_1"
  expect_error(valid(x),"unit mismatch")
  x <- fixture("zero_vs_missing"); x$neoantigen_results[[1]]$evidence_records_json[[2]]$target_id <- "other-candidate"
  expect_error(valid(x),"target foreign key")
  x <- fixture("multiple_registers")
  z <- x$neoantigen_results[[1]]$epitope_assignments_json[[1]]
  z$status <- "NOT_EVALUABLE"; z$reason <- "UNRESOLVED_MAPPING"
  z$binding_orientation <- "reverse"; z["epitope_key"] <- list(NULL)
  z$raw_offsets <- list(offset=99L, orientation="inverted")
  x$neoantigen_results[[1]]$epitope_assignments_json[[1]] <- z
  expect_true(valid(x))
  x$neoantigen_results[[1]]$epitope_assignments_json[[1]]$epitope_key <- "fake"
  expect_error(valid(x),"Unresolved assignment")
  x <- fixture("subject_roles"); x$run_manifest$analysis_contexts[[1]]$mode <- "allo_genome"
  x$peptide_origins[[1]]$origin_type <- "allo_genome"
  expect_true(valid(x)) # Storage contract only; no algorithm exists.
})

test_that("counts and conditional CCF preserve uncertainty without clipping", {
  x <- fixture("simple_somatic"); x$events[[1]]$tumor_vaf <- 0
  expect_true(valid(x))
  x$events[[1]]$tumor_vaf <- 1
  expect_true(valid(x))
  x$events[[1]]$tumor_vaf_denominator <- 0L
  expect_error(valid(x),"Zero denominator")
  x <- fixture("simple_somatic"); x$events[[1]]$ccf_raw <- 1.2
  x$events[[1]]$ccf_status <- "OUT_OF_RANGE"
  expect_true(valid(x)); expect_null(x$events[[1]]$ccf)
  x$events[[1]]$ccf_status <- "CONDITIONAL"; x$events[[1]]$ccf <- .8
  x$events[[1]]$ccf_identifiability <- "CONDITIONAL"
  expect_true(valid(x))
  x$events[[1]]$ccf <- 1.2
  expect_error(valid(x),"maximum")
})

test_that("legacy self fields and grades cannot imply evaluation when not run", {
  x <- fixture("simple_somatic"); x$neoantigen_results[[1]]$self_exact_match <- FALSE
  expect_error(valid(x),"evaluated search")
  x$neoantigen_results[[1]]$self_status <- "EVALUATED"
  expect_true(valid(x))
  x <- fixture("simple_somatic"); x$neoantigen_results[[1]]$final_grade <- "E"
  x$neoantigen_results[[1]]$grading_status <- "GRADED"
  expect_error(valid(x),"requires NOT_GRADED")
})

test_that("linked recognition cannot transfer B evidence to a different source antigen", {
  x <- fixture("linked_antigens")
  expect_true(valid(x))
  expect_length(x$neoantigen_results[[1]]$origin_ids,2)
  x$neoantigen_results[[1]]$linked_recognition_json[[1]]$t_cell_origin_ids <- list(x$peptide_origins[[2]]$origin_id)
  expect_error(valid(x),"crosses source antigens")
  x <- fixture("linked_antigens")
  x$neoantigen_results[[1]]$linked_recognition_json[[1]]$b_cell_evidence_ids <- list("unresolved")
  expect_error(valid(x),"foreign key")
})
