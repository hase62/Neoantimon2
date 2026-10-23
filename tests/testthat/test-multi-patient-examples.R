# Execute the user-facing script through the public API, using installed tiny
# references. It also runs under R CMD check, where docs/ is excluded.
cohort_example <- new.env(parent = globalenv())
invisible(capture.output(sys.source(system.file("examples", "multi-sample.R",
  package = "Neoantimon2", mustWork = TRUE), envir = cohort_example)))

expect_cohort_associations <- function(p, samples, hla) {
  events <- neo_events(p)
  origins <- neo_origins(p)
  candidates <- neo_results(p)
  for (rows in list(events, origins, candidates)) {
    expect_setequal(unique(rows$sample_id), samples$sample_id)
    expect_identical(rows$patient_id, samples$patient_id[match(rows$sample_id, samples$sample_id)])
    expect_identical(rows$analysis_context, paste0("somatic:", rows$patient_id))
  }
  expect_equal(anyDuplicated(events$event_id), 0L)
  expect_equal(anyDuplicated(origins$origin_id), 0L)
  expect_equal(anyDuplicated(candidates$candidate_id), 0L)
  # Identical variants have shared biological keys, distinct occurrence IDs.
  expect_length(unique(events$event_key), 1L)
  for (s in samples$sample_id) {
    patient <- samples$patient_id[match(s, samples$sample_id)]
    ev <- neo_events(p, sample_id = s)
    ori <- neo_origins(p, sample_id = s)
    res <- neo_results(p, sample_id = s)
    expect_equal(nrow(ev), 1L)
    expect_equal(nrow(ori), 7L)
    expect_equal(nrow(res), 4L * length(hla[[patient]]))
    expect_setequal(res$mhc_key, hla[[patient]])
    for (rows in list(ev, ori, res)) {
      expect_identical(unique(rows$sample_id), s)
      expect_identical(unique(rows$patient_id), patient)
    }
    sources <- ev$event_details_json[[1L]]$source_records
    expect_true(length(sources) > 0L)
    for (source in sources) {
      expect_identical(source$sample_id, s)
      expect_identical(source$input_id, paste0(s, "-somatic"))
    }
    for (i in seq_len(nrow(ori))) {
      ids <- unlist(ori$source_event_ids[[i]])
      expect_true(length(ids) > 0L)
      expect_true(all(ids %in% ev$event_id))
    }
    for (i in seq_len(nrow(res))) {
      ids <- unlist(res$origin_ids[[i]])
      expect_equal(length(ids), res$n_origins[[i]])
      expect_true(all(ids %in% ori$origin_id))
      linked <- ori[match(ids, ori$origin_id), , drop = FALSE]
      expect_true(all(linked$sample_peptide_id == res$sample_peptide_id[[i]]))
      expect_true(all(linked$ligand_sequence == res$ligand_sequence[[i]]))
      expect_true(all(linked$presenter_subject_id == patient))
    }
  }
  expect_identical(unique(candidates$prediction_status), "NOT_RUN")
  expect_identical(unique(candidates$grading_status), "NOT_GRADED")
  expect_identical(unique(candidates$final_grade), "NOT_GRADED")
  # This accessor is not a generic inventory of the project's samples.
  expect_equal(nrow(neo_multisample(p)), 0L)
  reopened <- neo_project(path = p$path)
  expect_identical(neo_events(reopened), events)
  expect_identical(neo_origins(reopened), origins)
  expect_identical(neo_results(reopened), candidates)
}

for (scenario in seq_len(4L)) {
  test_that(paste("documented cohort scenario", scenario, "preserves associations"), {
    expect_cohort_associations(cohort_example[[paste0("p", scenario)]],
      cohort_example[[paste0("samples_", scenario)]],
      cohort_example[[paste0("hla_", scenario)]])
  })
}

test_that("shared patient and cross-patient candidates are never deduplicated by key", {
  p <- cohort_example$p4
  a <- neo_results(p, sample_id = "S1")
  b <- neo_results(p, sample_id = "S2")
  c <- neo_results(p, sample_id = "S3")
  shared_c <- c[c$mhc_key == "HLA-A*02:01", , drop = FALSE]
  expect_setequal(a$candidate_key, b$candidate_key)
  expect_setequal(a$candidate_key, shared_c$candidate_key)
  expect_length(intersect(a$candidate_id, b$candidate_id), 0L)
  expect_length(intersect(a$candidate_id, c$candidate_id), 0L)
  expect_false("HLA-B*07:02" %in% c(a$mhc_key, b$mhc_key))
  # Ordering declarations cannot change project identity or HLA assignment.
  reordered <- neo_project(samples = cohort_example$samples_4[3:1, ],
    hla = rev(cohort_example$hla_4), inputs = rev(cohort_example$inputs_4),
    config = cohort_example$cfg)
  expect_identical(reordered$analysis_id, p$analysis_id)
  expect_identical(reordered$hla[sort(names(reordered$hla))], p$hla[sort(names(p$hla))])
})

test_that("asymmetric sample variants agree with separate runs without cross-sample reconstruction", {
  other <- file.path(cohort_example$input_root, "S2-asymmetric.vcf")
  writeLines(sub("chr1\t4\t.\tG\tA", "chr1\t5\t.\tA\tC",
    readLines(file.path(cohort_example$input_root, "S2.vcf")), fixed = TRUE), other)
  inputs <- cohort_example$inputs_2
  inputs[[2L]]$path <- other
  joint <- neo_run(neo_project(samples = cohort_example$samples_2,
    hla = cohort_example$hla_2, inputs = inputs, config = cohort_example$cfg))
  for (s in c("S1", "S2")) {
    single <- neo_run(neo_project(
      samples = data.frame(sample_id = s, patient_id = "P1"),
      hla = cohort_example$hla_2,
      inputs = Filter(function(x) x$sample_id == s, inputs), config = cohort_example$cfg))
    for (spec in list(
      list(fun = neo_events, fields = c("event_id", "event_key", "sample_id", "patient_id", "pos", "ref", "alt")),
      list(fun = neo_origins, fields = c("origin_id", "sample_id", "patient_id", "peptide_sequence", "source_event_ids")),
      list(fun = neo_results, fields = c("candidate_id", "candidate_key", "sample_id", "patient_id", "mhc_key", "origin_ids")))) {
      a <- spec$fun(joint, sample_id = s)[, spec$fields, drop = FALSE]
      b <- spec$fun(single, sample_id = s)[, spec$fields, drop = FALSE]
      rownames(a) <- rownames(b) <- NULL
      expect_true(nrow(a) > 0L)
      expect_identical(a, b)
    }
  }
})

test_that("missing patient HLA never borrows another patient's presenters", {
  p <- neo_run(neo_project(samples = cohort_example$samples_3,
    hla = list(P1 = "HLA-A*02:01"), inputs = cohort_example$inputs_3,
    config = cohort_example$cfg))
  expect_equal(nrow(neo_events(p, sample_id = "S2")), 1L)
  expect_equal(nrow(neo_origins(p, sample_id = "S2")), 7L)
  expect_equal(nrow(neo_results(p, sample_id = "S2")), 0L)
  expect_identical(unique(neo_results(p)$patient_id), "P1")
})

test_that("invalid project declarations fail under the existing validation contract", {
  make <- function(samples = cohort_example$samples_2,
                   hla = cohort_example$hla_2, inputs = cohort_example$inputs_2) {
    neo_project(samples = samples, hla = hla, inputs = inputs, config = cohort_example$cfg)
  }
  duplicate <- cohort_example$samples_2
  duplicate$sample_id[[2L]] <- "S1"
  expect_error(make(samples = duplicate), "Duplicate sample_id")
  duplicate$patient_id[[2L]] <- "P2"
  expect_error(make(samples = duplicate), "Duplicate sample_id")
  expect_error(make(samples = data.frame(sample_id = "S1")), "patient_id must be a nonempty string")
  for (bad in list("", NA_character_, 1L)) {
    expect_error(make(samples = data.frame(sample_id = "S1", patient_id = bad)),
      "patient_id must be a nonempty string")
  }
  unknown <- cohort_example$inputs_2
  unknown[[1L]]$sample_id <- "unknown"
  expect_error(make(inputs = unknown), "Input sample binding missing")
  duplicate <- cohort_example$inputs_2
  duplicate[[2L]]$input_id <- duplicate[[1L]]$input_id
  expect_error(make(inputs = duplicate), "Duplicate input_id")
  expect_error(make(hla = list(P2 = "HLA-B*07:02")), "Unknown key in hla subject")
  expect_error(make(hla = list(S1 = "HLA-A*02:01")), "Unknown key in hla subject")
  expect_error(make(hla = list(P1 = NA_character_)), "declared allele strings")
  missing <- cohort_example$inputs_2
  missing[[1L]]$path <- "not-created.vcf"
  expect_error(make(inputs = missing), "Required input file missing")
  extra <- cohort_example$samples_2
  extra$patient <- "P1"
  expect_error(make(samples = extra), "Unknown key in sample")
  wrong_context <- cohort_example$p3$analysis_context
  wrong_context[[1L]]$sample_bindings <- list("S2")
  expect_error(neo_project(samples = cohort_example$samples_3,
    analysis_context = wrong_context, config = cohort_example$cfg),
    "Invalid context sample identity binding")
})

test_that("VCF header mismatches and multiple sample columns fail at execution", {
  inputs <- cohort_example$inputs_1
  inputs[[1L]]$path <- "S2.vcf"
  p <- neo_project(samples = cohort_example$samples_1, inputs = inputs,
    config = cohort_example$cfg)
  expect_error(neo_run(p), "VCF sample does not match authoritative input sample binding")
  path <- file.path(cohort_example$input_root, "two-columns.vcf")
  lines <- readLines(file.path(cohort_example$input_root, "S1.vcf"))
  lines[startsWith(lines, "#CHROM")] <- paste0(lines[startsWith(lines, "#CHROM")], "\tS2")
  lines[startsWith(lines, "chr1")] <- paste0(lines[startsWith(lines, "chr1")], "\t0/1")
  writeLines(lines, path)
  inputs[[1L]]$path <- path
  p <- neo_project(samples = cohort_example$samples_2, inputs = inputs,
    config = cohort_example$cfg)
  expect_error(neo_run(p), "Multi-sample VCF requires a single explicit sample input")
})

test_that("ordinary timepoint labels do not activate longitudinal semantics", {
  samples <- cohort_example$samples_2
  samples$timepoint <- c("baseline", "followup")
  p <- neo_run(neo_project(samples = samples, hla = cohort_example$hla_2,
    inputs = cohort_example$inputs_2, config = cohort_example$cfg))
  expect_equal(nrow(neo_multisample(p, patient_id = "P1")), 0L)
  expect_equal(nrow(neo_results(p)), 8L)
  expect_null(Neoantimon2:::longitudinal_settings(p))
  expect_false(dir.exists(file.path(p$path, "results", "longitudinal")))
  expect_setequal(neo_results(p)$candidate_id, neo_results(cohort_example$p2)$candidate_id)
})
