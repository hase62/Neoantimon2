test_that("hand-authored reference CDS, phase joins, and all primary isoforms agree", {
  p <- reconstruction_test_project()
  b <- Neoantimon2:::reference_bundle(p)
  r <- Neoantimon2:::event_reference_open(b)
  on.exit(Neoantimon2:::event_reference_close(r))
  models <- Neoantimon2:::reconstruction_models(b, c("TA.1", "TB.1", "TC.1"))
  catalog <- Neoantimon2:::reconstruction_fasta_select(b$resources$protein_fasta$resolved_path,
    c("PA.1", "PB.1", "PC.1"))
  transcript_catalog <- Neoantimon2:::reconstruction_fasta_select(
    b$resources$transcript_fasta$resolved_path, c("TA.1", "TB.1", "TC.1"))
  expected <- list(TA.1 = c("ATGGAACCTAAACCCGGG", "MEPKPG"),
                   TB.1 = c("ATGGATAAACCCGGG", "MDKPG"),
                   TC.1 = c("ATGAAACCTGGGAAACCC", "MKPGKP"))
  for (id in names(expected)) {
    x <- Neoantimon2:::reconstruction_reference(models[[id]], r, catalog, character(),
      transcript_catalog)
    expect_null(x$reason)
    expect_identical(c(x$cds, x$protein), expected[[id]])
    expect_identical(x$transcript, transcript_catalog[[id]])
    expect_identical(x$terminal_codon, "TAA")
  }
  # TA exon 1 contributes CC; exon 2 contributes T to the split CCT codon.
  expect_identical(substr(expected$TA.1[[1L]], 7L, 9L), "CCT")
  expect_identical(Neoantimon2:::reconstruction_translate("ATGGAACCGAAACCCGGG"), "MEPKPG")
})

test_that("plus, minus, MNV, junction, and isoform consequences use sequence", {
  p <- reconstruction_test_project(); b <- Neoantimon2:::reference_bundle(p)
  run <- function(...) Neoantimon2:::reconstruction_contexts(p, b, list(...))
  get <- function(xs, id) Filter(function(x) x$transcript_id == id, xs)[[1L]]
  plus <- get(run(reconstruction_test_row("s1", "chr1", 4L, "G", "A")), "TA.1")
  expect_identical(plus$tumor_cds, "ATGAAACCTAAACCCGGG")
  expect_identical(plus$tumor_protein, "MKPKPG")
  expect_identical(plus$tumor_transcript, "ATGAAACCTAAACCCGGGTAA")
  minus <- get(run(reconstruction_test_row("s2", "chr2", 25L, "T", "C")), "TC.1")
  expect_identical(minus$tumor_cds, "ATGGAACCTGGGAAACCC")
  expect_identical(minus$tumor_protein, "MEPGKP")
  expect_identical(minus$tumor_transcript, "ATGGAACCTGGGAAACCCTAA")
  mnv <- get(run(reconstruction_test_row("m", "chr1", 4L, "GA", "AC")), "TA.1")
  expect_identical(mnv$tumor_cds, "ATGACACCTAAACCCGGG")
  expect_identical(mnv$tumor_protein, "MTPKPG")
  junction <- run(reconstruction_test_row("j", "chr1", 21L, "T", "G"))
  expect_setequal(vapply(junction, `[[`, "", "transcript_id"), c("TA.1", "TB.1"))
  expect_identical(get(junction, "TA.1")$tumor_protein, "MEPKPG")
  expect_identical(get(junction, "TB.1")$tumor_protein, "MEKPG")
  minus_join <- get(run(reconstruction_test_row("mj", "chr2", 13L, "A", "C")), "TC.1")
  expect_identical(minus_join$tumor_cds, "ATGAAACCGGGGAAACCC")
  expect_identical(minus_join$tumor_protein, "MKPGKP")
  expect_identical(get(junction, "TA.1")$mane_status, "MANE Select")
  expect_length(get(junction, "TB.1")$mane_status, 0L)
})

test_that("hand-authored in-frame insertions, deletions, and replacement translate", {
  p <- reconstruction_test_project(); b <- Neoantimon2:::reference_bundle(p)
  run <- function(row, id) Filter(function(x) x$transcript_id == id,
    Neoantimon2:::reconstruction_contexts(p, b, list(row)))[[1L]]
  cases <- list(
    list(reconstruction_test_row("pi", "chr1", 6L, "A", "ACCC"), "TA.1",
         "ATGGAACCCCCTAAACCCGGG", "MEPPKPG"),
    list(reconstruction_test_row("pd", "chr1", 3L, "GGAA", "G"), "TA.1",
         "ATGCCTAAACCCGGG", "MPKPG"),
    list(reconstruction_test_row("pr", "chr1", 4L, "GAA", "GTT"), "TA.1",
         "ATGGTTCCTAAACCCGGG", "MVPKPG"),
    list(reconstruction_test_row("mi", "chr2", 25L, "T", "TCCC"), "TC.1",
         "ATGGGGAAACCTGGGAAACCC", "MGKPGKP"),
    list(reconstruction_test_row("md", "chr2", 22L, "GTTT", "G"), "TC.1",
         "ATGCCTGGGAAACCC", "MPGKP"))
  for (case in cases) {
    x <- run(case[[1L]], case[[2L]])
    expect_identical(x$status, "COMPLETE")
    expect_identical(x$tumor_cds, case[[3L]])
    expect_identical(x$tumor_protein, case[[4L]])
  }
  # GAA -> GTCGAC trims to AA -> TCGAC: a true +3 nt replacement.
  replacement <- reconstruction_test_row("true_replacement", "chr1", 4L,
    "GAA", "GTCGAC")
  edit <- Neoantimon2:::reconstruction_edit_core(replacement)
  expect_identical(replacement$event_type, "INDEL")
  expect_identical(edit$old, "AA")
  expect_identical(edit$new, "TCGAC")
  expect_identical(edit$delta, 3L)
  replaced <- run(replacement, "TA.1")
  expect_identical(replaced$status, "COMPLETE")
  expect_identical(replaced$tumor_cds, "ATGGTCGACCCTAAACCCGGG")
  expect_identical(replaced$tumor_protein, "MVDPKPG")
  inserted <- reconstruction_test_row("ins", "chr1", 6L, "A", "ACCC", side = "LEFT")
  adjacent <- reconstruction_test_row("sub", "chr1", 7L, "C", "T", side = "LEFT")
  combined <- Filter(function(x) x$transcript_id == "TA.1",
    Neoantimon2:::reconstruction_contexts(p, b, list(adjacent, inserted)))[[1L]]
  expect_identical(combined$tumor_cds, "ATGGAACCCTCTAAACCCGGG")
  expect_identical(combined$tumor_protein, "MEPSKPG")
  expect_identical(Filter(function(x) x$transcript_id == "TA.1",
    Neoantimon2:::reconstruction_contexts(p, b, list(inserted, adjacent)))[[1L]], combined)
})

test_that("cis events combine, trans and unknown remain independent", {
  p <- reconstruction_test_project(); b <- Neoantimon2:::reference_bundle(p)
  context <- function(rows) Filter(function(x) x$transcript_id == "TA.1",
    Neoantimon2:::reconstruction_contexts(p, b, rows))
  a <- reconstruction_test_row("a", "chr1", 4L, "G", "A", side = "LEFT")
  c <- reconstruction_test_row("c", "chr1", 5L, "A", "C", side = "LEFT")
  mnv <- reconstruction_test_row("m", "chr1", 4L, "GA", "AC")
  both <- context(list(a, c))
  expect_length(both, 1L)
  expect_identical(both[[1L]]$tumor_cds, "ATGACACCTAAACCCGGG")
  expect_identical(both[[1L]]$tumor_protein, "MTPKPG")
  expect_identical(both[[1L]]$tumor_cds, context(list(mnv))[[1L]]$tumor_cds)
  expect_identical(both[[1L]]$tumor_protein, context(list(mnv))[[1L]]$tumor_protein)
  c$event_details_json$phase_assertions[[1L]]$local_side <- "RIGHT"
  trans <- context(list(a, c))
  expect_length(trans, 2L)
  expect_setequal(vapply(trans, `[[`, "", "tumor_protein"), c("MKPKPG", "MAPKPG"))
  c$event_details_json$phase_assertions <- list()
  unknown <- context(list(a, c))
  expect_length(unknown, 2L)
  expect_setequal(vapply(unknown, `[[`, "", "tumor_protein"), c("MKPKPG", "MAPKPG"))
  expect_identical(context(list(c, a)), unknown)
  a$event_details_json$phase_assertions[[1L]]$local_side <- "RIGHT"
  c$event_details_json$phase_assertions[[1L]] <- list(status = "INTERPRETABLE",
    genotype_state = "LEFT", local_side = "LEFT", block_id = "block:fixture")
  expect_setequal(vapply(context(list(a, c)), `[[`, "", "tumor_cds"),
                  vapply(trans, `[[`, "", "tumor_cds"))
  minus_a <- reconstruction_test_row("ma", "chr2", 24L, "T", "G", side = "LEFT")
  minus_b <- reconstruction_test_row("mb", "chr2", 25L, "T", "C", side = "LEFT")
  minus_mnv <- reconstruction_test_row("mm", "chr2", 24L, "TT", "GC")
  minus_context <- function(rows) Filter(function(x) x$transcript_id == "TC.1",
    Neoantimon2:::reconstruction_contexts(p, b, rows))[[1L]]
  split_minus <- minus_context(list(minus_a, minus_b))
  combined_minus <- minus_context(list(minus_mnv))
  expect_identical(split_minus$tumor_cds, "ATGGCACCTGGGAAACCC")
  expect_identical(split_minus$tumor_protein, "MAPGKP")
  expect_identical(split_minus$tumor_cds, combined_minus$tumor_cds)
  expect_identical(split_minus$tumor_protein, combined_minus$tumor_protein)
})

test_that("patient-normal personalization precedes somatic edits", {
  p <- reconstruction_test_project(); b <- Neoantimon2:::reference_bundle(p)
  context <- function(rows) Filter(function(x) x$transcript_id == "TA.1",
    Neoantimon2:::reconstruction_contexts(p, b, rows))[[1L]]
  germ <- reconstruction_test_row("g", "chr1", 4L, "G", "A", "germline", "LEFT")
  som <- reconstruction_test_row("s", "chr1", 5L, "A", "C", "somatic", "LEFT")
  cis <- context(list(germ, som))
  expect_identical(cis$reference_protein, "MEPKPG")
  expect_identical(cis$patient_normal_cds, "ATGAAACCTAAACCCGGG")
  expect_identical(cis$patient_normal_protein, "MKPKPG")
  expect_identical(cis$tumor_protein, "MTPKPG")
  som$event_details_json$phase_assertions[[1L]]$local_side <- "RIGHT"
  trans <- context(list(germ, som))
  expect_identical(trans$patient_normal_protein, "MEPKPG")
  expect_identical(trans$tumor_protein, "MAPKPG")
  som$event_details_json$phase_assertions <- list()
  unknown <- context(list(germ, som))
  expect_identical(unknown$phase_state, "BACKGROUND_UNKNOWN")
  expect_identical(unknown$tumor_protein, "MAPKPG")
  germ$event_details_json$phase_assertions[[1L]]$genotype_state <- "BOTH"
  germ$event_details_json$phase_assertions[[1L]]$local_side <- NULL
  homo <- context(list(germ, som))
  expect_identical(homo$patient_normal_protein, "MKPKPG")
  expect_identical(homo$tumor_protein, "MTPKPG")
  same <- reconstruction_test_row("s_same", "chr1", 4L, "G", "A", "somatic", "LEFT")
  same_context <- context(list(germ, same))
  expect_identical(same_context$patient_normal_protein, same_context$tumor_protein)
  expect_identical(same_context$protein_change$removed, "")
  germ$event_details_json$phase_assertions[[1L]]$genotype_state <- "LEFT"
  germ$event_details_json$phase_assertions[[1L]]$local_side <- "LEFT"
  third <- reconstruction_test_row("s_third", "chr1", 4L, "G", "C", "somatic", "LEFT")
  layered <- context(list(germ, third))
  expect_identical(layered$patient_normal_protein, "MKPKPG")
  expect_identical(layered$tumor_protein, "MQPKPG")
})

test_that("orientation flips and cross-source anchor loss preserve uncertainty", {
  p <- reconstruction_test_project(); b <- Neoantimon2:::reference_bundle(p)
  contexts <- function(rows) Filter(function(x) x$transcript_id == "TA.1",
    Neoantimon2:::reconstruction_contexts(p, b, rows))
  a <- reconstruction_test_row("a", "chr1", 4L, "G", "A", side = "LEFT")
  c <- reconstruction_test_row("c", "chr1", 5L, "A", "C", side = "LEFT")
  original <- contexts(list(a, c))
  a$event_details_json$phase_assertions[[1L]]$local_side <- "RIGHT"
  c$event_details_json$phase_assertions[[1L]]$local_side <- "RIGHT"
  expect_identical(contexts(list(a, c)), original)
  c$event_details_json$phase_assertions[[1L]]$block_id <- "block:second"
  expect_length(contexts(list(a, c)), 2L)
  a$phase_links_json <- list(list(kind = "block_orientation", source_block = "block:fixture",
    target_block = "block:second", orientation = "SAME", evidence = "SHARED_MARKER"))
  expect_length(contexts(list(a, c)), 1L)
  a$phase_links_json[[1L]]$orientation <- "CONFLICT"
  expect_length(contexts(list(a, c)), 2L)
  local <- reconstruction_test_row("local", "chr1", 7L, "C", "T", side = "RIGHT")
  expect_true(any(vapply(contexts(list(a, local, c)), function(x)
    setequal(x$somatic_event_ids, c("a", "local")), TRUE)))
})

test_that("homozygous somatic event is carried on both local tumor sides", {
  p <- reconstruction_test_project(); b <- Neoantimon2:::reference_bundle(p)
  a <- reconstruction_test_row("a", "chr1", 4L, "G", "A", side = "LEFT")
  c <- reconstruction_test_row("c", "chr1", 5L, "A", "C", side = "RIGHT")
  both <- reconstruction_test_row("h", "chr1", 22L, "A", "G", both = TRUE)
  contexts <- Filter(function(x) x$transcript_id == "TA.1",
    Neoantimon2:::reconstruction_contexts(p, b, list(a, c, both)))
  expect_length(contexts, 2L)
  expect_true(all(vapply(contexts, function(x) "h" %in% x$somatic_event_ids, TRUE)))
  expect_setequal(vapply(contexts, `[[`, "", "tumor_cds"),
    c("ATGAAACCTGAACCCGGG", "ATGGCACCTGAACCCGGG"))
  expect_false(any(vapply(contexts, function(x)
    identical(x$somatic_event_ids, "h"), TRUE)))
  germ <- reconstruction_test_row("g", "chr1", 7L, "C", "T", "germline", "LEFT")
  with_germline <- Filter(function(x) x$transcript_id == "TA.1" &&
    "a" %in% x$somatic_event_ids, Neoantimon2:::reconstruction_contexts(p, b,
      list(a, c, both, germ)))[[1L]]
  expect_identical(with_germline$phase_state, "RESOLVED_OR_UNAMBIGUOUS")
  expect_true("g" %in% with_germline$germline_event_ids)
})

test_that("BOTH somatic overlays both normal and tumor backgrounds", {
  p <- reconstruction_test_project(); b <- Neoantimon2:::reference_bundle(p)
  contexts <- function(rows) Filter(function(x) x$transcript_id == "TA.1",
    Neoantimon2:::reconstruction_contexts(p, b, rows))
  flip <- function(rows) lapply(rows, function(row) {
    for (i in seq_along(row$event_details_json$phase_assertions)) {
      assertion <- row$event_details_json$phase_assertions[[i]]
      if (!is.null(assertion$local_side)) {
        assertion$local_side <- if (assertion$local_side == "LEFT") "RIGHT" else "LEFT"
        assertion$genotype_state <- assertion$local_side
        row$event_details_json$phase_assertions[[i]] <- assertion
      }
    }
    row
  })
  check_flip <- function(rows, original) {
    flipped <- contexts(flip(rows))
    expect_setequal(vapply(flipped, `[[`, "", "tumor_cds"),
      vapply(original, `[[`, "", "tumor_cds"))
    expect_setequal(vapply(flipped, `[[`, "", "reconstruction_id"),
      vapply(original, `[[`, "", "reconstruction_id"))
  }
  h <- reconstruction_test_row("h", "chr1", 22L, "A", "G", both = TRUE)
  g <- reconstruction_test_row("g", "chr1", 4L, "G", "A", "germline", "LEFT")
  sole <- contexts(list(h, g))
  expect_length(sole, 2L)
  expected <- list(
    list(normal = "ATGGAACCTAAACCCGGG", normal_protein = "MEPKPG",
      tumor = "ATGGAACCTGAACCCGGG", tumor_protein = "MEPEPG"),
    list(normal = "ATGAAACCTAAACCCGGG", normal_protein = "MKPKPG",
      tumor = "ATGAAACCTGAACCCGGG", tumor_protein = "MKPEPG"))
  for (case in expected) {
    x <- Filter(function(x) identical(x$patient_normal_cds, case$normal), sole)
    expect_length(x, 1L)
    expect_identical(x[[1L]]$status, "COMPLETE")
    expect_identical(x[[1L]]$patient_normal_protein, case$normal_protein)
    expect_identical(x[[1L]]$tumor_cds, case$tumor)
    expect_identical(x[[1L]]$tumor_protein, case$tumor_protein)
    expect_identical(x[[1L]]$somatic_event_ids, "h")
  }
  check_flip(list(h, g), sole)
  h_only <- contexts(list(h))
  expect_length(h_only, 1L)
  expect_identical(h_only[[1L]]$tumor_cds, "ATGGAACCTGAACCCGGG")
  expect_identical(h_only[[1L]]$tumor_protein, "MEPEPG")

  a <- reconstruction_test_row("a", "chr1", 4L, "G", "A", side = "LEFT")
  one <- contexts(list(a, h))
  expect_length(one, 2L)
  expect_setequal(vapply(one, `[[`, "", "tumor_cds"),
    c("ATGAAACCTGAACCCGGG", "ATGGAACCTGAACCCGGG"))
  expect_setequal(vapply(one, `[[`, "", "tumor_protein"), c("MKPEPG", "MEPEPG"))
  expect_setequal(vapply(one, function(x) paste(x$somatic_event_ids, collapse = "/"), ""),
    c("a/h", "h"))
  check_flip(list(a, h), one)

  c <- reconstruction_test_row("c", "chr1", 5L, "A", "C", side = "LEFT")
  cis <- contexts(list(a, c, h))
  expect_length(cis, 2L)
  expect_setequal(vapply(cis, `[[`, "", "tumor_cds"),
    c("ATGACACCTGAACCCGGG", "ATGGAACCTGAACCCGGG"))
  expect_setequal(vapply(cis, `[[`, "", "tumor_protein"), c("MTPEPG", "MEPEPG"))
  expect_setequal(vapply(cis, function(x) paste(x$somatic_event_ids, collapse = "/"), ""),
    c("a/c/h", "h"))
  check_flip(list(a, c, h), cis)

  c$event_details_json$phase_assertions[[1L]]$local_side <- "RIGHT"
  c$event_details_json$phase_assertions[[1L]]$genotype_state <- "RIGHT"
  trans <- contexts(list(a, c, h))
  expect_length(trans, 2L)
  expect_setequal(vapply(trans, `[[`, "", "tumor_cds"),
    c("ATGAAACCTGAACCCGGG", "ATGGCACCTGAACCCGGG"))
  expect_setequal(vapply(trans, `[[`, "", "tumor_protein"), c("MKPEPG", "MAPEPG"))
  expect_setequal(vapply(trans, function(x) paste(x$somatic_event_ids, collapse = "/"), ""),
    c("a/h", "c/h"))
  expect_false(any(vapply(trans, function(x) identical(x$somatic_event_ids, "h"), TRUE)))
  check_flip(list(a, c, h), trans)
  c$event_details_json$phase_assertions <- list()
  unknown <- contexts(list(a, c, h))
  expect_length(unknown, 2L)
  expect_setequal(vapply(unknown, function(x)
    paste(x$somatic_event_ids, collapse = "/"), ""), c("a/h", "c/h"))
})

test_that("BOTH uses resolved germline pairs and retains unknown backgrounds", {
  p <- reconstruction_test_project(); b <- Neoantimon2:::reference_bundle(p)
  contexts <- function(rows) Filter(function(x) x$transcript_id == "TA.1",
    Neoantimon2:::reconstruction_contexts(p, b, rows))
  h <- reconstruction_test_row("h", "chr1", 22L, "A", "G", both = TRUE)
  g1 <- reconstruction_test_row("g1", "chr1", 4L, "G", "A", "germline", "LEFT")
  g2 <- reconstruction_test_row("g2", "chr1", 7L, "C", "T", "germline", "LEFT")
  cis <- contexts(list(h, g1, g2))
  expect_length(cis, 2L)
  expect_setequal(vapply(cis, `[[`, "", "tumor_cds"),
    c("ATGGAACCTGAACCCGGG", "ATGAAATCTGAACCCGGG"))
  expect_setequal(vapply(cis, `[[`, "", "tumor_protein"), c("MEPEPG", "MKSEPG"))
  expect_setequal(lapply(cis, `[[`, "germline_event_ids"), list(character(), c("g1", "g2")))
  g2$event_details_json$phase_assertions[[1L]]$local_side <- "RIGHT"
  g2$event_details_json$phase_assertions[[1L]]$genotype_state <- "RIGHT"
  trans <- contexts(list(h, g1, g2))
  expect_length(trans, 2L)
  expect_setequal(vapply(trans, `[[`, "", "tumor_cds"),
    c("ATGAAACCTGAACCCGGG", "ATGGAATCTGAACCCGGG"))
  expect_setequal(lapply(trans, `[[`, "germline_event_ids"), list("g1", "g2"))
  g2$event_details_json$phase_assertions <- list()
  unknown <- contexts(list(h, g1, g2))
  expect_length(unknown, 1L)
  expect_identical(unknown[[1L]]$phase_state, "BACKGROUND_UNKNOWN")
  expect_length(unknown[[1L]]$germline_event_ids, 0L)
  expect_identical(unknown[[1L]]$tumor_cds, "ATGGAACCTGAACCCGGG")

  a <- reconstruction_test_row("a", "chr1", 4L, "G", "A", side = "LEFT")
  g <- reconstruction_test_row("g", "chr1", 7L, "C", "T", "germline", "RIGHT")
  opposite <- contexts(list(a, h, g))
  expect_length(opposite, 2L)
  h_side <- Filter(function(x) identical(x$somatic_event_ids, "h"), opposite)
  expect_length(h_side, 1L)
  expect_identical(h_side[[1L]]$patient_normal_cds, "ATGGAATCTAAACCCGGG")
  expect_identical(h_side[[1L]]$tumor_cds, "ATGGAATCTGAACCCGGG")
  expect_identical(h_side[[1L]]$tumor_protein, "MESEPG")
})

test_that("patient-normal backgrounds stay within their sample", {
  p <- reconstruction_test_project(); b <- Neoantimon2:::reference_bundle(p)
  first <- reconstruction_test_row("s1", "chr1", 4L, "G", "A")
  second <- reconstruction_test_row("s2", "chr1", 4L, "G", "A")
  second$sample_id <- "sample-2"
  germ <- reconstruction_test_row("g2", "chr1", 7L, "C", "T", "germline", both = TRUE)
  germ$sample_id <- "sample-2"
  contexts <- Filter(function(x) x$transcript_id == "TA.1",
    Neoantimon2:::reconstruction_contexts(p, b, list(first, second, germ)))
  expect_length(contexts, 2L)
  by_sample <- stats::setNames(contexts, vapply(contexts, `[[`, "", "sample_id"))
  expect_identical(by_sample[["sample-1"]]$patient_normal_cds, "ATGGAACCTAAACCCGGG")
  expect_length(by_sample[["sample-1"]]$germline_event_ids, 0L)
  expect_identical(by_sample[["sample-2"]]$patient_normal_cds, "ATGGAATCTAAACCCGGG")
  expect_identical(by_sample[["sample-2"]]$germline_event_ids, "g2")
})

test_that("conflicts and T009 boundaries defer without invented protein", {
  p <- reconstruction_test_project(); b <- Neoantimon2:::reference_bundle(p)
  one <- function(...) Filter(function(x) x$transcript_id == "TA.1",
    Neoantimon2:::reconstruction_contexts(p, b, list(...)))[[1L]]
  a <- reconstruction_test_row("a", "chr1", 4L, "G", "A", side = "LEFT")
  conflict <- reconstruction_test_row("b", "chr1", 4L, "G", "C", side = "LEFT")
  expect_identical(one(a, conflict)$reason, "OVERLAPPING_EDIT_CONFLICT")
  compatible <- reconstruction_test_row("m", "chr1", 4L, "GA", "AC", side = "LEFT")
  expect_identical(one(a, compatible)$tumor_protein, "MTPKPG")
  expect_identical(one(reconstruction_test_row("fs", "chr1", 6L, "A", "AT"))$reason,
                   "FRAMESHIFT")
  compensating <- one(reconstruction_test_row("fi", "chr1", 6L, "A", "AT", side = "LEFT"),
    reconstruction_test_row("fd", "chr1", 21L, "TA", "T", side = "LEFT"))
  expect_identical(compensating$reason, "FRAMESHIFT")
  expect_equal(nchar(compensating$tumor_cds), nchar(compensating$reference_cds))
  expect_null(compensating$tumor_protein)
  expect_identical(one(reconstruction_test_row("sg", "chr1", 4L, "G", "T"))$reason,
                   "STOP_GAIN")
  expect_identical(one(reconstruction_test_row("sl", "chr1", 1L, "A", "G"))$reason,
                   "START_LOSS")
  stoploss <- one(reconstruction_test_row("st", "chr1", 31L, "T", "C"))
  expect_identical(stoploss$reason, "STOP_LOSS")
  expect_identical(stoploss$tumor_terminal_codon, "CAA")
  coding <- reconstruction_test_row("coding", "chr1", 4L, "G", "A", side = "LEFT")
  stop_cis <- reconstruction_test_row("stop", "chr1", 31L, "T", "C", side = "LEFT")
  combined_stop <- one(coding, stop_cis)
  expect_identical(combined_stop$reason, "STOP_LOSS")
  expect_identical(combined_stop$tumor_cds, "ATGAAACCTAAACCCGGG")
  expect_identical(combined_stop$tumor_terminal_codon, "CAA")
  sec <- Filter(function(x) x$transcript_id == "TD.1",
    Neoantimon2:::reconstruction_contexts(p, b,
      list(reconstruction_test_row("ss", "chrS", 7L, "C", "T"))))[[1L]]
  expect_identical(sec$reason, "SELENOCYSTEINE")
  expect_identical(sec$reference_cds, "ATGTGACCTAAACCCGGG")
  expect_identical(sec$tumor_cds, "ATGTGATCTAAACCCGGG")
  mito <- Neoantimon2:::reconstruction_contexts(p, b,
    list(reconstruction_test_row("mm", "chrM", 4L, "A", "G")))[[1L]]
  expect_identical(mito$reason, "MITOCHONDRIAL_CODE")
  expect_identical(mito$tumor_cds, "ATGGAA")
  boundary <- one(reconstruction_test_row("bd", "chr1", 8L, "CA", "T"))
  expect_identical(boundary$reason, "CDS_BOUNDARY_OR_SPLICE_INTERSECTION")
  gap <- one(reconstruction_test_row("gap", "chr1", 8L, "C", "CTTT"))
  expect_identical(gap$reason, "CDS_BOUNDARY_OR_SPLICE_INTERSECTION")
})

test_that("incomplete first CDS and inconsistent later phase are not translated", {
  p <- reconstruction_test_project(); b <- Neoantimon2:::reference_bundle(p)
  model <- Neoantimon2:::reconstruction_models(b, "TA.1")[[1L]]
  protein <- Neoantimon2:::reconstruction_fasta_select(
    b$resources$protein_fasta$resolved_path, "PA.1")
  ref <- Neoantimon2:::event_reference_open(b)
  on.exit(Neoantimon2:::event_reference_close(ref))
  model$segments$phase[[1L]] <- "1"
  expect_identical(Neoantimon2:::reconstruction_reference(model, ref, protein,
    character())$reason, "INCOMPLETE_CDS")
  model$segments$phase[[1L]] <- "0"
  model$segments$phase[[2L]] <- "2"
  expect_identical(Neoantimon2:::reconstruction_reference(model, ref, protein,
    character())$reason, "GTF_PHASE_INCONSISTENT")
  model$model_reason <- "EXON_STRUCTURE_UNASSESSABLE"
  expect_identical(Neoantimon2:::reconstruction_reference(model, ref, protein,
    character())$reason, "EXON_STRUCTURE_UNASSESSABLE")
})

test_that("noncoding, zero events, and canonical row objects are preserved", {
  p <- reconstruction_test_project(); b <- Neoantimon2:::reference_bundle(p)
  expect_length(Neoantimon2:::reconstruction_contexts(p, b, list()), 0L)
  outside <- reconstruction_test_row("o", "chr1", 12L, "A", "G")
  expect_length(Neoantimon2:::reconstruction_contexts(p, b, list(outside)), 0L)
  symbolic <- outside; symbolic$event_type <- "BND"
  source <- outside; source$event_type <- "SEQUENCE"; source$contig <- NULL
  expect_length(Neoantimon2:::reconstruction_contexts(p, b, list(symbolic, source)), 0L)
  before <- serialize(outside, NULL)
  Neoantimon2:::reconstruction_contexts(p, b, list(outside))
  expect_identical(serialize(outside, NULL), before)
})

test_that("real runtime writes a checksum-bound internal reconstruction checkpoint", {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=chr1,length=33>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample-1",
    "chr1\t4\t.\tG\tA\t.\tPASS\t.\tGT\t0/1"), vcf)
  p <- reconstruction_test_project(list(list(input_id = "example", input_role = "somatic_vcf",
    sample_id = "sample-1", path = vcf)))
  imported <- Neoantimon2:::event_import_project(p)
  p <- neo_run(p)
  contexts <- readRDS(Neoantimon2:::reconstruction_stage_path(p))
  expect_true(any(vapply(contexts, function(x) identical(x$status, "COMPLETE"), TRUE)))
  after <- Neoantimon2:::read_runtime_snapshot(p)$events
  keys <- c("event_id", "event_key", "canonical_event", "event_role", "contig",
            "pos", "ref", "alt")
  expect_identical(imported[[1L]][keys], after[[1L]][keys])
  cp <- Filter(function(x) x$step_id == "sequence_reconstruction",
    Neoantimon2:::read_manifest(p$path)$checkpoints)[[1L]]
  expect_identical(cp$status, "COMPLETED")
  expect_silent(Neoantimon2:::verify_checkpoint_outputs(cp, p$path))
  expect_length(Neoantimon2:::reconstruction_stage_run(p,
    Neoantimon2:::read_runtime_snapshot(p)$events), length(contexts))
})

test_that("valid no-event production input stores zero reconstruction contexts", {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=chr1,length=33>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample-1",
    "chr1\t4\t.\tG\tA\t.\tPASS\t.\tGT\t0/0"), vcf)
  p <- reconstruction_test_project(list(list(input_id = "zero", input_role = "germline_vcf",
    sample_id = "sample-1", path = vcf)))
  p <- neo_run(p)
  expect_length(readRDS(Neoantimon2:::reconstruction_stage_path(p)), 0L)
  expect_equal(nrow(neo_events(p)), 0L)
  expect_equal(nrow(neo_origins(p)), 0L)
  expect_equal(nrow(neo_results(p)), 0L)
  expect_equal(nrow(neo_multisample(p)), 0L)
})

test_that("coding-free somatic input completes with one event and no contexts", {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=chr1,length=33>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t12\t.\tA\tG\t.\tPASS\t."), vcf)
  p <- reconstruction_test_project(list(list(input_id = "noncoding", input_role = "somatic_vcf",
    sample_id = "sample-1", path = vcf)))
  p <- neo_run(p)
  expect_equal(nrow(neo_events(p)), 1L)
  expect_length(readRDS(Neoantimon2:::reconstruction_stage_path(p)), 0L)
})
