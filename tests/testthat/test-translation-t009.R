t009_contexts <- function(rows) {
  p <- reconstruction_test_project()
  b <- Neoantimon2:::reference_bundle(p)
  old <- Neoantimon2:::reconstruction_contexts(p, b, rows)
  list(old = old, new = Neoantimon2:::translation_contexts(p, b, rows, old),
       project = p, bundle = b)
}

t009_one <- function(rows, transcript) {
  x <- t009_contexts(rows)
  Filter(function(y) identical(y$transcript_id, transcript), x$new)[[1L]]
}

test_that("translation spec maps GENCODE CDS and separate annotated stop", {
  row <- reconstruction_test_row("s", "chrF", 38L, "T", "C")
  x <- t009_one(list(row), "TF.1")
  s <- x$translation$spec
  expect_identical(s$status, "READY")
  expect_identical(s$stop_relation, "AFTER_CDS")
  expect_identical(s$coding_start_tx, 1L)
  expect_identical(s$coding_end_tx, 27L)
  expect_identical(s$annotated_stop_tx, 28L)
  expect_identical(s$transcript_length, 42L)
  expect_identical(s$genetic_code_id, 1L)
  expect_identical(s$start_source, "GTF_START_CODON")
  expect_identical(x$translation$reference_oracle, "EXACT")
  expect_identical(x$translation$reference_protein, "MKPGKPGKP")
  expect_identical(x$translation$reference_stop_aa, 10L)
})

test_that("T008 complete nuclear contexts retain all frozen biological fields", {
  rows <- list(reconstruction_test_row("syn", "chr1", 6L, "A", "G"),
               reconstruction_test_row("mis", "chr1", 4L, "G", "A"),
               reconstruction_test_row("ins", "chr1", 6L, "A", "ACCC"))
  for (row in rows) {
    x <- t009_contexts(list(row))
    for (i in seq_along(x$old)) {
      old <- x$old[[i]]; new <- x$new[[i]]
      expect_identical(new$reconstruction_id, old$reconstruction_id)
      expect_identical(new$reference_cds, old$reference_cds)
      expect_identical(new$patient_normal_cds, old$patient_normal_cds)
      expect_identical(new$tumor_cds, old$tumor_cds)
      expect_identical(new$patient_normal_protein, old$patient_normal_protein)
      expect_identical(new$tumor_protein, old$tumor_protein)
      expect_identical(new$somatic_event_ids, old$somatic_event_ids)
      expect_identical(new$status, old$status)
      expect_identical(new$translation$patient_normal_protein, old$patient_normal_protein)
      expect_identical(new$translation$mutant_protein, old$tumor_protein)
    }
  }
  syn <- t009_one(list(rows[[1L]]), "TA.1")$translation
  expect_identical(syn$protein_effect, "identical")
  expect_length(syn$normal_tumor_difference_intervals, 0L)
  expect_null(syn$first_changed_aa)
  miss <- t009_one(list(rows[[2L]]), "TA.1")$translation
  expect_identical(miss$protein_effect, "substitution")
  expect_identical(miss$mutant_protein, "MKPKPG")
  ins <- t009_one(list(rows[[3L]]), "TA.1")$translation
  expect_identical(ins$protein_effect, "inframe_indel")
  expect_identical(ins$mutant_protein, "MEPPKPG")
})

test_that("plus one and minus one frameshifts terminate at observed downstream stops", {
  plus <- t009_one(list(reconstruction_test_row("plus", "chr1", 6L, "A", "AT")), "TA.1")
  expect_identical(plus$tumor_cds, "ATGGAATCCTAAACCCGGG")
  expect_identical(plus$translation$patient_normal_protein, "MEPKPG")
  expect_identical(plus$translation$mutant_protein, "MES")
  expect_identical(plus$translation$frame_effect, "shifted_unrestored")
  expect_identical(plus$translation$termination_effect, "premature")
  expect_identical(plus$translation$protein_effect, "altered_tail")
  expect_identical(plus$translation$mutant_boundary, "STOP_CODON")
  expect_identical(plus$translation$mutant_stop_aa, 4L)
  minus <- t009_one(list(reconstruction_test_row("minus", "chrD", 5L, "AA", "A")), "TU.1")
  expect_identical(minus$patient_normal_cds, "ATGAAACCCGGGCTAAACCCGGGT")
  expect_identical(minus$tumor_cds, "ATGAACCCGGGCTAAACCCGGGT")
  expect_identical(minus$translation$patient_normal_protein, "MKPGLNPG")
  expect_identical(minus$translation$mutant_protein, "MNPG")
  expect_identical(minus$translation$frame_effect, "shifted_unrestored")
  expect_identical(minus$translation$termination_effect, "premature")
  expect_identical(minus$translation$protein_effect, "altered_tail")
  expect_identical(minus$translation$mutant_stop_aa, 5L)
})

test_that("unrestored frameshift reaches a different stop after the normal stop", {
  x <- t009_one(list(reconstruction_test_row("late", "chrW", 6L, "A", "AC")), "TW.1")
  expect_identical(x$translation$frame_effect, "shifted_unrestored")
  expect_identical(x$translation$termination_effect, "lost_with_recovered_stop")
  expect_identical(x$translation$protein_effect, "altered_tail")
  expect_identical(x$translation$mutant_boundary, "STOP_CODON")
  expect_gt(x$translation$mutant_stop_aa, x$translation$normal_stop_aa)
  expect_false(is.null(x$translation$mutant_stop_origin))
})

test_that("unrestored frameshift stops at known transcript end without padding", {
  x <- t009_one(list(reconstruction_test_row("open", "chrF", 6L, "A", "AC")), "TF.1")
  expect_identical(x$translation$frame_effect, "shifted_unrestored")
  expect_identical(x$translation$termination_effect, "lost_open_end")
  expect_identical(x$translation$protein_effect, "altered_tail")
  expect_identical(x$translation$status, "BOUNDED")
  expect_identical(x$translation$mutant_protein, "MKPRETRETLTRVR")
  expect_identical(x$translation$mutant_boundary, "INCOMPLETE_TERMINAL_CODON")
  expect_null(x$translation$mutant_stop_aa)
  expect_identical(x$translation$residual_terminal_nt, 1L)
})

test_that("cis plus one then minus one restores frame after an exact altered interval", {
  rows <- list(reconstruction_test_row("i", "chrF", 6L, "A", "AC", side = "LEFT"),
               reconstruction_test_row("d", "chrF", 26L, "CC", "C", side = "LEFT"))
  x <- t009_one(rows, "TF.1")
  expect_identical(x$patient_normal_cds, "ATGAAACCCGGGAAACCCGGGAAACCC")
  expect_identical(x$tumor_cds, "ATGAAACCCCGGGAAACCGGGAAACCC")
  expect_identical(x$translation$patient_normal_protein, "MKPGKPGKP")
  expect_identical(x$translation$mutant_protein, "MKPRETGKP")
  expect_identical(x$translation$frame_effect, "shifted_restored")
  expect_identical(x$translation$mutant_frame$final_frame, 0L)
  expect_identical(x$translation$frame_restore_aa, 7L)
  expect_identical(x$translation$frame_restore_nt, 18L)
  expect_identical(x$translation$mutant_frame$shifted_intervals[[1L]]$altered_aa_segment,
                   "PRET")
  expect_identical(substring(x$translation$mutant_protein, 7L), "GKP")
  expect_identical(x$translation$normal_tumor_difference_intervals[[1L]]$normal_sequence,
                   "GKP")
  expect_identical(x$translation$normal_tumor_difference_intervals[[1L]]$mutant_sequence,
                   "RET")
})

test_that("minus one then plus one and three cumulative shifts restore frame", {
  minus_plus <- list(reconstruction_test_row("d", "chrF", 5L, "AA", "A", side = "LEFT"),
    reconstruction_test_row("i", "chrF", 27L, "C", "CC", side = "LEFT"))
  x <- t009_one(minus_plus, "TF.1")
  expect_identical(x$tumor_cds, "ATGAACCCGGGAAACCCCGGGAAACCC")
  expect_identical(x$translation$mutant_protein, "MNPGNPGKP")
  expect_identical(x$translation$frame_effect, "shifted_restored")
  expect_identical(vapply(x$translation$mutant_frame$changes, `[[`, 0L, "frame_after"),
                   c(2L, 0L))
  triple <- list(reconstruction_test_row("a", "chrF", 6L, "A", "AC", side = "LEFT"),
    reconstruction_test_row("b", "chrF", 12L, "G", "GC", side = "LEFT"),
    reconstruction_test_row("c", "chrF", 28L, "C", "CC", side = "LEFT"))
  y <- t009_one(triple, "TF.1")
  expect_identical(y$tumor_cds, "ATGAAACCCCGGGCAAACCCCGGGAAACCC")
  expect_identical(y$translation$mutant_protein, "MKPRANPGKP")
  expect_identical(vapply(y$translation$mutant_frame$changes, `[[`, 0L, "frame_after"),
                   c(1L, 2L, 0L))
  expect_identical(y$translation$frame_effect, "shifted_restored")
  expect_identical(y$translation$frame_restore_aa, 8L)
})

test_that("an upstream stop prevents genomic frame restoration from rescuing protein", {
  rows <- list(reconstruction_test_row("i", "chr1", 6L, "A", "AT", side = "LEFT"),
               reconstruction_test_row("d", "chr1", 21L, "TA", "T", side = "LEFT"))
  x <- t009_one(rows, "TA.1")
  expect_identical(nchar(x$tumor_cds), nchar(x$reference_cds))
  expect_identical(x$translation$mutant_frame$final_frame, 0L)
  expect_identical(x$translation$mutant_protein, "MES")
  expect_identical(x$translation$termination_effect, "premature")
  expect_identical(x$translation$frame_effect, "shifted_unrestored")
  expect_false(x$translation$mutant_frame$shifted_intervals[[1L]]$restoration_reached_by_translation)
})

test_that("minus strand coding frameshift is oriented on the spliced transcript", {
  x <- t009_one(list(reconstruction_test_row("minus", "chrR", 46L, "GT", "G")), "TR.1")
  expect_identical(x$reference_cds, "ATGAAACCCGGGAAACCCGGGAAACCC")
  expect_identical(x$tumor_cds, "ATGAACCCGGGAAACCCGGGAAACCC")
  expect_identical(x$translation$frame_effect, "shifted_unrestored")
  expect_identical(x$translation$mutant_protein, "MNPGNPGNPNPGK")
})

test_that("direct and compound stop gains use the final codon", {
  direct <- t009_one(list(reconstruction_test_row("stop", "chr1", 4L, "G", "T")), "TA.1")
  expect_identical(direct$translation$termination_effect, "premature")
  expect_identical(direct$translation$frame_effect, "unchanged")
  expect_identical(direct$translation$mutant_protein, "M")
  expect_identical(direct$translation$mutant_stop_aa, 2L)
  expect_identical(direct$translation$protein_effect, "truncated")
  expect_identical(direct$translation$normal_tumor_difference_intervals[[1L]]$mutant_sequence, "")
  two <- list(reconstruction_test_row("a", "chrF", 10L, "G", "T", side = "LEFT"),
              reconstruction_test_row("b", "chrF", 12L, "G", "A", side = "LEFT"))
  a <- t009_one(list(two[[1L]]), "TF.1")
  b <- t009_one(list(two[[2L]]), "TF.1")
  together <- t009_one(two, "TF.1")
  expect_identical(a$translation$termination_effect, "unchanged")
  expect_identical(b$translation$termination_effect, "unchanged")
  expect_identical(together$translation$termination_effect, "premature")
  expect_identical(together$translation$mutant_protein, "MKP")
  expect_identical(together$translation$mutant_stop_aa, 4L)
  mnv <- t009_one(list(reconstruction_test_row("mnv", "chrF", 10L, "GGG", "TGA")), "TF.1")
  expect_identical(mnv$translation$mutant_protein, together$translation$mutant_protein)
})

test_that("a phased second SNV rescues an intermediate stop", {
  a <- reconstruction_test_row("stop", "chr1", 4L, "G", "T", side = "LEFT")
  b <- reconstruction_test_row("rescue", "chr1", 5L, "A", "C", side = "LEFT")
  expect_identical(t009_one(list(a), "TA.1")$translation$mutant_protein, "M")
  combined <- t009_one(list(a, b), "TA.1")
  expect_identical(combined$tumor_cds, "ATGTCACCTAAACCCGGG")
  expect_identical(combined$translation$mutant_protein, "MSPKPG")
  expect_identical(combined$translation$termination_effect, "unchanged")
})

test_that("stop loss uses spliced plus and minus transcript extension", {
  plus <- t009_one(list(reconstruction_test_row("plus", "chrF", 38L, "T", "C")), "TF.1")
  minus <- t009_one(list(reconstruction_test_row("minus", "chrR", 15L, "A", "G")), "TR.1")
  for (x in list(plus, minus)) {
    expect_identical(x$translation$patient_normal_protein, "MKPGKPGKP")
    expect_identical(x$translation$mutant_protein, "MKPGKPGKPQPG")
    expect_identical(x$translation$extension_sequence, "QPG")
    expect_identical(x$translation$termination_effect, "lost_with_recovered_stop")
    expect_identical(x$translation$frame_effect, "unchanged")
    expect_identical(x$translation$protein_effect, "extended")
    expect_identical(x$translation$mutant_stop_aa, 13L)
    expect_identical(x$translation$normal_stop_origin, 28:30)
    expect_identical(x$translation$mutant_stop_origin, 37:39)
    expect_identical(x$translation$spec$stop_relation, "AFTER_CDS")
  }
  expect_identical(plus$translation$mutant_boundary, "STOP_CODON")
  expect_identical(minus$translation$mutant_boundary, "STOP_CODON")
  expect_identical(minus$strand, "-")
  # chrR's intron contains only T; its genomic sequence cannot encode QPG here.
  expect_identical(minus$translation$spec$transcript_length, 42L)
})

test_that("stop loss without another known stop remains bounded", {
  x <- t009_one(list(reconstruction_test_row("open", "chrO", 28L, "T", "C")), "TO.1")
  expect_identical(x$translation$mutant_protein, "MKPGKPGKPQPGK")
  expect_identical(x$translation$extension_sequence, "QPGK")
  expect_identical(x$translation$termination_effect, "lost_open_end")
  expect_identical(x$translation$status, "BOUNDED")
  expect_null(x$translation$mutant_stop_aa)
  expect_identical(x$translation$mutant_boundary, "TRANSCRIPT_END_NO_STOP")
})

test_that("extension variants outside T008 context are recorded, not ignored", {
  stop <- reconstruction_test_row("stop", "chrF", 38L, "T", "C", side = "LEFT")
  downstream <- reconstruction_test_row("extension", "chrF", 41L, "C", "T",
                                        side = "LEFT")
  same <- t009_one(list(stop, downstream), "TF.1")
  expect_identical(same$translation$status, "UNRESOLVED")
  expect_identical(same$translation$reason,
                   "EXTENSION_EVENT_OUTSIDE_T008_HAPLOTYPE")
  expect_identical(same$translation$extension_event_ids, "extension")
  expect_identical(same$translation$extension_event_statuses[["extension"]],
                   "KNOWN_SAME_HAPLOTYPE")
  downstream$event_details_json$phase_assertions <- list()
  unknown <- t009_one(list(stop, downstream), "TF.1")
  expect_identical(unknown$translation$reason, "EXTENSION_PHASE_UNRESOLVED")
  downstream$event_details_json$phase_assertions <-
    reconstruction_test_row("extension", "chrF", 41L, "C", "T",
      side = "RIGHT")$event_details_json$phase_assertions
  opposite <- t009_one(list(stop, downstream), "TF.1")
  expect_identical(opposite$translation$status, "EVALUABLE")
  expect_identical(opposite$translation$mutant_protein, "MKPGKPGKPQPG")
})

test_that("annotated start loss cannot use a downstream ATG", {
  x <- t009_one(list(reconstruction_test_row("lost", "chrA", 1L, "A", "G")), "TV.1")
  expect_identical(x$tumor_cds, "GTGAAAATGCCC")
  expect_identical(x$translation$spec$start_codon, "ATG")
  expect_identical(x$translation$start_effect, "lost")
  expect_null(x$translation$mutant_protein)
  expect_identical(x$translation$protein_effect, "unavailable")
  expect_identical(substring(x$tumor_cds, 7L, 9L), "ATG")
})

test_that("mitochondrial codons, non-ATG start and partial stop are catalog validated", {
  x <- t009_one(list(reconstruction_test_row("codon", "chrM", 13L, "T", "C")), "TM2.1")
  expect_identical(x$translation$translation_system, "vertebrate_mitochondrial")
  expect_identical(x$translation$genetic_code_id, 2L)
  expect_identical(x$translation$spec$start_codon, "ATA")
  expect_identical(x$translation$reference_protein, "MW")
  expect_identical(x$translation$mutant_protein, "MR")
  expect_identical(x$translation$reference_oracle, "EXACT")
  partial <- t009_one(list(reconstruction_test_row("partial", "chrM", 22L, "T", "C")), "TM3.1")
  expect_identical(partial$translation$reference_protein, "MW")
  expect_identical(partial$translation$spec$terminal_completion,
                   "CATALOG_VALIDATED_PARTIAL_TAA")
  expect_identical(partial$translation$mutant_boundary, "CATALOG_COMPLETED_PARTIAL_STOP")
  expect_identical(partial$translation$mutant_protein, "MR")
  expect_identical(partial$translation$reference_oracle, "EXACT")
  start <- t009_one(list(reconstruction_test_row("start", "chrM", 19L, "A", "G")), "TM3.1")
  expect_identical(start$translation$start_effect, "lost")
  expect_null(start$translation$mutant_protein)
})

test_that("Sec applies only at exact annotated UGA codon and restored frame", {
  sec <- t009_one(list(reconstruction_test_row("syn", "chrS2", 4L, "A", "G")), "TS.1")
  expect_identical(sec$translation$reference_protein, "MKPGKPGUP")
  expect_identical(sec$translation$mutant_protein, "MEPGKPGUP")
  expect_identical(sec$translation$sec_sites_used, "GENCODE_SEC:TS.1:8")
  changed <- t009_one(list(reconstruction_test_row("changed", "chrS2", 22L, "T", "C")), "TS.1")
  expect_identical(changed$translation$mutant_protein, "MKPGKPGRP")
  expect_length(changed$translation$sec_sites_used, 0L)
  shifted <- t009_one(list(reconstruction_test_row("shift", "chrS2", 6L, "A", "AC")), "TS.1")
  expect_length(shifted$translation$sec_sites_used, 0L)
  rows <- list(reconstruction_test_row("shift", "chrS2", 6L, "A", "AC", side = "LEFT"),
               reconstruction_test_row("restore", "chrS2", 16L, "CC", "C", side = "LEFT"))
  restored <- t009_one(rows, "TS.1")
  expect_identical(restored$translation$frame_effect, "shifted_restored")
  expect_identical(restored$translation$sec_sites_used, "GENCODE_SEC:TS.1:8")
  expect_true(grepl("U", restored$translation$mutant_protein, fixed = TRUE))
  code <- Neoantimon2:::translation_codons("nuclear_standard")
  expect_identical(unname(code[["TGA"]]), "*")
})

test_that("partial boundaries, absent sequence and oracle mismatch stay explicit", {
  rows <- list(reconstruction_test_row("m", "chrF", 4L, "A", "G"))
  x <- t009_contexts(rows); context <- Filter(function(y) y$transcript_id == "TF.1", x$old)[[1L]]
  model <- Neoantimon2:::reconstruction_models(x$bundle, "TF.1")[[1L]]
  features <- Neoantimon2:::reference_transcript_features(x$bundle, "TF.1")
  proteins <- list("PF.1" = "MKPGKPGKP")
  transcripts <- list("TF.1" = context$reference_transcript)
  five <- model; five$segments$phase[[1L]] <- "1"
  spec5 <- Neoantimon2:::translation_specification(five, context, x$bundle,
    proteins, transcripts, features)
  expect_true(spec5$partial_5)
  expect_identical(Neoantimon2:::translation_reference_oracle(spec5)$status, "PARTIAL")
  three <- model; three$stop_segments <- three$stop_segments[0L, , drop = FALSE]
  spec3 <- Neoantimon2:::translation_specification(three, context, x$bundle,
    proteins, transcripts, features)
  expect_true(spec3$partial_3)
  expect_identical(Neoantimon2:::translation_reference_oracle(spec3)$status, "PARTIAL")
  wrong <- Neoantimon2:::translation_specification(model, context, x$bundle,
    list("PF.1" = "WRONG"), transcripts, features)
  expect_error(Neoantimon2:::translation_reference_oracle(wrong, strict = TRUE),
               "REFERENCE_PROTEIN_MISMATCH")
  missing <- context; missing$tumor_transcript <- NULL
  unresolved <- Neoantimon2:::translation_context(missing, model, rows, x$bundle,
    proteins, transcripts, features)
  expect_identical(unresolved$translation$status, "UNRESOLVED")
  expect_identical(unresolved$translation$reason, "T008_TUMOR_TRANSCRIPT_DISAGREEMENT")
})

test_that("two separated substitutions remain two mapped protein intervals", {
  rows <- list(reconstruction_test_row("a", "chrF", 4L, "A", "G", side = "LEFT"),
               reconstruction_test_row("b", "chrF", 29L, "G", "A", side = "LEFT"))
  x <- t009_one(rows, "TF.1")
  expect_identical(x$translation$patient_normal_protein, "MKPGKPGKP")
  expect_identical(x$translation$mutant_protein, "MEPGKPRKP")
  expect_length(x$translation$normal_tumor_difference_intervals, 2L)
  expect_identical(vapply(x$translation$normal_tumor_difference_intervals,
    `[[`, 0L, "normal_aa_start"), c(2L, 7L))
  expect_identical(x$translation$first_changed_aa, 2L)
  expect_identical(x$translation$last_changed_aa, 7L)
})

test_that("patient-normal, isoforms, phase orientation and source order are retained", {
  som <- reconstruction_test_row("s", "chr1", 22L, "A", "G", both = TRUE)
  germ <- reconstruction_test_row("g", "chr1", 4L, "G", "A", "germline", "LEFT")
  x <- t009_contexts(list(som, germ))$new
  ta <- Filter(function(z) z$transcript_id == "TA.1", x)
  expect_length(ta, 2L)
  expect_setequal(vapply(ta, function(z) z$translation$patient_normal_protein, ""),
                  c("MEPKPG", "MKPKPG"))
  expect_setequal(vapply(ta, function(z) z$translation$mutant_protein, ""),
                  c("MEPEPG", "MKPEPG"))
  expect_true(any(vapply(x, function(z) z$transcript_id == "TB.1", TRUE)))
  reverse <- t009_contexts(list(germ, som))$new
  key <- function(xs) sort(vapply(xs, function(z) paste(z$reconstruction_id,
    z$translation$translation_id, z$translation$mutant_protein, sep = "/"), ""))
  expect_identical(key(x), key(reverse))
  flipped <- germ
  flipped$event_details_json$phase_assertions[[1L]]$local_side <- "RIGHT"
  flipped$event_details_json$phase_assertions[[1L]]$genotype_state <- "RIGHT"
  expect_identical(key(x), key(t009_contexts(list(som, flipped))$new))
  a <- reconstruction_test_row("a", "chr1", 4L, "G", "A", side = "LEFT")
  b <- reconstruction_test_row("b", "chr1", 5L, "A", "C", side = "RIGHT")
  trans <- t009_contexts(list(a, b))$new
  expect_length(Filter(function(z) z$transcript_id == "TA.1", trans), 2L)
  expect_false(any(vapply(trans, function(z) identical(z$somatic_event_ids, c("a", "b")), TRUE)))
})

test_that("no event and no coding context are successful empty translations", {
  p <- reconstruction_test_project(); b <- Neoantimon2:::reference_bundle(p)
  expect_length(Neoantimon2:::translation_contexts(p, b, list(), list()), 0L)
  outside <- reconstruction_test_row("outside", "chr1", 12L, "A", "G")
  expect_length(t009_contexts(list(outside))$new, 0L)
})

test_that("premature-stop geometry remains structural with no NMD prediction", {
  upstream <- t009_one(list(reconstruction_test_row("stop", "chrF", 10L, "GGG", "TGA")), "TF.1")
  g <- upstream$translation$termination_geometry
  expect_identical(g$status, "STRUCTURAL_ONLY_NO_BINARY_POLICY")
  expect_identical(g$ptc_tx, 10L)
  expect_identical(g$final_junction_tx, 15L)
  expect_identical(g$distance_to_final_junction_nt, 5L)
  expect_false(g$ptc_in_last_exon)
  expect_false(g$intronless)
  expect_null(g$nmd_prediction)
  last <- t009_one(list(reconstruction_test_row("stop", "chrF", 29L, "GGG", "TGA")), "TF.1")
  expect_true(last$translation$termination_geometry$ptc_in_last_exon)
  single <- t009_one(list(reconstruction_test_row("stop", "chrD", 13L, "CTA", "TAA")), "TU.1")
  expect_true(single$translation$termination_geometry$intronless)
  expect_null(single$translation$termination_geometry$final_junction_tx)
  spec <- upstream$translation$spec; spec$exon_segments <- NULL
  insufficient <- Neoantimon2:::translation_termination_geometry(spec,
    list(stop_nt = 10L), "premature")
  expect_identical(insufficient$status, "JUNCTION_INFORMATION_INSUFFICIENT")
  expect_null(insufficient$nmd_prediction)
})

test_that("runtime stores translation checkpoint and does not generate peptides", {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=chr1,length=33>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t4\t.\tG\tA\t.\tPASS\t."), vcf)
  p <- reconstruction_test_project(list(list(input_id = "t009", input_role = "somatic_vcf",
    sample_id = "sample-1", path = vcf)))
  p <- neo_run(p)
  artifact <- readRDS(Neoantimon2:::translation_stage_path(p))
  expect_true(length(artifact) >= 1L)
  cp <- Filter(function(z) z$step_id == "translation_resolution",
    Neoantimon2:::read_manifest(p$path)$checkpoints)[[1L]]
  expect_identical(cp$status, "COMPLETED")
  expect_silent(Neoantimon2:::verify_checkpoint_outputs(cp, p$path))
  expect_equal(nrow(neo_origins(p)), 0L)
  expect_equal(nrow(neo_results(p)), 0L)
  expect_identical(length(Neoantimon2:::translation_stage_run(p,
    Neoantimon2:::read_runtime_snapshot(p)$events,
    readRDS(Neoantimon2:::reconstruction_stage_path(p)))), length(artifact))
})
