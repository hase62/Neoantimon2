source_fixture_input <- function(name, role) list(input_id = name,
  input_role = role, sample_id = "sample-1",
  path = normalizePath(testthat::test_path("..", "fixtures", "source", name),
    winslash = "/"), checksum = NULL, version = "fixture-v1")

source_fixture_project <- function(inputs, short = FALSE, hla = NULL) {
  base <- reconstruction_test_project()
  cfg <- base$config
  if (short) {
    values <- unclass(cfg)
    values$peptides <- list(class_i_lengths = list(6L),
      class_ii_lengths = list(13L))
    cfg <- do.call(neo_config, c(values, list(env = character())))
  }
  neo_project(samples = data.frame(sample_id = "sample-1",
    patient_id = "patient-1"), inputs = inputs, config = cfg, hla = hla)
}

test_that("STAR-Fusion and Arriba preserve observations and share breakpoint identity", {
  inputs <- list(source_fixture_input("star-fusion.tsv", "star_fusion"),
    source_fixture_input("arriba.tsv", "arriba"))
  p <- source_fixture_project(inputs)
  x <- Neoantimon2:::source_extension_construct(p)
  expect_length(x$observations, 2L)
  expect_length(x$canonical_events, 1L)
  expect_identical(x$observations[[1L]]$source_event_id,
    x$observations[[2L]]$source_event_id)
  expect_length(x$canonical_events[[1L]]$observation_ids, 2L)
  expect_true(all(vapply(x$observations, function(z) length(z$raw_fields) > 10L, TRUE)))
  expect_true(length(x$hypotheses) >= 1L)
  reversed <- source_fixture_project(rev(inputs))
  y <- Neoantimon2:::source_extension_construct(reversed)
  expect_identical(x$canonical_events, y$canonical_events)
  expect_identical(vapply(x$hypotheses, `[[`, "", "transcript_structure_id"),
    vapply(y$hypotheses, `[[`, "", "transcript_structure_id"))
  alias_path <- tempfile("arriba-alias-", fileext = ".tsv")
  alias_lines <- readLines(inputs[[2L]]$path)
  alias_lines[[2L]] <- sub("^G1\\t", "G1_ALIAS\t", alias_lines[[2L]])
  writeLines(alias_lines, alias_path)
  alias <- inputs[[2L]]
  alias$input_id <- "arriba-alias"
  alias$path <- alias_path
  z <- Neoantimon2:::source_extension_construct(source_fixture_project(
    list(inputs[[1L]], alias)))
  expect_length(z$canonical_events, 1L)
})

test_that("fusion adapters fail unknown required headers", {
  path <- tempfile("unsupported-fusion-", fileext = ".tsv")
  writeLines(c("FusionName\tGene1\tGene2", "X--Y\tX\tY"), path)
  input <- list(input_id = "bad", input_role = "star_fusion",
    sample_id = "sample-1", resolved_path = path, version = "fixture")
  expect_error(Neoantimon2:::source_fusion_rows(input, "GRCh38.p14"),
    "Unsupported STAR-Fusion header")
})

test_that("fusion and explicit splice structures reach candidate results", {
  hla <- list("patient-1" = "HLA-A*02:01")
  fusion <- source_fixture_project(list(source_fixture_input("star-fusion.tsv",
    "star_fusion")), hla = hla)
  fusion <- neo_project(samples = fusion$samples, hla = hla, config = fusion$config,
    inputs = list(source_fixture_input("star-fusion.tsv", "star_fusion")),
    path = file.path(fusion$config$paths$run_root, "fusion-source-e2e"))
  fusion <- neo_run(fusion)
  rows <- neo_results(fusion)
  expect_gt(nrow(rows), 0L)
  expect_true(all(rows$source_type == "fusion"))
  expect_true(all(rows$fusion_partners == "G1--GS"))
  expect_true(all(grepl("chr1:8", rows$fusion_breakpoints, fixed = TRUE)))
  expect_true(all(rows$grading_status == "NOT_GRADED"))
  expect_true(length(neo_explain(fusion, candidate_id = rows$candidate_id[[1L]])$
    candidates[[1L]]$source_extension$provider_observations) >= 1L)
  source <- Neoantimon2:::source_read(fusion)
  prepared <- Neoantimon2:::predictor_prepare(source$candidates,
    "fixture", "mhc1_binding", predictor_test_capability("fixture"))
  peptide <- source$candidates[[1L]]$ligand_sequence
  allele <- names(prepared$batches)[[1L]]
  raw <- paste("Pos MHC Peptide Of Core Score_EL %Rank_EL Score_BA Affinity(nM) %Rank_BA",
    paste("1", allele, peptide, "1", peptide, "0.5", "1.2", "0.4", "123", "2.0"),
    sep = "\n")
  prediction <- Neoantimon2:::predictor_run(prepared,
    Neoantimon2:::predictor_artifact_root(fusion),
    fixture_raw = stats::setNames(list(raw), allele))
  expect_true(length(prediction$evidence) > 0L)
  expect_true(any(neo_results(fusion)$prediction_status == "SUCCESS"))
  explained <- neo_explain(fusion, candidate_id = prediction$evidence[[1L]]$candidate_id)
  expect_true(length(explained$candidates[[1L]]$source_extension$
    prediction_evidence$evidence) > 0L)

  splice <- source_fixture_project(list(source_fixture_input("splice-structure.tsv",
    "splice_structure")), short = TRUE, hla = hla)
  splice <- neo_project(samples = splice$samples, hla = hla, config = splice$config,
    inputs = list(source_fixture_input("splice-structure.tsv", "splice_structure")),
    path = file.path(splice$config$paths$run_root, "splice-source-e2e"))
  splice <- neo_run(splice)
  result <- neo_results(splice)
  expect_gt(nrow(result), 0L)
  expect_true(all(result$source_type == "novel_junction"))
  expect_true(all(result$splice_donor == "7"))
  expect_true(all(result$grading_status == "NOT_GRADED"))
  splice_source <- Neoantimon2:::source_read(splice)
  splice_prepared <- Neoantimon2:::predictor_prepare(splice_source$candidates,
    "fixture", "mhc1_binding", predictor_test_capability("fixture"))
  splice_allele <- names(splice_prepared$batches)[[1L]]
  splice_peptide <- splice_source$candidates[[1L]]$ligand_sequence
  splice_raw <- paste("Pos MHC Peptide Of Core Score_EL %Rank_EL Score_BA Affinity(nM) %Rank_BA",
    paste("1", splice_allele, splice_peptide, "1", splice_peptide,
      "0.5", "1.2", "0.4", "123", "2.0"), sep = "\n")
  splice_prediction <- Neoantimon2:::predictor_run(splice_prepared,
    Neoantimon2:::predictor_artifact_root(splice),
    fixture_raw = stats::setNames(list(splice_raw), splice_allele))
  expect_true(length(splice_prediction$evidence) > 0L)
  expect_identical(neo_results(splice)$prediction_status[[1L]], "SUCCESS")
  expect_true(length(neo_explain(splice, candidate_id = result$candidate_id[[1L]])$
    candidates[[1L]]$source_extension$prediction_evidence$evidence) > 0L)
})

test_that("SpliceAI prediction alone never reconstructs sequence", {
  p <- source_fixture_project(list(source_fixture_input("spliceai.vcf",
    "spliceai_vcf")))
  x <- Neoantimon2:::source_extension_construct(p)
  expect_length(x$splice_predictions, 1L)
  expect_false(x$splice_predictions[[1L]]$sequence_generating)
  expect_length(x$hypotheses, 0L)
  expect_length(x$origins, 0L)
  expect_length(x$candidates, 0L)
  bound <- neo_project(samples = p$samples, config = p$config,
    inputs = list(source_fixture_input("spliceai.vcf", "spliceai_vcf")),
    path = file.path(p$config$paths$run_root, "spliceai-evidence-only"))
  bound <- neo_run(bound)
  expect_equal(nrow(neo_results(bound)), 0L)
  expect_length(Neoantimon2:::source_read(bound)$splice_predictions, 1L)
})

test_that("STAR junction observations remain evidence only", {
  p <- source_fixture_project(list(source_fixture_input("SJ.out.tab",
    "star_junction")))
  x <- Neoantimon2:::source_extension_construct(p)
  expect_length(x$junction_observations, 1L)
  expect_identical(x$junction_observations[[1L]]$intron_first_base, 9L)
  expect_false(x$junction_observations[[1L]]$sequence_generating)
  expect_length(x$candidates, 0L)
})

test_that("fusion mapping follows coordinates and sequence, not caller frame", {
  p <- source_fixture_project(list(source_fixture_input("star-fusion-cases.tsv",
    "star_fusion")))
  x <- Neoantimon2:::source_extension_construct(p)
  expect_length(x$observations, 6L)
  no_model <- Filter(function(e) e$structure$genes[[2L]]$id == "GU.9",
    x$canonical_events)
  expect_identical(no_model[[1L]]$mapping_status, "NO_COMPATIBLE_TRANSCRIPT")
  inside <- Filter(function(e) e$structure$left$position == 7L &&
      e$structure$right$contig == "chrS",
    x$canonical_events)[[1L]]
  hypothesis <- Filter(function(h) h$source_event_id == inside$source_event_id,
    x$hypotheses)[[1L]]
  expect_identical(hypothesis$orf$frame_state, "out_of_frame")
  provider <- Filter(function(o) o$source_event_id == inside$source_event_id,
    x$observations)[[1L]]
  expect_identical(provider$provider_frame, "INFRAME")
  expect_true(length(hypothesis$orf$novel_aa) > 0L)
  expect_true(length(x$peptides) > 0L)
  intra <- Filter(function(e) e$structure$left$contig ==
    e$structure$right$contig, x$canonical_events)
  expect_true(any(vapply(intra, function(e) e$mapping_status == "MAPPED", TRUE)))
  ambiguous <- Filter(function(e) e$structure$left$position == 5L,
    x$canonical_events)[[1L]]
  expect_identical(ambiguous$transcript_hypothesis_status, "AMBIGUOUS")
  hypotheses <- Filter(function(h) h$source_event_id == ambiguous$source_event_id,
    x$hypotheses)
  expect_length(hypotheses, 2L)
  expect_false(identical(hypotheses[[1L]]$transcript_structure_id,
    hypotheses[[2L]]$transcript_structure_id))
  expect_false(identical(hypotheses[[1L]]$altered_transcript_id,
    hypotheses[[2L]]$altered_transcript_id))
  for (field in c("orf_hypothesis_id", "protein_origin_id"))
    expect_false(identical(hypotheses[[1L]]$orf[[field]], hypotheses[[2L]]$orf[[field]]))
  shared <- Filter(function(p) p$source_event_id == ambiguous$source_event_id,
    x$peptides)
  expect_gte(length(shared), 2L)
  expect_identical(shared[[1L]]$peptide_key, shared[[2L]]$peptide_key)
  expect_false(identical(shared[[1L]]$origin_id, shared[[2L]]$origin_id))
})

test_that("fusion without retained annotated start has no ORF or peptide", {
  p <- source_fixture_project(list(source_fixture_input("star-fusion.tsv",
    "star_fusion")))
  bundle <- Neoantimon2:::reference_bundle(p, retain_indexes = TRUE)
  reference <- Neoantimon2:::event_reference_open(bundle)
  on.exit(Neoantimon2:::event_reference_close(reference), add = TRUE)
  models <- Neoantimon2:::source_model_catalog(bundle,
    Neoantimon2:::source_fusion_rows(p$inputs[[1L]], "GRCh38.p14"))
  sequences <- lapply(models, Neoantimon2:::source_model_sequence,
    bundle = bundle, reference = reference)
  left <- sequences[["TA.1"]]
  right <- sequences[["TD.1"]]
  left$coding_start <- 5L
  left$coding_start_coordinate <- left$coordinates[[5L]]
  event <- Neoantimon2:::source_fusion_rows(p$inputs[[1L]],
    "GRCh38.p14")[[1L]]
  event$structure$left$position <- 3L
  altered <- Neoantimon2:::source_fusion_structure(event, left, right, bundle)
  expect_true(nzchar(altered$sequence))
  expect_identical(altered$status, "ORF_UNAVAILABLE_START_LOST")
  orf <- Neoantimon2:::source_orf(event, altered, "lost-start")
  expect_identical(orf$status, "UNAVAILABLE")
  expect_length(Neoantimon2:::source_peptide_windows(orf,
    Neoantimon2:::peptide_policy()), 0L)
})

test_that("explicit splice coordinates control donor, acceptor and intron", {
  p <- source_fixture_project(list(source_fixture_input("splice-cases.tsv",
    "splice_structure")))
  x <- Neoantimon2:::source_extension_construct(p)
  by_type <- function(type) {
    e <- Filter(function(z) z$source_type == type, x$observations)[[1L]]
    Filter(function(z) z$source_event_id == e$source_event_id,
      x$hypotheses)[[1L]]
  }
  expect_identical(by_type("alternative_donor")$orf$frame_state, "out_of_frame")
  expect_identical(by_type("alternative_acceptor")$orf$frame_state, "out_of_frame")
  expect_identical(by_type("retained_intron")$orf$frame_state, "in_frame")
  retained <- Filter(function(z) z$structure$source_type == "retained_intron" &&
    z$structure$strand == "-", x$canonical_events)[[1L]]
  retained_orf <- Filter(function(z) z$source_event_id == retained$source_event_id,
    x$hypotheses)[[1L]]$orf
  expect_identical(retained_orf$frame_state, "out_of_frame")
  minus <- Filter(function(z) z$structure$strand == "-", x$canonical_events)
  expect_identical(minus[[1L]]$mapping_status, "MAPPED")
  absent <- Filter(function(z) z$structure$gene_id == "GU.9", x$canonical_events)
  expect_identical(absent[[1L]]$mapping_status, "NO_COMPATIBLE_TRANSCRIPT")
  ambiguous <- Filter(function(z) z$structure$donor == 4L,
    x$canonical_events)
  expect_identical(ambiguous[[1L]]$transcript_hypothesis_status, "AMBIGUOUS")
})

test_that("exon skip uses annotated flanks and translated frame", {
  p <- source_fixture_project(list(source_fixture_input("splice-structure.tsv",
    "splice_structure")))
  bundle <- Neoantimon2:::reference_bundle(p)
  reference <- Neoantimon2:::event_reference_open(bundle)
  on.exit(Neoantimon2:::event_reference_close(reference), add = TRUE)
  make_model <- function(middle_end) {
    exons <- data.frame(start = c(1, 16, 26), end = c(9, middle_end, 52))
    coordinates <- unlist(Map(seq.int, exons$start, exons$end), use.names = FALSE)
    sequence <- paste(vapply(seq_len(nrow(exons)), function(i)
      Neoantimon2:::event_reference_fetch(reference, "chrF",
        exons$start[[i]], exons$end[[i]]), ""), collapse = "")
    list(model = list(gene_id = "GF.1", contig = "chrF", strand = "+",
      transcript_id = paste0("manual-", middle_end), exons = exons),
      sequence = sequence, coordinates = coordinates, coding_start = 1L,
      coding_start_coordinate = 1L, sec_sites = list())
  }
  event <- list(source_type = "exon_skip", source_event_id = "manual-skip",
    structure = list(source_type = "exon_skip", gene_id = "GF.1",
      contig = "chrF", strand = "+", donor = 9L, acceptor = 26L,
      transcript_id = NULL))
  in_frame <- Neoantimon2:::source_splice_structure(event, make_model(24), reference)
  shifted <- Neoantimon2:::source_splice_structure(event, make_model(23), reference)
  in_orf <- Neoantimon2:::source_orf(event, in_frame, "in")
  expect_identical(in_orf$frame_state,
    "in_frame")
  expect_identical(Neoantimon2:::source_orf(event, shifted, "shift")$frame_state,
    "out_of_frame")
  expect_true(length(in_frame$coordinates) < length(make_model(24)$coordinates))
  expect_length(in_frame$novel_nt_intervals, 0L)
  expect_length(in_orf$novel_aa, 0L)
  expect_length(in_orf$inserted_sequence_aa, 0L)
  skip_peptides <- Neoantimon2:::source_peptide_windows(in_orf,
    Neoantimon2:::peptide_policy(list(class_i_lengths = list(3L),
      class_ii_lengths = list(4L))))
  expect_true(length(skip_peptides) > 0L)
  expect_true(all(vapply(skip_peptides, function(x)
    x$junction_crossing || x$novel_residue, TRUE)))
  expect_true(all(vapply(skip_peptides, `[[`, TRUE, "junction_crossing")))
  expect_false(any(vapply(skip_peptides, function(x) x$start_aa >= 4L, TRUE)))
  lost <- make_model(24)
  lost$coding_start <- match(16L, lost$coordinates)
  lost$coding_start_coordinate <- 16L
  lost_structure <- Neoantimon2:::source_splice_structure(event, lost, reference)
  expect_identical(lost_structure$status, "ORF_UNAVAILABLE_START_LOST")
  expect_length(Neoantimon2:::source_orf(event, lost_structure, "lost")$peptides, 0L)
})

test_that("retained intron translates real stop and shifted open tail", {
  make_case <- function(intron) {
    upstream <- paste0("ATG", paste(rep("AAA", 8L), collapse = ""))
    downstream <- paste(rep("CCC", 9L), collapse = "")
    genome <- paste0(upstream, intron, downstream)
    path <- tempfile("mini-intron-", fileext = ".fa")
    writeLines(c(">chrMini", genome), path)
    reference <- list(con = file(path, "rb"), fai = data.frame(
      contig = "chrMini", length = nchar(genome), offset = 9,
      line_bases = nchar(genome), line_width = nchar(genome) + 1L))
    exons <- data.frame(start = c(1L, 28L + nchar(intron)),
      end = c(27L, 27L + nchar(intron) + 27L))
    model <- list(model = list(gene_id = "GMINI", contig = "chrMini",
      strand = "+", transcript_id = "TXMINI", exons = exons),
      sequence = paste0(upstream, downstream),
      coordinates = c(1:27, seq.int(exons$start[[2L]], exons$end[[2L]])),
      coding_start = 1L, coding_start_coordinate = 1L, sec_sites = list())
    event <- list(source_type = "retained_intron", source_event_id = "mini",
      structure = list(source_type = "retained_intron", gene_id = "GMINI",
        contig = "chrMini", strand = "+", donor = 27L,
        acceptor = exons$start[[2L]], transcript_id = "TXMINI",
        retained_start = 28L, retained_end = exons$start[[2L]] - 1L))
    altered <- Neoantimon2:::source_splice_structure(event, model, reference)
    close(reference$con)
    Neoantimon2:::source_orf(event, altered, paste0("mini-", intron))
  }
  stopped <- make_case("TAA")
  expect_identical(stopped$termination, "STOP_CODON")
  expect_identical(stopped$frame_state, "in_frame")
  expect_equal(nchar(stopped$protein), 9L)
  expect_length(stopped$inserted_sequence_aa, 0L)
  expect_length(stopped$novel_aa, 0L)
  stopped_shifted <- make_case("TAAA")
  expect_identical(stopped_shifted$termination, "STOP_CODON")
  expect_length(stopped_shifted$novel_aa, 0L)
  expect_length(Neoantimon2:::source_peptide_windows(stopped_shifted,
    Neoantimon2:::peptide_policy()), 0L)
  shifted <- make_case("AAAA")
  expect_identical(shifted$frame_state, "out_of_frame")
  expect_true(shifted$termination %in% c("TRANSCRIPT_END_NO_STOP",
    "INCOMPLETE_TERMINAL_CODON"))
  expect_gt(length(shifted$novel_aa), 0L)
})

test_that("fusion in-frame novelty remains junction-specific", {
  protein <- "MKPGKPGKPGKPGKPGKPGK"
  in_frame <- list(status = "TRANSLATED", protein = protein,
    junction_aa = list(10L), novel_aa = list(), frame_state = "in_frame")
  policy <- Neoantimon2:::peptide_policy(list(class_i_lengths = list(3L),
    class_ii_lengths = list(4L)))
  windows <- Neoantimon2:::source_peptide_windows(in_frame, policy)
  expect_true(length(windows) > 0L)
  expect_true(all(vapply(windows, function(p)
    p$start_aa <= 10L && p$end_aa >= 11L, TRUE)))
  shifted <- in_frame
  shifted$frame_state <- "out_of_frame"
  shifted$novel_aa <- as.list(11:20)
  tail <- Neoantimon2:::source_peptide_windows(shifted, policy)
  expect_true(any(vapply(tail, function(p)
    p$start_aa >= 12L && !p$junction_crossing, TRUE)))
})

# Real-shaped deterministic genome fixture; these are specified structures,
# not patient observations. The retained region encodes 14 complete residues.
source_inserted_fixture <- function(phase = 0L, strand = "+") {
  path <- testthat::test_path("..", "fixtures", "source", "inserted-splice.fa")
  fai <- Neoantimon2:::reference_fai(paste0(path, ".fai"), path)
  if (strand == "-") {
    genome <- Neoantimon2:::reconstruction_reverse_complement(readLines(path)[[2L]])
    path <- tempfile(fileext = ".fa")
    writeLines(c(">chrInserted", genome), path)
  }
  reference <- list(con = file(path, "rb"), fai = fai)
  exons <- data.frame(start = c(1L, 79L + phase), end = c(36L + phase, 117L))
  if (strand == "-") exons <- data.frame(start = 118L - rev(exons$end),
    end = 118L - rev(exons$start))
  coordinates <- Neoantimon2:::translation_coordinates(exons, strand)
  sequence <- paste0(vapply(coordinates, function(p)
    Neoantimon2:::event_reference_fetch(reference, "chrInserted", p, p), ""),
    collapse = "")
  if (strand == "-") sequence <- chartr("ACGT", "TGCA", sequence)
  start <- 1L
  model <- list(model = list(gene_id = "GINSERT", contig = "chrInserted",
    strand = strand, transcript_id = "TXINSERT", exons = exons),
    sequence = sequence, coordinates = coordinates, coding_start = start,
    coding_start_coordinate = coordinates[[start]], sec_sites = list())
  list(model = model, reference = reference)
}

test_that("in-frame inserted splice sequence emits internal novel peptides", {
  fixture <- source_inserted_fixture()
  on.exit(close(fixture$reference$con), add = TRUE)
  events <- Neoantimon2:::source_splice_rows(list(input_id = "inserted",
    sample_id = "sample-1", version = "fixture-v1", resolved_path =
      testthat::test_path("..", "fixtures", "source", "inserted-splice.tsv")),
    "GRCh38.p14")
  policy <- Neoantimon2:::peptide_policy(list(class_i_lengths = list(8L),
    class_ii_lengths = list(13L)))
  for (event in events) {
    altered <- Neoantimon2:::source_splice_structure(event, fixture$model,
      fixture$reference)
    orf <- Neoantimon2:::source_orf(event, altered, event$source_event_id)
    interval <- altered$novel_nt_intervals[[1L]]
    retained <- event$source_type == "retained_intron"
    expect_identical(interval$start_nt, 37L)
    expect_identical(interval$end_nt, if (retained) 78L else 72L)
    expect_identical(interval$reference_coordinates, as.list(if (retained)
      37:78 else if (event$source_type == "alternative_donor") 37:72 else 43:78))
    expect_identical(interval$contig, "chrInserted")
    expect_identical(interval$strand, "+")
    expect_identical(orf$frame_state, "in_frame")
    expect_identical(orf$termination, "STOP_CODON")
    expect_identical(unlist(orf$inserted_sequence_aa), if (retained) 13:26 else 13:24)
    expect_identical(orf$novel_aa, orf$inserted_sequence_aa)
    windows <- Neoantimon2:::source_peptide_windows(orf, policy)
    internal <- Filter(function(p) p$start_aa == 14L && p$end_aa == 21L, windows)
    expect_length(internal, 1L)
    expect_true(internal[[1L]]$novel_residue)
    expect_false(internal[[1L]]$junction_crossing)
    expect_identical(internal[[1L]]$reference_comparison,
      "inserted_sequence_no_single_reference")
    expect_true(any(vapply(windows, `[[`, TRUE, "junction_crossing")))
    expect_true(any(vapply(windows, function(p)
      p$start_aa < 13L && p$end_aa >= 13L && p$novel_residue, TRUE)))
    expect_false(any(vapply(windows, function(p)
      p$end_aa <= 12L || p$start_aa > max(unlist(orf$novel_aa)), TRUE)))
    # Exact self matching is independent of structural novelty.
    proteome <- tempfile(fileext = ".fa")
    writeLines(c(">unrelated_reference_protein", internal[[1L]]$sequence), proteome)
    self <- Neoantimon2:::peptide_self_search(internal[[1L]]$sequence,
      proteome, list(), sec_encoding = "LITERAL_U")
    expect_true(self[[internal[[1L]]$sequence]]$exact_match)
    project <- source_fixture_project(list(), hla = list("patient-1" = "HLA-A*02:01"))
    origin <- Neoantimon2:::source_origin(project, event, "inserted-fixture",
      altered, orf, internal[[1L]])
    expect_null(origin$wt_peptide_sequence)
    expect_identical(origin$wt_comparison_status, "inserted_sequence_no_single_reference")
  }
})

test_that("inserted nucleotide overlap includes both boundary codons", {
  for (phase in 1:2) {
    fixture <- source_inserted_fixture(phase)
    event <- list(source_type = "retained_intron", structure = list(
      source_type = "retained_intron", gene_id = "GINSERT", contig = "chrInserted",
      strand = "+", transcript_id = "TXINSERT", donor = 36L + phase,
      acceptor = 79L + phase, retained_start = 37L + phase,
      retained_end = 78L + phase))
    altered <- Neoantimon2:::source_splice_structure(event, fixture$model,
      fixture$reference)
    close(fixture$reference$con)
    orf <- Neoantimon2:::source_orf(event, altered, paste0("boundary-", phase))
    expect_identical(orf$frame_state, "in_frame")
    expect_identical(unlist(orf$inserted_sequence_aa), 13:27)
    expect_identical(unlist(orf$novel_aa), 13:27)
    expect_false(28L %in% unlist(orf$novel_aa))
  }
})

test_that("minus-strand inserted intervals retain oriented genomic provenance", {
  fixture <- source_inserted_fixture(strand = "-")
  on.exit(close(fixture$reference$con), add = TRUE)
  event <- list(source_type = "retained_intron", structure = list(
    source_type = "retained_intron", gene_id = "GINSERT", contig = "chrInserted",
    strand = "-", transcript_id = "TXINSERT", donor = 82L, acceptor = 39L,
    retained_start = 40L, retained_end = 81L))
  altered <- Neoantimon2:::source_splice_structure(event, fixture$model,
    fixture$reference)
  interval <- altered$novel_nt_intervals[[1L]]
  expect_identical(interval$reference_coordinates, as.list(81:40))
  expect_identical(interval$start_nt, 37L)
  expect_identical(interval$end_nt, 78L)
  orf <- Neoantimon2:::source_orf(event, altered, "minus-inserted")
  expect_identical(unlist(orf$novel_aa), 13:26)
  expect_identical(substring(orf$protein, 13L, 26L), "AGSVLDEFRYTHQN")
})

test_that("deletion-only donor and acceptor do not invent inserted novelty", {
  fixture <- source_inserted_fixture()
  on.exit(close(fixture$reference$con), add = TRUE)
  for (type in c("alternative_donor", "alternative_acceptor")) {
    event <- list(source_type = type, structure = list(source_type = type,
      gene_id = "GINSERT", contig = "chrInserted", strand = "+",
      transcript_id = "TXINSERT", donor = if (type == "alternative_donor") 33L else 36L,
      acceptor = if (type == "alternative_acceptor") 82L else 79L))
    altered <- Neoantimon2:::source_splice_structure(event, fixture$model,
      fixture$reference)
    orf <- Neoantimon2:::source_orf(event, altered, type)
    expect_length(altered$novel_nt_intervals, 0L)
    expect_identical(orf$frame_state, "in_frame")
    expect_length(orf$novel_aa, 0L)
    windows <- Neoantimon2:::source_peptide_windows(orf, Neoantimon2:::peptide_policy())
    expect_true(length(windows) > 0L)
    expect_true(all(vapply(windows, `[[`, TRUE, "junction_crossing")))
  }
})

test_that("source biological IDs survive additive overlay evolution", {
  kinds <- c("source_event_id", "fusion_event_id", "transcript_structure_id",
    "altered_transcript_id", "orf_hypothesis_id", "protein_origin_id")
  payload <- list(structure = "unchanged", transcripts = list("A", "B"))
  before <- lapply(kinds, Neoantimon2:::source_id, payload = payload)
  expected <- lapply(kinds, function(kind) paste0(kind, ":",
    Neoantimon2:::infra_hash(list(version = "source-biological-id-v1",
      kind = kind, payload = payload))))
  expect_identical(before, expected)
  expect_identical(Neoantimon2:::source_identity_schema_version(), "source-biological-id-v1")
  p <- source_fixture_project(list(source_fixture_input("star-fusion.tsv", "star_fusion")),
    hla = list("patient-1" = "HLA-A*02:01"))
  old <- Neoantimon2:::source_extension_construct(p)
  testthat::local_mocked_bindings(source_extension_version = function() "source-event-2.1.0",
    .package = "Neoantimon2")
  expect_identical(lapply(kinds, Neoantimon2:::source_id, payload = payload), before)
  new <- Neoantimon2:::source_extension_construct(p)
  expect_identical(new$schema_version, "source-event-2.1.0")
  expect_identical(new$identity_schema_version, old$identity_schema_version)
  for (field in c("canonical_events", "hypotheses", "peptides", "origins", "candidates"))
    expect_identical(new[[field]], old[[field]])
})

test_that("provider evidence and field order do not perturb biological IDs", {
  input <- source_fixture_input("arriba.tsv", "arriba")
  original <- source_fixture_project(list(input), hla = list("patient-1" = "HLA-A*02:01"))
  before <- Neoantimon2:::source_extension_construct(original)
  fields <- lapply(readLines(input$path), function(line)
    strsplit(line, "\t", fixed = TRUE)[[1L]])
  header <- sub("^#", "", fields[[1L]])
  fields[[1L]] <- header
  fields[[2L]][match("split_reads1", header)] <- "101"
  fields[[2L]][match("confidence", header)] <- "low"
  changed <- tempfile(fileext = ".tsv")
  writeLines(vapply(fields, function(row)
    paste(c(rev(row), "additive_metadata"), collapse = "\t"), ""), changed)
  input$path <- changed
  after <- Neoantimon2:::source_extension_construct(source_fixture_project(list(input),
    hla = list("patient-1" = "HLA-A*02:01")))
  expect_false(identical(before$observations[[1L]]$observation_id,
    after$observations[[1L]]$observation_id))
  expect_false(identical(before$observations[[1L]]$support, after$observations[[1L]]$support))
  expect_identical(before$canonical_events[[1L]]$source_event_id,
    after$canonical_events[[1L]]$source_event_id)
  for (field in c("transcript_structure_id", "altered_transcript_id"))
    expect_identical(before$hypotheses[[1L]][[field]], after$hypotheses[[1L]][[field]])
  expect_identical(before$hypotheses[[1L]]$orf, after$hypotheses[[1L]]$orf)
  for (field in c("peptide_key", "sample_peptide_id"))
    expect_identical(vapply(before$origins, `[[`, "", field),
      vapply(after$origins, `[[`, "", field))
  for (field in c("candidate_key", "candidate_id"))
    expect_identical(vapply(before$candidates, `[[`, "", field),
      vapply(after$candidates, `[[`, "", field))
  # Candidate IDs from the approved pre-fix fixture output remain exact.
  expect_setequal(vapply(after$candidates, `[[`, "", "candidate_id"), c(
    "candidate_id:69f27013477f4c1d5d5e24b5c42b616966e7377cb2f3041174eda5da723f4718",
    "candidate_id:0398205cb2a1b1571a688d0f1b1965161a296bf313d0795e26c7242eaeaf752a",
    "candidate_id:a3bd6cbc37be2941f6b887a11d585d980e40b9888dba4c10d17eced7f5908818"))
})

test_that("frozen peptide grouping retains independent source origins", {
  contract <- Neoantimon2:::peptide_contract()
  sequence <- "MKPGKPGK"
  key <- Neoantimon2:::biological_id("peptide_key",
    list(ligand_sequence = sequence), contract)
  sample <- Neoantimon2:::biological_id("sample_peptide_id",
    list(context_id = "somatic:patient-1", patient_id = "patient-1",
      sample_id = "sample-1", peptide_key = key), contract)
  template <- list(transcript_id = "TX", transcript_version = NULL,
    reconstructed_sequence_checksum = "checksum", peptide_start_aa = 1L,
    peptide_end_aa = 8L, source_event_keys = list("small-variant-event"),
    phase_namespace = NULL, phase_set = NULL, haplotype_id = "haplotype",
    source_subject_id = "patient-1", source_antigen_key = NULL,
    scenario_id = "default")
  origins <- lapply(c("small_variant", "fusion", "splice"), function(type) {
    x <- template
    x$source_event_keys <- list(type)
    x$source_antigen_key <- type
    Neoantimon2:::biological_id("origin_id",
      list(sample_peptide_id = sample, origin_identity = x), contract)
  })
  expect_length(unique(unlist(origins)), 3L)
  candidate <- function() Neoantimon2:::biological_id("candidate_id",
    list(context_id = "somatic:patient-1", sample_id = "sample-1",
      candidate_key = Neoantimon2:::biological_id("candidate_key",
        list(species = "NCBITaxon:9606", hla_class = "I",
          ligand_sequence = sequence, mhc_key = "HLA-A*02:01"), contract)), contract)
  expect_identical(candidate(), candidate())
})
