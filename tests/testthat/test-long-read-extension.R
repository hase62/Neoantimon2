# Independent release gates exercise observable structure, coordinates and IDs.
lr_cache <- new.env(parent = emptyenv())
lr_all <- function() {
  if (is.null(lr_cache$all)) lr_cache$all <- long_read_fixture_construct()
  lr_cache$all
}
lr_denovo <- function(file = "models.gtf") long_read_fixture_construct(list(
  long_read_fixture_input(file), long_read_fixture_input("settings-denovo.json", "long_read_settings", "settings")))

test_that("generic GTF and GFF share exact oriented structure", {
  root <- test_path("..", "fixtures", "long_read")
  gtf <- Neoantimon2:::long_read_models(file.path(root, "fsm.gtf"), "GRCh38.p14")[[1L]]
  gff <- Neoantimon2:::long_read_models(file.path(root, "models.gff3"), "GRCh38.p14")[[1L]]
  expect_identical(gtf$structure, gff$structure)
  expect_identical(gtf$structure$transcript_start, 1L)
  expect_identical(gtf$structure$transcript_end, 195L)
  expect_identical(gtf$structure$junctions, list(list(donor = 36L, acceptor = 79L),
    list(donor = 114L, acceptor = 157L)))
  shuffled <- tempfile(fileext = ".gtf")
  writeLines(rev(readLines(file.path(root, "fsm.gtf"))), shuffled)
  expect_identical(Neoantimon2:::long_read_models(shuffled, "GRCh38.p14")[[1L]]$structure, gtf$structure)
})

test_that("malformed models fail rather than invent exon sequence", {
  bad <- tempfile(fileext = ".gtf")
  writeLines('chrLR\tx\ttranscript\t1\t30\t.\t+\t.\ttranscript_id "bad";', bad)
  expect_error(Neoantimon2:::long_read_models(bad, "GRCh38.p14"), "Missing required")
  writeLines('chrLR\tx\texon\t1\t30\t.\t+\t.\tgene_id "G";', bad)
  expect_error(Neoantimon2:::long_read_models(bad, "GRCh38.p14"), "Missing transcript")
  expect_error(Neoantimon2:::long_read_structure("GRCh38", "chrLR", "+",
    data.frame(start = c(1L, 10L), end = c(12L, 20L))), "Overlapping")
  writeLines('chrLR\tx\texon\t0\t30\t.\t+\t.\ttranscript_id "bad";', bad)
  expect_error(Neoantimon2:::long_read_models(bad, "GRCh38"), "Invalid")
})

test_that("SQANTI3 FLAIR and IsoQuant observations share one structure", {
  inputs <- unlist(lapply(c("sqanti3", "flair", "isoquant"), function(p) list(
    long_read_fixture_input(paste0(p, ".gtf"), paste0(p, "_gtf"), p),
    long_read_fixture_input(paste0(p, ".tsv"), paste0(p, "_metadata"), paste0(p, "-meta")))), recursive = FALSE)
  x <- long_read_fixture_construct(inputs)
  y <- long_read_fixture_construct(rev(inputs))
  expect_length(x$structures, 1L)
  expect_length(x$observations, 3L)
  expect_length(unique(vapply(x$observations, `[[`, "", "isoform_structure_id")), 1L)
  expect_length(unique(vapply(x$observations, `[[`, "", "long_read_observation_id")), 3L)
  expect_identical(x, y)
  sq <- Filter(function(o) o$provider == "sqanti3", x$observations)[[1L]]
  expect_identical(sq$metadata[[1L]]$raw_fields$unknown_future_column, "retained")
  expect_identical(sq$metadata[[1L]]$raw_fields$RTS_stage, "FALSE")
  expect_identical(sq$provider_version, "synthetic-v1")
  flair <- Filter(function(o) o$provider == "flair", x$observations)[[1L]]
  expect_identical(flair$metadata[[1L]]$raw_fields$assigned_reads, "0")
  expect_identical(flair$metadata[[1L]]$raw_fields$TPM, "NA")
  iso <- Filter(function(o) o$provider == "isoquant", x$observations)[[1L]]
  expect_identical(iso$metadata[[1L]]$raw_fields$count, "14")
  expect_true(length(x$candidates) > 0L)
  expect_true(any(vapply(x$hypotheses, function(h) h$orf$evidence_type == "provider_supplied_orf", TRUE)))
})

test_that("provider metadata uses headers and preserves group-specific support", {
  x <- long_read_fixture_construct(list(long_read_fixture_input("isoquant.gtf", "isoquant_gtf"),
    long_read_fixture_input("isoquant-grouped.tsv", "isoquant_metadata", "metadata")))
  expect_length(x$observations[[1L]]$metadata, 2L)
  expect_identical(vapply(x$observations[[1L]]$metadata, function(m) m$raw_fields$count, ""), c("0", "3"))
  bad <- tempfile(fileext = ".tsv")
  writeLines(c("isoform\tchrom\tstrand\tstructural_category", "x\tchrLR\t+\tNNC"), bad)
  expect_error(Neoantimon2:::long_read_metadata(list(input_role = "sqanti3_metadata", resolved_path = bad)), "Unsupported")
})

test_that("biological IDs ignore provider evidence and additive schema fields", {
  s <- Filter(function(o) o$provider_transcript_id == "fsm", lr_all()$observations)[[1L]]$model$structure
  original <- Neoantimon2:::long_read_structure_id(s)
  for (field in c("provider", "transcript_id", "read_count", "TPM", "classification", "schema_version")) {
    s[[field]] <- "changed"
    expect_identical(Neoantimon2:::long_read_structure_id(s), original)
  }
  for (field in c("strand", "reference_build")) {
    changed <- s; changed[[field]] <- if (field == "strand") if (s$strand == "+") "-" else "+" else "GRCh37"
    expect_false(identical(Neoantimon2:::long_read_structure_id(changed), original))
  }
  changed <- s; changed$exons[[1L]]$end <- changed$exons[[1L]]$end + 1L
  expect_false(identical(Neoantimon2:::long_read_structure_id(changed), original))
  changed <- s; changed$exons <- rev(changed$exons)
  expect_false(identical(Neoantimon2:::long_read_structure_id(changed), original))
})

test_that("reference comparison preserves all hypotheses and provider labels stay evidence", {
  x <- lr_all(); fsm <- long_read_test_event(x, "fsm")
  expect_identical(fsm$structural_comparison$status, "REFERENCE_EQUIVALENT")
  expect_gte(length(fsm$structural_comparison$comparisons), 4L)
  expect_identical(sum(vapply(fsm$structural_comparison$comparisons, `[[`, TRUE, "exact_structure")), 2L)
  expect_length(long_read_test_orfs(x, "fsm"), 4L)
  for (orf in long_read_test_orfs(x, "fsm")) expect_false(orf$peptide_eligible)
  expect_identical(long_read_test_event(x, "no-reference")$structural_comparison$status, "UNRESOLVED_COMPARISON")
  expect_identical(long_read_test_event(x, "antisense")$structural_comparison$status, "UNRESOLVED_COMPARISON")
  expect_length(long_read_test_orfs(x, "intergenic"), 0L)
})

test_that("reference genome reconstruction maps plus minus single and multi exon", {
  x <- lr_all()
  for (tx in c("fsm", "minus", "single", "retained", "novel-exon", "terminal")) {
    e <- long_read_test_event(x, tx); r <- e$altered_transcript
    expect_identical(nchar(r$sequence), length(r$coordinates))
    expect_identical(r$sequence_authority, "REFERENCE_GENOME")
    for (b in r$blocks) expect_identical(b$end_nt - b$start_nt + 1L, b$end - b$start + 1L)
  }
  expect_identical(long_read_test_event(x, "fsm")$altered_transcript$sequence,
    long_read_test_event(x, "minus")$altered_transcript$sequence)
  expect_identical(long_read_test_event(x, "minus")$altered_transcript$coordinates[1:3], 300:298)
})

test_that("structural classes derive independently of labels", {
  x <- lr_all()
  cls <- function(tx) unlist(lapply(long_read_test_event(x, tx)$structural_comparison$comparisons, `[[`, "classes"))
  expect_true("known_junction_combination" %in% cls("nic"))
  expect_true("novel_donor" %in% cls("donor"))
  expect_true("novel_acceptor" %in% cls("acceptor"))
  expect_true("retained_intronic_block" %in% cls("retained"))
  expect_true("exon_insertion" %in% cls("novel-exon"))
  expect_true("exon_skipping" %in% cls("skip"))
  expect_true("alternative_terminal_structure" %in% cls("terminal"))
  expect_true("reference-compatible_structure" %in% cls("fsm"))
})

test_that("raw FASTA one-base and length differences never generate mutation truth", {
  baseline <- long_read_fixture_construct(list(long_read_fixture_input("fsm.gtf")))
  for (name in c("exact", "mismatch", "length")) {
    x <- long_read_fixture_construct(list(long_read_fixture_input("fsm.gtf"),
      long_read_fixture_input(paste0(name, ".fa"), "long_read_fasta", "fasta")))
    expected <- c(exact = "exact_match", mismatch = "sequence_difference", length = "length_difference")[[name]]
    expect_identical(x$observations[[1L]]$fasta_comparison$status, expected)
    expect_false(x$observations[[1L]]$fasta_comparison$mutation_inference)
    expect_identical(x$events[[1L]]$altered_transcript$sequence, baseline$events[[1L]]$altered_transcript$sequence)
    expect_identical(lapply(x$hypotheses, `[[`, "orf"), lapply(baseline$hypotheses, `[[`, "orf"))
    expect_false(identical(x$observations[[1L]]$long_read_observation_id,
      baseline$observations[[1L]]$long_read_observation_id))
    expect_length(x$peptides, 0L)
    expect_length(x$candidates, 0L)
  }
  expect_identical(baseline$observations[[1L]]$fasta_comparison$status, "not_supplied")
})

test_that("annotated start retention and loss remain explicit", {
  x <- lr_all()
  for (orf in long_read_test_orfs(x, "novel-exon")) {
    expect_identical(orf$start_nt, 1L)
    expect_identical(orf$start_codon, "ATG")
    expect_identical(orf$policy, "long-read-reference-orf-v1")
    expect_identical(orf$termination, "STOP_CODON")
  }
  lost <- long_read_test_event(x, "lost-start")
  expect_length(long_read_test_orfs(x, "lost-start"), 0L)
  expect_true(all(vapply(lost$orf_unavailable, `[[`, "", "reason") == "ANNOTATED_START_NOT_RETAINED"))
})

test_that("provider CDS validates coordinates start frame bounds and stop", {
  rec <- list(sequence = "ATGGCTGCTTAA", strand = "+")
  valid <- list(start = 1L, end = 9L, strand = "+", frame = 0L)
  expect_identical(Neoantimon2:::long_read_provider_orf_validate(valid, rec)$status, "valid")
  for (field in c("start", "end", "strand", "frame")) {
    p <- valid; p[[field]] <- switch(field, start = 2L, end = 8L, strand = "-", frame = 1L)
    expect_identical(Neoantimon2:::long_read_provider_orf_validate(p, rec)$status, "inconsistent")
  }
  p <- valid; p$start <- 0L
  expect_identical(Neoantimon2:::long_read_provider_orf_validate(p, rec)$reason, "START_OUT_OF_BOUNDS")
  p <- valid; p$end <- 99L
  expect_identical(Neoantimon2:::long_read_provider_orf_validate(p, rec)$reason, "CDS_END_INCONSISTENT")
})

test_that("de novo AUG enumeration is opt in and FSM safety gate is absolute", {
  expect_false(any(vapply(lr_all()$hypotheses, function(h) h$orf$evidence_type == "de_novo_aug_orf", TRUE)))
  x <- lr_denovo("fsm.gtf")
  orfs <- long_read_test_orfs(x, "fsm", "de_novo_aug_orf")
  expect_gt(length(orfs), 1L)
  starts <- gregexpr("ATG", x$events[[1L]]$altered_transcript$sequence, fixed = TRUE)[[1L]]
  expect_setequal(vapply(orfs, `[[`, 0L, "start_nt"), as.integer(starts))
  expect_length(unique(vapply(orfs, `[[`, "", "orf_hypothesis_id")), length(starts))
  expect_true(all(vapply(orfs, function(o) o$start_codon == "ATG" && !o$peptide_eligible, TRUE)))
  expect_length(x$peptides, 0L)
  expect_length(x$candidates, 0L)
})

test_that("AUG only scanner preserves stop boundary frame and Kozak evidence", {
  sequence <- "CTGGTGTTGACGATGAAATAAATGGCTATGGCTGCT"
  rec <- list(sequence = sequence, coordinates = seq_len(nchar(sequence)), strand = "+")
  cmp <- list(reference_equivalent = FALSE, status = "UNRESOLVED_COMPARISON",
    comparisons = list(), novel_junction_nt = integer())
  starts <- gregexpr("ATG", sequence, fixed = TRUE)[[1L]]
  orfs <- lapply(starts, function(p) Neoantimon2:::long_read_translate_orf("iso", rec, cmp,
    p, "de_novo_aug_orf", "noncanonical-orf-v1"))
  expect_identical(vapply(orfs, `[[`, 0L, "start_nt"), c(13L, 22L, 28L))
  expect_identical(orfs[[1L]]$protein, "MK")
  expect_identical(orfs[[1L]]$stop_nt, 19L)
  expect_identical(orfs[[1L]]$stop_status, "first_in_frame_stop")
  expect_identical(orfs[[2L]]$boundary_status, "transcript_boundary")
  expect_identical(orfs[[1L]]$start_context$minus_3, "A")
  expect_identical(orfs[[1L]]$start_context$plus_4, "A")
  expect_false(orfs[[1L]]$start_context$hard_filter)
  changed <- Neoantimon2:::long_read_translate_orf("iso", rec, cmp, 13L, "de_novo_aug_orf", "noncanonical-orf-v2")
  expect_false(identical(changed$orf_hypothesis_id, orfs[[1L]]$orf_hypothesis_id))
})

test_that("internal novel exon and retained intron peptides need no junction", {
  x <- lr_all()
  for (tx in c("novel-exon", "retained")) {
    e <- long_read_test_event(x, tx)
    expect_true(all(vapply(e$altered_transcript$novel_nt_intervals, function(nt)
      !is.null(nt$start_nt) && !is.null(nt$end_nt) && length(nt$reference_coordinates) > 0L, TRUE)))
    expect_gt(length(e$altered_transcript$novel_nt_intervals), 0L)
    peptides <- Filter(function(p) p$long_read_source_event_id == e$long_read_source_event_id, x$peptides)
    internal <- Filter(function(p) p$peptide$novel_residue && !p$peptide$junction_crossing, peptides)
    expect_gt(length(internal), 0L)
    expect_true(any(vapply(peptides, function(p) p$peptide$junction_crossing, TRUE)))
    orf <- long_read_test_orfs(x, tx)[[1L]]
    for (path in orf$novelty_paths) for (nt in path$novel_nt_intervals) {
      codons <- orf$start_nt + (seq_len(nchar(orf$protein)) - 1L) * 3L
      expected <- which(codons <= nt$end_nt & codons + 2L >= nt$start_nt)
      expect_true(all(expected %in% unlist(orf$novel_aa)))
    }
  }
})

test_that("in frame deletion preserves downstream residues while frameshifts retain tails", {
  x <- lr_all()
  skip <- long_read_test_orfs(x, "skip")
  for (o in skip) {
    expect_identical(o$frame_state, "in_frame")
    expect_length(o$inserted_sequence_aa, 0L)
    expect_length(o$novel_frame_aa, 0L)
    expect_true(length(o$junction_aa) > 0L)
  }
  shift <- long_read_test_orfs(x, "frameshift")
  expect_true(any(vapply(shift, function(o) length(o$novel_frame_aa) > 0L, TRUE)))
  expect_true(any(vapply(shift, function(o) o$frame_state == "out_of_frame", TRUE)))
})

test_that("de novo ORFs need structural anchors, including for intergenic isoforms", {
  x <- lr_denovo()
  novel <- long_read_test_orfs(x, "novel-exon", "de_novo_aug_orf")
  expect_true(any(vapply(novel, `[[`, TRUE, "peptide_eligible")))
  intergenic <- long_read_test_orfs(x, "intergenic", "de_novo_aug_orf")
  expect_gt(length(intergenic), 0L)
  expect_true(all(!vapply(intergenic, `[[`, TRUE, "peptide_eligible")))
  relation <- unique(unlist(lapply(novel, function(o) lapply(o$reference_orf_relationships, `[[`, "relationship"))))
  expect_true("downstream_ORF" %in% relation)
})

test_that("matched normal preserves shared structure and independent observations", {
  tumor <- long_read_fixture_input("novel-exon.gtf", id = "tumor")
  normal <- long_read_fixture_input("novel-exon.gtf", id = "normal", sample = "normal-1")
  settings <- tempfile(fileext = ".json")
  jsonlite::write_json(list(version = "long-read-run-settings-v1", bindings = list(
    list(model_input_id = "normal", observation_role = "normal", region = "normal-region"),
    list(model_input_id = "tumor", observation_role = "tumor", region = "tumor-region"))), settings, auto_unbox = TRUE)
  opts <- list(input_id = "settings", input_role = "long_read_settings", sample_id = "sample-1",
    path = settings, version = "v1")
  p <- long_read_fixture_project(list(tumor, normal, opts), samples = data.frame(
    sample_id = c("sample-1", "normal-1"), patient_id = "patient-1"))
  x <- Neoantimon2:::long_read_extension_construct(p)
  expect_length(x$structures, 1L); expect_length(x$observations, 2L); expect_length(x$events, 2L)
  expect_length(unique(vapply(x$events, `[[`, "", "long_read_source_event_id")), 2L)
  expect_true(all(vapply(x$events, `[[`, "", "matched_normal_status") == "normal_observed"))
  expect_gt(length(x$candidates), 0L)
  expect_true(all(vapply(x$candidates, `[[`, "", "sample_id") == "sample-1"))
  expect_true(all(vapply(x$events, function(e) length(e$normal_observation_ids) == 1L, TRUE)))
})

test_that("chimeric models remain evidence and never use invented single gene sequence", {
  path <- tempfile(fileext = ".gtf")
  writeLines(c('chr1\tx\texon\t1\t8\t.\t+\t.\tgene_id "G1"; transcript_id "fusion";',
    'chr2\tx\texon\t21\t28\t.\t-\t.\tgene_id "G2"; transcript_id "fusion";'), path)
  model <- Neoantimon2:::long_read_models(path, "GRCh38.p14")[[1L]]
  expect_identical(model$status, "CHIMERIC_STRUCTURE_DEFERRED")
  expect_null(model$structure)
})

test_that("exact self evidence and frozen candidate IDs remain source independent", {
  x <- lr_all()
  expect_true(all(vapply(x$candidates, `[[`, "", "final_grade") == "NOT_GRADED"))
  expect_true(any(vapply(x$candidates, function(c) c$hla_class == "II", TRUE)))
  expect_true(all(vapply(x$candidates, function(c) !is.null(c$self_status), TRUE)))
  expect_true(any(vapply(x$candidates, `[[`, TRUE, "self_exact_match")))
  expect_true(all(vapply(x$candidates, function(c) grepl("long_read_extension.rds",
    c$self_match_details_ref, fixed = TRUE), TRUE)))
  p <- long_read_fixture_project()
  c <- x$candidates[[1L]]
  key <- Neoantimon2:::biological_id("candidate_key", list(species = "NCBITaxon:9606",
    hla_class = c$hla_class, ligand_sequence = c$ligand_sequence, mhc_key = c$mhc_key), Neoantimon2:::peptide_contract())
  expect_identical(c$candidate_key, key)
  expect_identical(Neoantimon2:::biological_id("candidate_id", list(context_id = c$analysis_context,
    sample_id = c$sample_id, candidate_key = key), Neoantimon2:::peptide_contract()), c$candidate_id)
})

test_that("long-read artifact and ORF option changes invalidate fingerprints", {
  p <- long_read_fixture_project(list(long_read_fixture_input("fsm.gtf")))
  fp <- Neoantimon2:::long_read_fingerprints(p)
  expect_length(Neoantimon2:::long_read_canonical_input_ids(p), 1L)
  opts <- long_read_fixture_input("settings-denovo.json", "long_read_settings", "settings")
  q <- long_read_fixture_project(list(long_read_fixture_input("fsm.gtf"), opts))
  expect_false(identical(fp, Neoantimon2:::long_read_fingerprints(q)))
  expect_identical(Neoantimon2:::long_read_version(), "long-read-source-1.0.0")
  expect_identical(Neoantimon2:::long_read_identity_version(), "long-read-isoform-id-v1")
})

test_that("provider CDS failures remain evidence and GTF CDS shares ORF validation", {
  rec <- list(sequence = "ATGGCTGCTTAA", coordinates = 1:12, strand = "+")
  s <- Neoantimon2:::long_read_structure("GRCh38", "chrTest", "+", data.frame(start = 1L, end = 12L))
  cmp <- list(reference_equivalent = FALSE, status = "UNRESOLVED_COMPARISON",
    comparisons = list(), novel_junction_nt = integer())
  observation <- list(long_read_observation_id = "obs", model = list(structure = s,
    cds_features = list(list(start = 1L, end = 9L, phase = "0"))),
    metadata = list(list(input_id = "meta", raw_fields = list(CDS_start = "2", CDS_end = "9", strand = "+"))),
    de_novo_aug = FALSE)
  orfs <- Neoantimon2:::long_read_orfs("iso", rec, cmp, list(observation))
  expect_length(orfs$hypotheses, 1L)
  expect_identical(orfs$hypotheses[[1L]]$evidence_type, "provider_supplied_orf")
  expect_identical(orfs$hypotheses[[1L]]$protein, "MAA")
  expect_false(orfs$hypotheses[[1L]]$peptide_eligible)
  expect_identical(orfs$unavailable[[1L]]$provider_orf_status, "inconsistent")
  expect_identical(orfs$unavailable[[1L]]$reason, "START_CODON_INCONSISTENT")
})

test_that("uORF overlapping and downstream relationships remain unranked hypotheses", {
  bases <- rep("C", 60L)
  for (p in c(1L, 13L, 20L, 31L)) bases[p + 0:2] <- c("A", "T", "G")
  bases[7:9] <- c("T", "A", "A"); bases[58:60] <- c("T", "A", "A")
  seq <- paste0(bases, collapse = "")
  rec <- list(sequence = seq, coordinates = seq_len(nchar(seq)), strand = "+")
  r <- list(transcript_id = "ref", reconstructed = rec, coding_start = 13L,
    coding_start_coordinates = 13:15)
  cmp <- list(reference_equivalent = TRUE, status = "REFERENCE_EQUIVALENT",
    comparisons = list(list(reference = r, exact_structure = TRUE)), novel_junction_nt = integer())
  starts <- gregexpr("ATG", seq, fixed = TRUE)[[1L]]
  orfs <- lapply(starts, function(p) Neoantimon2:::long_read_translate_orf("iso", rec, cmp,
    p, "de_novo_aug_orf", "noncanonical-orf-v1"))
  classes <- vapply(orfs, function(o) o$reference_orf_relationships[[1L]]$relationship, "")
  expect_true("uORF" %in% classes)
  expect_true("overlapping_alternative_frame_ORF" %in% classes)
  expect_true("downstream_ORF" %in% classes)
  expect_true(all(!vapply(orfs, `[[`, TRUE, "peptide_eligible")))
  expect_length(unique(vapply(orfs, `[[`, "", "orf_hypothesis_id")), length(starts))
  expect_gt(length(unique(vapply(orfs, `[[`, 0L, "frame"))), 1L)
})

test_that("provider metadata reordering count and label changes preserve structure", {
  input <- long_read_fixture_input("sqanti3.gtf", "sqanti3_gtf")
  meta <- long_read_fixture_input("sqanti3.tsv", "sqanti3_metadata", "meta")
  p <- long_read_fixture_project(list(input, meta), hla = NULL)
  a <- Neoantimon2:::long_read_observations(p, Neoantimon2:::reference_bundle(p, retain_indexes = TRUE))[[1L]]
  lines <- readLines(meta$path)
  fields <- lapply(lines, function(l) strsplit(l, "\t", fixed = TRUE)[[1L]])
  fields[[2L]][match("FL", fields[[1L]])] <- "999"
  fields[[2L]][match("structural_category", fields[[1L]])] <- "FSM"
  fields[[2L]][match("iso_exp", fields[[1L]])] <- "88"
  path <- tempfile(fileext = ".tsv")
  writeLines(vapply(fields, function(row) paste(rev(row), collapse = "\t"), ""), path)
  meta$path <- path
  q <- long_read_fixture_project(list(input, meta), hla = NULL)
  b <- Neoantimon2:::long_read_observations(q, Neoantimon2:::reference_bundle(q, retain_indexes = TRUE))[[1L]]
  expect_identical(a$isoform_structure_id, b$isoform_structure_id)
  expect_identical(a$long_read_source_event_id, b$long_read_source_event_id)
  expect_false(identical(a$long_read_observation_id, b$long_read_observation_id))
  expect_identical(b$metadata[[1L]]$raw_fields$FL, "999")
})

test_that("AA20 plus authoritative U support does not infer novel Sec", {
  p <- long_read_fixture_project()
  bundle <- Neoantimon2:::reference_bundle(p, retain_indexes = TRUE)
  reference <- Neoantimon2:::event_reference_open(bundle)
  on.exit(Neoantimon2:::event_reference_close(reference), add = TRUE)
  s <- Neoantimon2:::long_read_structure("GRCh38.p14", "chrS", "+", data.frame(start = 1L, end = 21L))
  rec <- Neoantimon2:::long_read_reconstruct(s, reference)
  refs <- Neoantimon2:::long_read_reference_models(bundle, list(s), reference)
  cmp <- Neoantimon2:::long_read_compare(s, rec, refs)
  o <- Neoantimon2:::long_read_translate_orf("iso", rec, cmp, 1L, "reference_anchored_orf",
    "long-read-reference-orf-v1", "TD.1")
  expect_identical(o$protein, "MUPKPG")
  expect_length(o$translation_exception_ids, 1L)
  d <- Neoantimon2:::long_read_translate_orf("iso", rec, cmp, 1L, "de_novo_aug_orf", "noncanonical-orf-v1")
  expect_identical(d$protein, "M")
  expect_length(d$translation_exception_ids, 0L)
})

test_that("long-read and frozen splice retain shared candidate and distinct origins", {
  model <- tempfile(fileext = ".gtf")
  writeLines(c('chr1\tx\texon\t1\t7\t.\t+\t.\ttranscript_id "shared";',
    'chr1\tx\texon\t21\t33\t.\t+\t.\ttranscript_id "shared";'), model)
  lr <- long_read_fixture_input(); lr$path <- model
  src <- list(input_id = "splice", input_role = "splice_structure", sample_id = "sample-1",
    path = normalizePath(test_path("..", "fixtures", "source", "splice-structure.tsv")), version = "v1")
  p <- long_read_fixture_project(list(lr, src), hla = list("patient-1" = "HLA-A*02:01"))
  values <- unclass(p$config); values$peptides <- list(class_i_lengths = list(6L), class_ii_lengths = list(13L))
  cfg <- do.call(neo_config, c(values, list(env = character())))
  p <- neo_project(samples = p$samples, inputs = list(lr, src), hla = list("patient-1" = "HLA-A*02:01"), config = cfg,
    path = file.path(cfg$paths$run_root, "shared-origins"))
  p <- neo_run(p)
  long <- Neoantimon2:::long_read_read(p); source <- Neoantimon2:::source_read(p)
  shared <- intersect(vapply(long$candidates, `[[`, "", "candidate_id"),
    vapply(source$candidates, `[[`, "", "candidate_id"))
  expect_gt(length(shared), 0L)
  explained <- neo_explain(p, candidate_id = shared[[1L]])$candidates[[1L]]
  expect_true(length(explained$source_extension$provider_observations) > 0L)
  expect_true(length(explained$long_read_extension$provider_observations) > 0L)
  expect_length(unique(vapply(explained$origins, `[[`, "", "origin_id")), length(explained$origins))
  rows <- neo_results(p, candidate_id = shared[[1L]])
  expect_identical(rows$candidate_id, shared[[1L]])
  expect_gte(length(rows$origin_ids[[1L]]), 2L)
  expect_true(grepl("long_read_isoform", rows$source_type, fixed = TRUE))
  expect_gt(nrow(neo_origins(p)), 0L)
})

test_that("fractional CDS inputs fail and nested metadata stays outside identity", {
  expect_identical(Neoantimon2:::long_read_integer_evidence("1.5"), NA_integer_)
  expect_identical(Neoantimon2:::long_read_integer_evidence("NA"), NA_integer_)
  expect_identical(Neoantimon2:::long_read_integer_evidence("0"), 0L)
  s <- Filter(function(o) o$provider_transcript_id == "fsm", lr_all()$observations)[[1L]]$model$structure
  id <- Neoantimon2:::long_read_structure_id(s)
  s$exons[[1L]]$additive_qc <- "future"
  s$junctions[[1L]]$support <- "future"
  expect_identical(Neoantimon2:::long_read_structure_id(s), id)
})
