test_that("synonymous and unavailable proteins generate explicit zero results", {
  syn <- peptide_test_contexts(list(reconstruction_test_row("s", "chr1", 6L, "A", "G")), "TA.1")[[1L]]
  x <- Neoantimon2:::peptide_enumerate(syn)
  expect_identical(x$status, "PROTEIN_IDENTICAL"); expect_length(x$peptides, 0L)
  lost <- peptide_test_contexts(list(reconstruction_test_row("s", "chrA", 1L, "A", "G")), "TV.1")[[1L]]
  y <- Neoantimon2:::peptide_enumerate(lost)
  expect_identical(y$status, "UNAVAILABLE"); expect_length(y$peptides, 0L)
  expect_identical(y$provenance$reconstruction_id, lost$reconstruction_id)
})

test_that("missense windows and WT sequences have independent gold coordinates", {
  ct <- peptide_test_contexts(list(reconstruction_test_row("s", "chr1", 4L, "G", "A")), "TA.1")[[1L]]
  p <- Neoantimon2:::peptide_enumerate(ct, peptide_test_policy(3L, 3L))$peptides
  expect_identical(vapply(p, `[[`, "", "peptide_sequence"), c("MKP", "KPK"))
  expect_identical(vapply(p, function(x) x$normal_correspondence$sequence, ""), c("MEP", "EPK"))
  expect_identical(vapply(p, `[[`, 0L, "peptide_start_aa"), 1:2)
  expect_true(p[[2L]]$altered_at_n_terminus)
  expect_true(p[[1L]]$altered_internal)
  for (p0 in c(2L, 12L, 28L)) {
    ct <- peptide_test_long(p0)
    x <- Neoantimon2:::peptide_enumerate(ct, peptide_test_policy(9L, 15L))$peptides
    for (len in c(9L, 15L)) {
      got <- Filter(function(z) z$peptide_length == len, x)
      expected <- seq.int(max(1L, p0-len+1L), min(p0, nchar(ct$translation$mutant_protein)-len+1L))
      expect_identical(vapply(got, `[[`, 0L, "peptide_start_aa"), expected)
      expect_true(all(vapply(got, function(z) z$overlaps_altered_sequence, TRUE)))
    }
  }
  x <- Neoantimon2:::peptide_enumerate(peptide_test_long(), peptide_test_policy(9L, 15L))$peptides
  expect_true(any(vapply(x, `[[`, TRUE, "altered_at_c_terminus")))
  expect_true(any(vapply(x, `[[`, TRUE, "altered_at_n_terminus")))
  expect_true(any(vapply(x, function(z) "II" %in% unlist(z$classes), TRUE)))
})

test_that("MNV and cis SNVs give equal peptides while trans stays separate", {
  a <- reconstruction_test_row("a", "chr1", 4L, "G", "A", side = "LEFT")
  b <- reconstruction_test_row("b", "chr1", 5L, "A", "C", side = "LEFT")
  m <- reconstruction_test_row("m", "chr1", 4L, "GA", "AC")
  cis <- peptide_test_contexts(list(a,b), "TA.1")
  mnv <- peptide_test_contexts(list(m), "TA.1")
  expect_identical(peptide_test_sequences(cis[[1L]]), c("MTP", "TPK"))
  expect_identical(peptide_test_sequences(cis[[1L]]), peptide_test_sequences(mnv[[1L]]))
  b$event_details_json$phase_assertions[[1L]]$local_side <- "RIGHT"
  trans <- peptide_test_contexts(list(a,b), "TA.1")
  expect_length(trans, 2L)
  expect_false(any(vapply(trans, function(x) "MTP" %in% peptide_test_sequences(x), TRUE)))
  flip <- function(x) { z <- x$event_details_json$phase_assertions[[1L]]$local_side
    x$event_details_json$phase_assertions[[1L]]$local_side <- if (z == "LEFT") "RIGHT" else "LEFT"; x }
  reversed <- peptide_test_contexts(list(flip(b),flip(a)), "TA.1")
  key <- function(xs) sort(vapply(xs, function(x) paste(peptide_test_sequences(x), collapse="/"), ""))
  expect_identical(key(trans), key(reversed))
})

test_that("in-frame insertion deletion and replacement retain junction correspondence", {
  cases <- list(list("A", "ACCC", 6L, c("MEP","EPP","PPK")),
    list("GGAA", "G", 3L, "MPK"), list("GAA", "GTCGAC", 4L, c("MVD","VDP","DPK")))
  for (z in cases) {
    ct <- peptide_test_contexts(list(reconstruction_test_row("s", "chr1", z[[3L]], z[[1L]], z[[2L]])), "TA.1")[[1L]]
    expect_identical(peptide_test_sequences(ct), z[[4L]])
  }
  deletion <- peptide_test_contexts(list(reconstruction_test_row("d", "chr1", 3L, "GGAA", "G")), "TA.1")[[1L]]
  p <- Neoantimon2:::peptide_enumerate(deletion, peptide_test_policy(3L,3L))$peptides[[1L]]
  expect_true(p$overlaps_altered_junction); expect_false(p$overlaps_altered_sequence)
  expect_identical(p$normal_correspondence$status, "JUNCTION_INTERVAL")
  expect_identical(p$normal_correspondence$reference_span, "MEPK")
  expect_null(p$normal_correspondence$sequence)
})

test_that("frameshift termination axes bound peptides without being conflated", {
  cases <- list(list("chr1", "TA.1", "AT", "premature", "STOP_CODON"),
    list("chrW", "TW.1", "AC", "lost_with_recovered_stop", "STOP_CODON"),
    list("chrF", "TF.1", "AC", "lost_open_end", "INCOMPLETE_TERMINAL_CODON"))
  for (z in cases) {
    ct <- peptide_test_contexts(list(reconstruction_test_row("s", z[[1L]], 6L, "A", z[[3L]])), z[[2L]])[[1L]]
    x <- Neoantimon2:::peptide_enumerate(ct, peptide_test_policy(2L,3L))
    expect_identical(x$provenance$frame_effect, "shifted_unrestored")
    expect_identical(x$provenance$protein_effect, "altered_tail")
    expect_identical(x$provenance$termination_effect, z[[4L]])
    expect_identical(x$provenance$mutant_boundary, z[[5L]])
    expect_gt(length(x$peptides), 0L)
    for (p in x$peptides) {
      expect_lte(p$peptide_end_aa, nchar(ct$translation$mutant_protein))
      expect_identical(p$peptide_sequence, substring(ct$translation$mutant_protein, p$peptide_start_aa,p$peptide_end_aa))
      expect_null(p$normal_correspondence$sequence)
    }
  }
})

test_that("restored frames and stop loss use observed sequence and pure truncation yields zero", {
  ct <- peptide_test_contexts(list(reconstruction_test_row("i", "chrF", 6L, "A", "AC", side="LEFT"),
    reconstruction_test_row("d", "chrF", 26L, "CC", "C", side="LEFT")), "TF.1")[[1L]]
  expect_identical(peptide_test_sequences(ct), c("KPR","PRE","RET","ETG","TGK"))
  expect_identical(ct$translation$frame_effect, "shifted_restored")
  plus <- peptide_test_contexts(list(reconstruction_test_row("s", "chrF", 38L, "T", "C")), "TF.1")[[1L]]
  minus <- peptide_test_contexts(list(reconstruction_test_row("s", "chrR", 15L, "A", "G")), "TR.1")[[1L]]
  expect_identical(peptide_test_sequences(plus), c("KPQ","PQP","QPG"))
  expect_identical(peptide_test_sequences(plus), peptide_test_sequences(minus))
  expect_true(all(vapply(Neoantimon2:::peptide_enumerate(plus, peptide_test_policy(3L,3L))$peptides,
    function(p) is.null(p$normal_correspondence$sequence), TRUE)))
  trunc <- peptide_test_contexts(list(reconstruction_test_row("s", "chr1", 4L, "G", "T")), "TA.1")[[1L]]
  x <- Neoantimon2:::peptide_enumerate(trunc, peptide_test_policy(1L,1L))
  expect_identical(x$status, "NO_ALTERED_LIGAND"); expect_length(x$peptides, 0L)
})

test_that("multiple intervals union windows and preserve germline and BOTH contexts", {
  a <- reconstruction_test_row("a", "chrF", 4L, "A", "G", side="LEFT")
  b <- reconstruction_test_row("b", "chrF", 29L, "G", "A", side="LEFT")
  ct <- peptide_test_contexts(list(a,b), "TF.1")[[1L]]
  expect_identical(peptide_test_sequences(ct), c("MEP","EPG","KPR","PRK","RKP"))
  x <- Neoantimon2:::peptide_enumerate(ct, peptide_test_policy(8L,8L))$peptides
  expect_length(x, 2L); expect_identical(x[[1L]]$difference_interval_indices, list(1L,2L))
  germ <- reconstruction_test_row("g", "chr1", 7L, "C", "T", "germline", both=TRUE)
  som <- reconstruction_test_row("s", "chr1", 4L, "G", "A", both=TRUE)
  ct <- peptide_test_contexts(list(germ,som), "TA.1")[[1L]]
  expect_identical(ct$germline_event_ids, "g"); expect_identical(ct$somatic_event_ids, "s")
  p <- Neoantimon2:::peptide_enumerate(ct, peptide_test_policy(3L,3L))$peptides
  expect_identical(p[[1L]]$normal_correspondence$sequence, "MES")
  expect_identical(p[[1L]]$reference_correspondence$sequence, "MEP")
  expect_identical(p[[1L]]$peptide_sequence, "MKS")
})

test_that("Sec and mitochondrial peptides consume T009 translation unchanged", {
  sec <- peptide_test_contexts(list(reconstruction_test_row("s", "chrS2", 4L, "A", "G")), "TS.1")[[1L]]
  p <- Neoantimon2:::peptide_enumerate(sec, peptide_test_policy(8L,9L))$peptides
  expect_identical(vapply(p, `[[`, "", "peptide_sequence"), c("MEPGKPGU","EPGKPGUP","MEPGKPGUP"))
  expect_identical(p[[1L]]$normal_correspondence$sequence, "MKPGKPGU")
  expect_identical(p[[1L]]$provenance$translation_exception_ids, list("GENCODE_SEC:TS.1:8"))
  expect_true("II" %in% unlist(p[[3L]]$classes))
  mito <- peptide_test_contexts(list(reconstruction_test_row("s", "chrM", 13L, "T", "C")), "TM2.1")[[1L]]
  expect_identical(peptide_test_sequences(mito,1L,2L), c("R","MR"))
  expect_identical(mito$translation$genetic_code_id, 2L)
})

test_that("exact self search is literal, overlapping, and never a filter", {
  path <- tempfile(); writeLines(c(">one|tx1", "AAAUAAAU", ">two|tx2", "ACDEFG"), path)
  x <- Neoantimon2:::peptide_self_search(c("AA","AU","ACD","YYY"),path,list(version="test"),"LITERAL_U")
  expect_true(x$AU$exact_match); expect_identical(vapply(x$AU$matches, `[[`, 0L, "start_aa"), c(3L,7L))
  expect_length(x$AA$matches, 4L); expect_true(x$ACD$exact_match); expect_false(x$YYY$exact_match)
  expect_null(Neoantimon2:::peptide_self_search("AU",path)$AU$exact_match)
  expect_identical(Neoantimon2:::peptide_self_search("ACD")$ACD$status,"UNAVAILABLE")
  expect_identical(x, Neoantimon2:::peptide_self_search(c("YYY","ACD","AU","AA"),path,list(version="test"),"LITERAL_U"))
})

test_that("presenter registry retains explicit DR DQ DP pairs and unresolved input", {
  raw <- c("a*02:01", "HLA-A*02:01:01N", "HLA-DRB1*04:01",
    "HLA-DRA*01:01~HLA-DRB3*01:01", "HLA-DRA*01:01~HLA-DRB4*01:03",
    "HLA-DRA*01:01~HLA-DRB5*01:01", "DQB1*02:01~DQA1*05:01",
    "DPA1*01:03~DPB1*04:01", "DQA1*05:01", "HLA-DQA1*05:01~HLA-DPB1*04:01")
  x <- Neoantimon2:::peptide_presenters(raw)
  expect_length(x,10L); expect_equal(sum(vapply(x,function(z) z$status=="RESOLVED",TRUE)),7L)
  resolved <- Filter(function(z) z$status=="RESOLVED",x)
  keys <- vapply(resolved,`[[`,"","mhc_key")
  expect_true("HLA-A*02:01:01N" %in% keys)
  expect_true("HLA-DQA1*05:01~HLA-DQB1*02:01" %in% keys)
  expect_true("HLA-DPA1*01:03~HLA-DPB1*04:01" %in% keys)
})

test_that("runtime publishes provenance-preserving candidates without fabricated predictions", {
  p <- peptide_test_run()
  snapshot <- Neoantimon2:::read_runtime_snapshot(p)
  artifact <- readRDS(Neoantimon2:::peptide_stage_path(p))
  expect_gt(length(snapshot$peptide_origins),0L); expect_gt(length(snapshot$neoantigen_results),0L)
  expect_silent(Neoantimon2:::validate_snapshot(snapshot, Neoantimon2:::peptide_contract()))
  expect_error(Neoantimon2:::validate_snapshot(snapshot), "version|enum")
  expect_equal(nrow(neo_results(p)),length(snapshot$neoantigen_results))
  expect_equal(nrow(neo_origins(p)),length(snapshot$peptide_origins))
  expect_identical(neo_results(neo_project(path=p$path)),neo_results(p))
  expect_equal(nrow(artifact$prediction_evidence),0L)
  for (c in snapshot$neoantigen_results) {
    expect_null(c$epitope_key); expect_null(c$binding_core); expect_null(c$mhc_rank_pct)
    expect_identical(c$prediction_status,"NOT_RUN")
    expect_true(length(neo_explain(p,candidate_id=c$candidate_id)$candidates)==1L)
  }
  expected <- Neoantimon2:::peptide_construct(p, readRDS(Neoantimon2:::translation_stage_path(p)), snapshot$events)
  expect_identical(vapply(expected$candidates,`[[`,"","candidate_id"),vapply(artifact$candidates,`[[`,"","candidate_id"))
  fp <- Neoantimon2:::peptide_stage_fingerprints(p)
  q <- p; q$config$peptides$class_i_lengths <- list(4L)
  expect_false(identical(fp,Neoantimon2:::peptide_stage_fingerprints(q)))
  q <- p; q$hla <- list("patient-1"=list("HLA-A*01:01"))
  expect_false(identical(fp,Neoantimon2:::peptide_stage_fingerprints(q)))
})

test_that("Sec reaches both candidate classes with stable literal sequence identities", {
  p <- peptide_test_run(TRUE); artifact <- readRDS(Neoantimon2:::peptide_stage_path(p))
  expect_gt(length(artifact$candidates),0L)
  expect_setequal(vapply(artifact$candidates,`[[`,"","hla_class"),c("I","II"))
  for (c in artifact$candidates) {
    expect_match(c$ligand_sequence,"U",fixed=TRUE)
    expect_null(c$mhc_raw_score); expect_null(c$binding_core); expect_null(c$epitope_key)
    expect_identical(c$candidate_key,Neoantimon2:::biological_id("candidate_key",
      c[c("species","hla_class","ligand_sequence","mhc_key")]))
  }
  for (pep in artifact$peptides) expect_identical(pep$provenance$sec_sites_used,list("GENCODE_SEC:TS.1:8"))
  expect_silent(Neoantimon2:::validate_snapshot(Neoantimon2:::read_runtime_snapshot(p),Neoantimon2:::peptide_contract()))
})

test_that("same sequence in different biological origins remains represented", {
  p <- peptide_test_run(); snapshot <- Neoantimon2:::read_runtime_snapshot(p)
  contexts <- readRDS(Neoantimon2:::translation_stage_path(p))
  x <- contexts[[1L]]; y <- x; y$transcript_id <- "OTHER.1"
  y$reconstruction_id <- "reconstruction:other"; y$translation$translation_id <- "translation:other"
  a <- Neoantimon2:::peptide_construct(p,list(x,y),snapshot$events)
  keys <- vapply(a$origins,`[[`,"","sample_peptide_id")
  expect_true(anyDuplicated(keys)>0L)
  expect_false(anyDuplicated(vapply(a$origins,`[[`,"","origin_id"))>0L)
  expect_true(all(vapply(a$candidates,function(c) c$n_origins==2L,TRUE)))
  expect_identical(sort(vapply(a$candidates,`[[`,"","candidate_id")),
    sort(vapply(Neoantimon2:::peptide_construct(p,list(y,x),snapshot$events)$candidates,`[[`,"","candidate_id")))
  q <- p; q$hla <- list("patient-1"=list("HLA-A*01:01"))
  changed_hla <- Neoantimon2:::peptide_construct(q,list(x,y),snapshot$events)
  expect_identical(vapply(a$origins,`[[`,"","origin_id"),vapply(changed_hla$origins,`[[`,"","origin_id"))
  expect_false(identical(vapply(a$candidates,`[[`,"","candidate_id"),vapply(changed_hla$candidates,`[[`,"","candidate_id")))
})

test_that("mapping distinguishes shifted counterpart and reference spans", {
  normal <- "MACDEFGHIKLMNPQRSTVWY"; mutant <- "MADEFGHIALMNPQRSTVWY"
  intervals <- list(list(normal_aa_start=3L,normal_aa_end=3L,mutant_aa_start=3L,mutant_aa_end=2L,
    normal_sequence="C",mutant_sequence=""),
    list(normal_aa_start=10L,normal_aa_end=10L,mutant_aa_start=9L,mutant_aa_end=9L,
      normal_sequence="K",mutant_sequence="A"))
  map <- Neoantimon2:::peptide_map(normal,mutant,intervals)
  x <- Neoantimon2:::peptide_correspondence(normal,map,intervals,8L,10L,"unchanged","IAL")
  expect_identical(x$status,"INDEL_SHIFTED"); expect_identical(x$sequence,"IKL")
  expect_identical(x$reference_start_aa,9L)
  ct <- peptide_test_long(); ct$translation$frame_effect <- "shifted_unrestored"
  ct$translation$frameshift_start_aa <- 25L
  early <- Neoantimon2:::peptide_enumerate(ct,peptide_test_policy(3L,3L))$peptides
  expect_true(all(vapply(early,function(x)!is.null(x$normal_correspondence$sequence),TRUE)))
  broken <- intervals; broken[[2L]]$mutant_sequence <- "Y"
  expect_error(Neoantimon2:::peptide_map(normal,mutant,broken),"mapping")
})

test_that("configuration and cached inputs cannot silently reuse a different policy", {
  for (bad in list(0L,-1L,2.5,NA_real_,Inf,"9"))
    expect_error(neo_config(peptides=list(class_i_lengths=bad)),"length")
  expect_identical(Neoantimon2:::peptide_policy(list(class_i_lengths=c(11L,9L,9L)))$class_i_lengths,list(9L,11L))
  p <- peptide_test_run(); before <- Neoantimon2:::peptide_stage_fingerprints(p)
  q <- p
  ix <- which(vapply(q$resources,function(x)x$logical_name=="reference.protein_fasta",TRUE))
  file <- tempfile(); writeLines(c(">different","MACDEFGHIK"),file)
  q$resources[[ix]]$resolved_path <- file
  expect_false(identical(before$dependency_fingerprint,Neoantimon2:::peptide_stage_fingerprints(q)$dependency_fingerprint))
  expect_error(Neoantimon2:::checkpoint_reusable(q,"peptide_candidates",Neoantimon2:::peptide_stage_fingerprints(q)),"fingerprint")
  expect_identical(Neoantimon2:::peptide_stage_run(p,readRDS(Neoantimon2:::translation_stage_path(p)),
    Neoantimon2:::read_runtime_snapshot(p)$events),readRDS(Neoantimon2:::peptide_stage_path(p)))
})
