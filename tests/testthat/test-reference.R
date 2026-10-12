ref_project <- function(legacy = FALSE, mutate = NULL, checksum = "sha256") {
  root <- tempfile("reference-"); dir.create(root)
  files <- list.files(testthat::test_path("..", "fixtures", "reference"), full.names = TRUE)
  file.copy(files, root)
  if (!is.null(mutate)) mutate(root)
  roles <- if (legacy) c(refFlat = "legacy.refFlat", refMrna = "legacy.fa") else
    c(genome_fasta = "genome.fa", genome_fai = "genome.fa.fai", annotation_gtf = "annotation.gtf",
      transcript_fasta = "transcripts.fa", protein_coding_transcript_fasta = "coding.fa",
      protein_fasta = "proteins.fa", selenocysteine_metadata = "sec.tsv",
      mane_summary = "mane.tsv", mane_gtf = "mane.gtf")
  res <- lapply(names(roles), function(role) list(path = file.path(root, roles[[role]]),
    required = TRUE, resource_id = if (role == "annotation_gtf") Neoantimon2:::reference_scope_policy() else NULL, version = if (legacy) NULL else if (startsWith(role, "genome")) "GRCh38.p14"
      else if (startsWith(role, "mane")) "1.5" else "50", source = "synthetic"))
  names(res) <- paste0("reference.", names(roles))
  neo_project(config = neo_config(env = character(), resources = res, checksum = checksum,
    scientific = list(genome_assembly = if (legacy) NULL else "GRCh38.p14")))
}
ref_edit <- function(root, file, pattern, replacement) {
  path <- file.path(root, file)
  writeLines(gsub(pattern, replacement, readLines(path), fixed = TRUE), path)
}
ref_load <- function(p) Neoantimon2:::reference_bundle(p, retain_indexes = TRUE)

test_that("modern metadata is deterministic, complete and separate from runtime", {
  p <- ref_project(); b <- ref_load(p)
  expect_identical(b$family, "gencode")
  expect_identical(b, ref_load(p))
  expect_equal(b$genome$length, c(120, 60))
  expect_equal(nrow(b$transcripts), 3)
  expect_equal(nrow(b$proteins), 3)
  expect_equal(nrow(b$mane), 1)
  expect_true("ENST000002.1" %in% b$proteins$transcript_id)
  expect_true("1" %in% Neoantimon2:::reference_transcript_features(b, "ENST000001.2")$phase)
  expect_equal(b$sec$transcript_id, "ENST000001.2")
  expect_equal(b$mitochondrial_contigs, "chrM")
  expect_equal(b$coordinate_convention, "1-based inclusive")
  expect_false(any(c("sequence", "sequences") %in% names(p)))
  expect_true(all(b$checksum_state == "observed_sha256"))
  expect_identical(sort(getNamespaceExports("Neoantimon2")), sort(current_public_api))
  expect_identical(Neoantimon2:::schema_contract()$schema_version, "1.0.0")
  expect_identical(Neoantimon2:::schema_contract()$id_schema_version, "1.0.0")
})

test_that("legacy coordinates and unknown capabilities remain explicit", {
  b <- ref_load(ref_project(TRUE))
  expect_identical(b$family, "legacy")
  expect_identical(b$coordinate_convention, "0-based half-open")
  expect_equal(b$annotation$txStart, c(0, 60))
  expect_equal(b$annotation$exonStarts[[2]], c(60, 80))
  expect_equal(b$transcripts$id, c("NM_000001.2", "NM_000002.1"))
  expect_false(b$capabilities$protein_reference_available)
  expect_false(b$capabilities$mane_available)
  expect_true(is.na(b$capabilities$assembly_verified))
  expect_null(b$declared$assembly)
})

test_that("malformed modern resources fail without aliases or ID rewriting", {
  cases <- list(
    c("genome.fa.fai", "chr1\t120", "chr1\t0"),
    c("genome.fa.fai", "chrM", "chr1"),
    c("genome.fa.fai", "\t121", "\t0"),
    c("annotation.gtf", "chr1", "1"),
    c("annotation.gtf", "chr1", "chrMissing"),
    c("annotation.gtf", "\t+\t", "\t?\t"),
    c("annotation.gtf", "\t1\t12\t", "\t13\t12\t"),
    c("annotation.gtf", "\t0\tgene_id", "\t3\tgene_id"),
    c("annotation.gtf", 'gene_id "ENSG000001.3";', 'gene_id ENSG000001.3;'),
    c("transcripts.fa", "ENST000002.1", "ENST000001.2"),
    c("transcripts.fa", "ENST000001.2", "ENST000099.1"),
    c("proteins.fa", "ENST000001.2", "ENST000099.1"),
    c("mane.tsv", "ENST000001.2", "ENST000099.1"),
    c("mane.tsv", "MANE Select", "garbage"),
    c("sec.tsv", "ENST000001.2", "ENST000099.1"),
    c("annotation.gtf", "version 50", "version 49"))
  for (case in cases) expect_error(ref_load(ref_project(mutate = function(root)
    ref_edit(root, case[1], case[2], case[3]))), info = paste(case, collapse = "/"))
  expect_error(ref_project(mutate = function(root) unlink(file.path(root, "mane.tsv"))))
  expect_error(ref_load(ref_project(mutate = function(root)
    writeLines("chr1\t120", file.path(root, "genome.fa.fai")))))
})

test_that("legacy malformed records and missing sequence matches fail", {
  cases <- list(c("legacy.refFlat", "\t2\t", "\t3\t"),
    c("legacy.refFlat", "0,20,", "0,"), c("legacy.refFlat", "12,35,", "12,"),
    c("legacy.refFlat", "\t+\t", "\t?\t"), c("legacy.refFlat", "\t0\t35\t", "\t40\t35\t"),
    c("legacy.refFlat", "\t60\t95\t", "\tx\t95\t"),
    c("legacy.refFlat", "NM_000002.1", "NM_000001.2"),
    c("legacy.fa", "NM_000001.2", "ENST000001.2"),
    c("legacy.fa", "NM_000002.1", "NM_000001.2"))
  for (case in cases) expect_error(ref_load(ref_project(TRUE, function(root)
    ref_edit(root, case[1], case[2], case[3]))))
  expect_error(ref_load(ref_project(TRUE, function(root)
    writeLines("wrong\tfield\tcount", file.path(root, "legacy.refFlat")))))
})

test_that("reference identity uses frozen T003 resource checksums", {
  p <- ref_project(); b <- ref_load(p)
  role <- which(vapply(p$resources, function(r) r$logical_name == "reference.protein_fasta", TRUE))
  r <- p$resources[[role]]
  cfg <- p$config; cfg$resources[[r$logical_name]]$checksum <- r$integrity$observed
  verified <- ref_load(neo_project(config = cfg))
  expect_true("verified_sha256" %in% verified$checksum_state)
  cfg$resources[[r$logical_name]]$checksum$value <- paste(rep("0", 64), collapse = "")
  expect_error(neo_project(config = cfg), "checksum mismatch")
  cfg$checksum <- "metadata"
  expect_error(ref_load(neo_project(config = cfg)), "checksum mismatch")
  p2 <- ref_project(mutate = function(root) ref_edit(root, "proteins.fa", "MUAAAAAA", "MUAAAAAG"))
  expect_false(identical(p$analysis_id, p2$analysis_id))
  expect_identical(neo_project()$analysis_id, neo_project()$analysis_id)
})

test_that("dictionary comparison requires exact names and lengths", {
  b <- ref_load(ref_project()); compare <- Neoantimon2:::reference_compare_dictionary
  expect_true(compare(b$genome, b$genome)$compatible)
  expect_false(compare(b$genome[1, , drop = FALSE], b$genome)$compatible)
  expect_equal(compare(b$genome[1, , drop = FALSE], b$genome)$missing_from_supplied, "chrM")
  x <- b$genome; x$contig[1] <- "1"
  expect_false(compare(x, b$genome)$compatible)
  expect_equal(compare(x, b$genome)$alias_candidates, "1 -> chr1")
  x <- b$genome; x$length[1] <- 119
  expect_false(compare(x, b$genome)$compatible)
  expect_false(compare(data.frame(contig = "absent", length = 1), b$genome)$compatible)
})

ref_scope_fixture <- function(root) {
  append <- function(from, to) cat(paste0(paste(readLines(file.path(root, from)), collapse = "\n"), "\n"),
    file = file.path(root, to), append = TRUE)
  append("scope-extra.gtf", "mane.gtf")
  append("scope-extra-transcripts.fa", "transcripts.fa")
  append("scope-extra-transcripts.fa", "coding.fa")
  append("scope-extra-proteins.fa", "proteins.fa")
  append("scope-extra-mane.tsv", "mane.tsv")
  cat("ENST000004.1\t2\n", file = file.path(root, "sec.tsv"), append = TRUE)
}

ref_plus_clinical_fixture <- function(root) {
  ref_scope_fixture(root)
  clone_lines <- function(file, select_id, plus_id, select_protein, plus_protein,
                          select_refseq, plus_refseq, select_refseq_protein, plus_refseq_protein) {
    path <- file.path(root, file)
    lines <- readLines(path)
    positions <- which(grepl(select_id, lines, fixed = TRUE))
    if (endsWith(file, ".fa")) positions <- as.vector(rbind(positions, positions + 1L))
    selected <- lines[positions]
    selected <- gsub(select_id, plus_id, selected, fixed = TRUE)
    selected <- gsub(select_protein, plus_protein, selected, fixed = TRUE)
    selected <- gsub(select_refseq, plus_refseq, selected, fixed = TRUE)
    selected <- gsub(select_refseq_protein, plus_refseq_protein, selected, fixed = TRUE)
    selected <- gsub("MANE Select", "MANE Plus Clinical", selected, fixed = TRUE)
    cat(paste0(selected, "\n"), sep = "", file = path, append = TRUE)
  }
  for (file in c("annotation.gtf", "mane.gtf", "transcripts.fa", "coding.fa", "proteins.fa", "mane.tsv"))
    clone_lines(file, "ENST000001.2", "ENST000007.1", "ENSP000001.4", "ENSP000007.1",
      "NM_000001.2", "NM_000007.1", "NP_000001.2", "NP_000007.1")
  for (file in c("mane.gtf", "transcripts.fa", "coding.fa", "proteins.fa", "mane.tsv"))
    clone_lines(file, "ENST000004.1", "ENST000008.1", "ENSP000004.1", "ENSP000008.1",
      "NM_000004.1", "NM_000008.1", "NP_000004.1", "NP_000008.1")
}

test_that("Plus Clinical remains reconstructable in primary scope without becoming representative", {
  b <- ref_load(ref_project(mutate = ref_plus_clinical_fixture))
  select <- b$mane[b$mane$Ensembl_nuc == "ENST000001.2", ]
  plus <- b$mane[b$mane$Ensembl_nuc == "ENST000007.1", ]
  patch_plus <- b$mane[b$mane$Ensembl_nuc == "ENST000008.1", ]
  expect_equal(nrow(b$mane), 5)
  expect_identical(select$MANE_status, "MANE Select")
  expect_true(select$representative_eligible)
  expect_identical(plus$Ensembl_Gene, select$Ensembl_Gene)
  expect_identical(plus$MANE_status, "MANE Plus Clinical")
  expect_identical(plus$scope, "IN_SCOPE_PRIMARY")
  expect_true(plus$reconstruction_eligible)
  expect_false(plus$representative_eligible)
  expect_true("ENST000002.1" %in% b$annotation$transcript_id)
  expect_false("ENST000002.1" %in% b$mane$Ensembl_nuc)
  expect_identical(patch_plus$MANE_status, "MANE Plus Clinical")
  expect_identical(patch_plus$scope, "OUT_OF_SCOPE_REFERENCE_REGION")
  expect_identical(patch_plus$scope_reason, "MANE_PATCH_SEQUENCE")
  expect_false(patch_plus$reconstruction_eligible)
  expect_false(patch_plus$representative_eligible)
  counts <- Neoantimon2:::reference_summary(b)$counts
  expect_equal(unlist(counts[c("mane_records", "mane_select", "mane_plus_clinical",
    "mane_in_scope_select", "mane_in_scope_plus_clinical",
    "mane_out_of_scope_select", "mane_out_of_scope_plus_clinical")]),
    c(mane_records = 5, mane_select = 3, mane_plus_clinical = 2,
      mane_in_scope_select = 1, mane_in_scope_plus_clinical = 1,
      mane_out_of_scope_select = 2, mane_out_of_scope_plus_clinical = 1))
})

test_that("PRI reconstruction and ALL catalogs preserve every scope explicitly", {
  p <- ref_project(mutate = ref_scope_fixture); b <- ref_load(p)
  expect_equal(b$profile$coordinate_scope, "primary_assembly")
  expect_equal(b$profile$sequence_catalog_scope, "all_regions")
  expect_equal(nrow(b$mane), 3)
  expect_equal(sum(b$mane$scope == "IN_SCOPE_PRIMARY"), 1)
  expect_equal(sum(b$mane$scope == "OUT_OF_SCOPE_REFERENCE_REGION"), 2)
  expect_setequal(b$mane$scope_reason, c("PRIMARY_SEQUENCE", "MANE_PATCH_SEQUENCE", "MANE_ALTERNATE_LOCUS"))
  expect_equal(sum(b$mane$reconstruction_eligible), 1)
  expect_equal(sum(b$transcripts$reconstruction_scope == "PRIMARY_RECONSTRUCTION_ELIGIBLE"), 3)
  expect_equal(sum(b$transcripts$reconstruction_scope == "CATALOG_ONLY_OUT_OF_SCOPE"), 2)
  expect_equal(nrow(b$proteins), 5)
  expect_true(all(b$proteins$self_catalog_member))
  expect_equal(b$mane$seqname[b$mane$Ensembl_nuc == "ENST000004.1"], "chr1_GL000001v1_fix")
  expect_false("ENST000002.1" %in% b$mane$Ensembl_nuc)
  expect_true("ENST000002.1" %in% b$annotation$transcript_id)
  expect_equal(b$sec$reconstruction_scope[b$sec$transcript_id == "ENST000004.1"], "CATALOG_ONLY_OUT_OF_SCOPE")
  report <- Neoantimon2:::reference_summary(b)
  expect_equal(report$unexplained_mane_mismatches, 0)
  expect_equal(report$counts$mane_out_of_scope_feature_rows, 6)
  expect_equal(report$mane_reasons[[1]]$contigs, 1)
  expect_equal(report$mane_reasons[[2]]$contigs, 1)
  # Input metadata row order is irrelevant, while bound file content identity changes.
  for (role in c("annotation_gtf", "mane_gtf", "selenocysteine_metadata")) {
    file <- b$resources[[role]]$resolved_path; writeLines(rev(readLines(file)), file)
  }
  file <- b$resources$mane_summary$resolved_path; lines <- readLines(file)
  writeLines(c(lines[1], rev(lines[-1])), file)
  p2 <- neo_project(config = p$config); b2 <- ref_load(p2)
  for (n in c("annotation", "mane", "mane_annotation", "sec", "transcripts", "proteins"))
    expect_identical(b[[n]], b2[[n]], info = n)
  expect_false(identical(p$analysis_id, p2$analysis_id))
})

test_that("missing MANE contigs need corroborating metadata, never guessed aliases", {
  mutations <- list(c("mane.gtf", "chr1", "1"),
    c("mane.gtf", "chr1_GL000001v1_fix", "unknown"),
    c("mane.tsv", "NW_000001.1", "NC_000001.11"),
    c("mane.tsv", "\t1\t27\t+", "\t2\t27\t+"))
  for (m in mutations) expect_error(ref_load(ref_project(mutate = function(root) {
    ref_scope_fixture(root); ref_edit(root, m[1], m[2], m[3])
  })))
  expect_error(ref_load(ref_project(mutate = function(root) {
    ref_scope_fixture(root); ref_edit(root, "annotation.gtf", "chr1", "chr1_GL000001v1_fix")
  })), "GTF contig absent")
})

test_that("reference scope policy is included in frozen scientific identity", {
  p <- ref_project(); cfg <- p$config
  cfg$resources$reference.annotation_gtf$resource_id <- "future_all_regions_v2"
  changed <- neo_project(config = cfg)
  expect_false(identical(p$analysis_id, changed$analysis_id))
  expect_error(ref_load(changed), "scope policy")
  cfg$resources$reference.annotation_gtf$resource_id <- NULL
  expect_error(ref_load(neo_project(config = cfg)), "scope policy")
})

test_that("remaining file, provenance and legacy separation failures are explicit", {
  p <- ref_project(); cfg <- p$config
  cfg$scientific$genome_assembly <- "GRCh37"
  expect_error(ref_load(neo_project(config = cfg)), "assembly conflict")
  cfg <- p$config; cfg$resources$reference.mane_gtf$version <- "1.4"
  expect_error(ref_load(neo_project(config = cfg)), "release conflict")
  cfg <- p$config; cfg$resources$reference.refFlat <- cfg$resources$reference.annotation_gtf
  expect_error(ref_load(neo_project(config = cfg)), "Mixed modern/legacy")
  cfg <- p$config; cfg$resources$reference.protein_fasta <- NULL
  expect_error(ref_load(neo_project(config = cfg)), "reference role")
  cfg <- p$config; cfg$resources$reference.protein_fasta$type <- "directory"
  expect_error(neo_project(config = cfg), "type mismatch")
  r <- p$resources[[1]]; cat("changed\n", file = r$resolved_path, append = TRUE)
  expect_error(ref_load(p), "changed since project")
  expect_error(ref_load(ref_project(mutate = function(root)
    ref_edit(root, "transcripts.fa", "ATGTGA", "ATGZGA"))), "FASTA sequence")
  p <- ref_project(TRUE); cfg <- p$config; cfg$scientific$genome_assembly <- "GRCh37"
  b <- ref_load(neo_project(config = cfg))
  expect_identical(b$declared$assembly, "GRCh37")
  expect_true(is.na(b$capabilities$assembly_verified))
})

test_that("GTF duplicate attributes and unsupported compression fail", {
  expect_error(ref_load(ref_project(mutate = function(root)
    ref_edit(root, "annotation.gtf", 'gene_id "ENSG000001.3";',
      'gene_id "ENSG000001.3"; gene_id "ENSG000099.1";'))), "Duplicate GTF attribute")
  expect_error(ref_load(ref_project(mutate = function(root)
    ref_edit(root, "genome.fa.fai", "\t120\t121", "\t119\t120"))), "line-width")
  p <- ref_project(); cfg <- p$config
  path <- cfg$resources$reference.annotation_gtf$path
  gz <- paste0(path, ".gz"); con <- gzfile(gz, "wt"); writeLines(readLines(path), con); close(con)
  cfg$resources$reference.annotation_gtf$path <- gz
  expect_equal(nrow(ref_load(neo_project(config = cfg))$annotation), 3)
  file.copy(path, gz, overwrite = TRUE)
  expect_error(ref_load(neo_project(config = cfg)), "Compression/extension")
})

test_that("streaming FASTA handles record boundaries and retains exact headers", {
  file <- tempfile(); writeLines(c(">NM_1.2 exact header", rep("A", 20001), ">NM_2.3", "CC"), file)
  x <- Neoantimon2:::reference_fasta(file, "legacy")
  expect_equal(x$length, c(20001, 2))
  expect_equal(x$raw_header, c("NM_1.2 exact header", "NM_2.3"))
  writeLines(c(">NM_1", ">NM_2", "A"), file)
  expect_error(Neoantimon2:::reference_fasta(file, "legacy"), "Empty FASTA")
})

test_that("MANE noncoding records use RNA traceability without fabricated proteins", {
  p <- ref_project(mutate = function(root) {
    a <- 'gene_id "ENSG000006.1"; transcript_id "ENST000006.1"; gene_type "snRNA"; transcript_type "snRNA";'
    line <- paste("chr1", "synthetic", "transcript", 100, 110, ".", "+", ".", a, sep = "\t")
    for (file in c("annotation.gtf", "mane.gtf")) cat(line, "\n", sep = "", file = file.path(root, file), append = TRUE)
    cat(">ENST000006.1|ENSG000006.1|-|-|RNA-201|RNA|11|snRNA|\nACGTACGTACG\n",
      file = file.path(root, "transcripts.fa"), append = TRUE)
    cat("GeneID:6\tENSG000006.1\t\tRNA\tsynthetic RNA\tNR_000006.1\t\tENST000006.1\t\tMANE Select\tNC_000001.11\t100\t110\t+\n",
      file = file.path(root, "mane.tsv"), append = TRUE)
  })
  b <- ref_load(p)
  row <- b$mane[b$mane$Ensembl_nuc == "ENST000006.1", ]
  expect_true(row$representative_eligible)
  expect_false(row$protein_reference_available)
  expect_equal(nrow(b$proteins), 3)
  expect_equal(nrow(b$transcripts), 4)
})

test_that("unbound content and conflicting provenance never become validated references", {
  expect_error(ref_load(ref_project(checksum = "none")), "metadata identity")
  expect_error(ref_load(ref_project(mutate = function(root)
    ref_edit(root, "annotation.gtf", "##provider: GENCODE", "##provider: Other"))), "provider conflict")
  expect_error(ref_load(ref_project(mutate = function(root)
    ref_edit(root, "annotation.gtf", "(GRCh38)", "(GRCh37)"))), "assembly conflict")
  expect_error(ref_load(ref_project(mutate = function(root)
    ref_edit(root, "annotation.gtf", "stop_codon\t33\t35", "stop_codon\t36\t38"))), "outside transcript bounds")
  p <- ref_project(); b <- ref_load(p)
  expect_equal(b$annotation$strand[b$annotation$transcript_id == "ENST000002.1"], "-")
  f <- Neoantimon2:::reference_transcript_features(b, "ENST000001.2")
  expect_equal(sum(f$feature == "exon"), 2)
  expect_true(all(c("CDS", "stop_codon", "Selenocysteine") %in% f$feature))
  expect_equal(f$start[f$feature == "exon"], c(1, 21))
  expect_equal(f$end[f$feature == "exon"], c(12, 35))
})


test_that("extra trailing empty fields are malformed, not silently truncated", {
  for (file in c("genome.fa.fai", "annotation.gtf", "mane.tsv"))
    expect_error(ref_load(ref_project(mutate = function(root) {
      path <- file.path(root, file); x <- readLines(path)
      row <- which(!startsWith(x, "#"))[1]; x[row] <- paste0(x[row], "\t")
      writeLines(x, path)
    })), "field count")
  expect_error(ref_load(ref_project(TRUE, function(root) {
    path <- file.path(root, "legacy.refFlat"); writeLines(paste0(readLines(path), "\t"), path)
  })), "field count")
})


test_that("MANE GRCh38 metadata cannot be declared as a GRCh37 modern bundle", {
  p <- ref_project(); cfg <- p$config
  cfg$scientific$genome_assembly <- "GRCh37"
  cfg$resources$reference.genome_fasta$version <- "GRCh37"
  cfg$resources$reference.genome_fai$version <- "GRCh37"
  expect_error(ref_load(neo_project(config = cfg)), "MANE summary requires GRCh38")
})


test_that("default bundle is a compact descriptor with the same validated identity", {
  p <- ref_project(mutate = ref_scope_fixture)
  full <- ref_load(p); compact <- Neoantimon2:::reference_bundle(p)
  expect_false(compact$indexes_materialized)
  expect_null(compact$transcripts)
  expect_null(compact$annotation)
  expect_identical(compact$analysis_id, full$analysis_id)
  expect_identical(compact$resources, full$resources)
  expect_identical(compact$mane, full$mane)
  expect_identical(compact$sec, full$sec)
  expect_identical(Neoantimon2:::reference_summary(compact), Neoantimon2:::reference_summary(full))
  expect_identical(Neoantimon2:::reference_transcript_features(compact, "ENST000001.2"),
    Neoantimon2:::reference_transcript_features(full, "ENST000001.2"))
  expect_error(Neoantimon2:::reference_bundle(p, retain_indexes = NA), "retain_indexes")
})
