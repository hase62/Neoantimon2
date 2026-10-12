event_fixture <- function(name) normalizePath(testthat::test_path("..", "fixtures", "events", name))

event_test_reference <- function() {
  path <- event_fixture("repeat.fa")
  list(con = file(path, "rb"), fai = Neoantimon2:::reference_fai(event_fixture("repeat.fa.fai"), path))
}

event_test_project <- function(input_names, roles, ids = basename(input_names)) {
  root <- tempfile("t006-"); dir.create(root)
  files <- list.files(testthat::test_path("..", "fixtures", "reference"), full.names = TRUE)
  file.copy(files, root)
  ref_roles <- c(genome_fasta = "genome.fa", genome_fai = "genome.fa.fai",
    annotation_gtf = "annotation.gtf", transcript_fasta = "transcripts.fa",
    protein_coding_transcript_fasta = "coding.fa", protein_fasta = "proteins.fa",
    selenocysteine_metadata = "sec.tsv", mane_summary = "mane.tsv", mane_gtf = "mane.gtf")
  resources <- lapply(names(ref_roles), function(role) list(path = file.path(root, ref_roles[[role]]),
    required = TRUE, resource_id = if (role == "annotation_gtf") Neoantimon2:::reference_scope_policy() else NULL,
    version = if (startsWith(role, "genome")) "GRCh38.p14" else
      if (startsWith(role, "mane")) "1.5" else "50", source = "synthetic"))
  names(resources) <- paste0("reference.", names(ref_roles))
  inputs <- unname(Map(function(name, role, id) list(input_id = id, input_role = role,
    sample_id = "sample-1", path = if (file.exists(name)) normalizePath(name) else event_fixture(name)),
    input_names, roles, ids))
  config <- neo_config(env = character(), resources = resources, checksum = "sha256",
    paths = list(run_root = file.path(root, "runs")),
    scientific = list(genome_assembly = "GRCh38.p14"))
  neo_project(samples = data.frame(sample_id = "sample-1", patient_id = "patient-1"),
              inputs = inputs, config = config)
}

test_that("FAI coordinate conversion and independently specified normalizations", {
  ref <- event_test_reference(); on.exit(close(ref$con))
  expect_identical(Neoantimon2:::event_reference_fetch(ref, "chrR", 1, 3), "CCC")
  expect_identical(Neoantimon2:::event_reference_fetch(ref, "chrR", 15, 16), "GG")
  expect_error(Neoantimon2:::event_reference_fetch(ref, "chrR", 0, 1), "outside")
  gold <- read.delim(event_fixture("gold_cases.tsv"), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(gold))) {
    x <- gold[i, ]
    y <- Neoantimon2:::event_normalize(ref, x$chrom, x$pos, x$ref, x$alt)
    expect_equal(y$pos, x$expected_pos, info = x$case)
    expect_identical(y$ref, x$expected_ref, info = x$case)
    expect_identical(y$alt, x$expected_alt, info = x$case)
    z <- Neoantimon2:::event_normalize(ref, y$contig, y$pos, y$ref, y$alt)
    expect_identical(z[c("pos", "ref", "alt")], y[c("pos", "ref", "alt")])
  }
  a <- Neoantimon2:::event_normalize(ref, "chrR", 5, "A", "AA")
  b <- Neoantimon2:::event_normalize(ref, "chrR", 8, "A", "AA")
  expect_identical(a[c("pos", "ref", "alt")], b[c("pos", "ref", "alt")])
  expect_error(Neoantimon2:::event_normalize(ref, "chrR", 5, "T", "C"), "REF mismatch")
  expect_error(Neoantimon2:::event_normalize(ref, "chrR", 5, "A", "A"), "No-op")
  expect_error(Neoantimon2:::event_normalize(ref, "chrR", 5, "A", "N"), "Malformed")
  replacement <- Neoantimon2:::event_normalize(ref, "chrR", 4L, "A", "CT")
  expect_identical(replacement$event_type, "INDEL")
  expect_identical(replacement[c("pos", "ref", "alt")], list(pos = 4L, ref = "A", alt = "CT"))
})

test_that("VCF events retain provenance, allelic projections, and raw phase tokens", {
  p <- neo_run(event_test_project("somatic.vcf", "somatic_vcf"))
  rows <- Neoantimon2:::read_runtime_snapshot(p)$events
  expect_length(rows, 6L)
  expect_equal(nrow(neo_events(p)), 6L)
  expect_equal(nrow(neo_origins(p)), 0L)
  expect_equal(nrow(neo_results(p)), 0L)
  expect_equal(nrow(neo_multisample(p)), 0L)
  expect_equal(sum(vapply(rows, `[[`, "", "event_type") == "MNV"), 1L)
  multi <- Filter(function(x) x$event_details_json$source_records[[1]]$original$id == "multi", rows)
  expect_length(multi, 2L)
  expect_setequal(vapply(multi, function(x) x$event_details_json$source_records[[1]]$selected_alt_index, 0L), c(1L, 2L))
  imported <- jsonlite::read_json(file.path(p$path, "intermediate", "event_import.json"),
    simplifyVector = FALSE)
  imported_multi <- Filter(function(x) x$event_details_json$source_records[[1]]$original$id == "multi", imported)
  for (x in imported_multi) {
    src <- x$event_details_json$source_records[[1]]
    expect_identical(src$raw_format$sample, "1|2:10,3,4:0,1,2,3,4,5")
    expect_identical(src$original$alts, list("C", "G"))
    expect_identical(src$projected_info$AF, c("0.2", "0.3")[[src$selected_alt_index]])
    expect_identical(src$projected_format$AD$alt, c("3", "4")[[src$selected_alt_index]])
    expect_null(src$projected_format$PL)
    expect_null(x$phase_set)
    expect_length(x$phase_links_json, 0L)
    expect_null(x$tumor_alt_count)
  }
  expect_identical(Neoantimon2:::phase_relative_events(rows, multi[[1]]$event_id,
    multi[[2]]$event_id)$relation, "OPPOSITE_HAPLOTYPE")
  expect_identical(Neoantimon2:::read_manifest(p$path)$status, "COMPLETED")
  expect_true("event_import" %in% vapply(Neoantimon2:::read_manifest(p$path)$checkpoints, `[[`, "", "step_id"))
})

test_that("duplicate callers merge deterministically and somatic/germline differ", {
  p1 <- event_test_project(c("somatic.vcf", "duplicate.vcf"), c("somatic_vcf", "somatic_vcf"), c("one", "two"))
  p2 <- event_test_project(c("duplicate.vcf", "somatic.vcf"), c("somatic_vcf", "somatic_vcf"), c("two", "one"))
  a <- Neoantimon2:::event_import_project(p1)
  b <- Neoantimon2:::event_import_project(p2)
  expect_identical(sort(vapply(a, `[[`, "", "event_key")), sort(vapply(b, `[[`, "", "event_key")))
  expect_length(a, 6L)
  duplicate <- Filter(function(x) x$pos == 1L, a)[[1L]]
  expect_length(duplicate$event_details_json$source_records, 2L)
  expect_setequal(vapply(duplicate$caller_evidence_json, `[[`, "", "caller"), c("CallerOne", "CallerTwo"))
  p3 <- event_test_project(c("somatic.vcf", "duplicate.vcf"), c("somatic_vcf", "germline_vcf"), c("one", "two"))
  c <- Neoantimon2:::event_import_project(p3)
  expect_length(c, 7L)
  same_locus <- Filter(function(x) x$pos == 1L, c)
  expect_setequal(vapply(same_locus, `[[`, "", "event_role"), c("somatic", "germline"))
  expect_length(unique(vapply(same_locus, `[[`, "", "event_key")), 2L)
})

test_that("GDC/TCGA MAF alleles anchor on genome and preserve source values", {
  p <- event_test_project("somatic.maf", "somatic_maf")
  rows <- Neoantimon2:::event_import_project(p)
  expect_length(rows, 3L)
  expect_setequal(vapply(rows, `[[`, "", "event_type"), c("SNV", "INS", "DEL"))
  src <- lapply(rows, function(x) x$event_details_json$source_records[[1L]])
  expect_true(any(vapply(src, function(x) x$original$Strand == "-", TRUE)))
  expect_true(any(vapply(src, function(x) x$original$Reference_Allele == "-", TRUE)))
})

test_that("source insertion conventions produce one canonical event", {
  write_source <- function(lines, ext) {
    path <- tempfile(fileext = ext); writeLines(lines, path); path
  }
  cases <- list(
    VCF = list("somatic_vcf", write_source(c("##fileformat=VCFv4.2",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
      "chr1\t5\t.\tA\tAT\t.\tPASS\t."), ".vcf")),
    ANNOVAR = list("somatic_annovar", write_source(c("Chr\tStart\tEnd\tRef\tAlt",
      "chr1\t5\t5\t-\tT"), ".tsv")),
    VEP = list("somatic_vep", write_source(c("Uploaded_variation\tLocation\tAllele\tConsequence",
      "ins\tchr1:6-5\tT\tframeshift_variant"), ".tsv")),
    MAF = list("somatic_maf", write_source(c(
      "Chromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tVariant_Type",
      "chr1\t5\t6\t-\t-\tT\tINS"), ".maf")))
  rows <- lapply(cases, function(x) Neoantimon2:::event_import_project(event_test_project(x[[2]], x[[1]]))[[1L]])
  facts <- c("contig", "pos", "ref", "alt", "event_type", "event_key")
  expect_true(all(vapply(rows, function(x) identical(x[facts], rows$VCF[facts]), TRUE)))
  expect_identical(rows$VCF[facts[1:5]], list(contig = "chr1", pos = 5L,
    ref = "A", alt = "AT", event_type = "INS"))
  expect_identical(rows$VEP$event_details_json$source_records[[1L]]$reference_allele_origin,
    "selected_T005_reference_not_source_column")
})

test_that("ANNOVAR native alleles and invalid forms are distinguished", {
  write_annovar <- function(row) {
    path <- tempfile(fileext = ".tsv")
    writeLines(c("Chr\tStart\tEnd\tRef\tAlt", row), path)
    path
  }
  valid <- c("chr1\t1\t1\tA\tC", "chr1\t5\t6\tAC\t-",
    "chr1\t5\t5\t-\tT", "chr1\t5\t5\tA\tCT")
  types <- c("SNV", "DEL", "INS", "INDEL")
  for (i in seq_along(valid)) {
    row <- Neoantimon2:::event_import_project(event_test_project(write_annovar(valid[[i]]),
      "somatic_annovar"))[[1L]]
    expect_identical(row$event_type, types[[i]])
  }
  invalid <- list(c("chr1\t6\t5\t-\tT", "Start = End"),
    c("chr1\t1\t1\tT\tC", "REF mismatch"),
    c("chr1\t1\t1\t0\tC", "Ref=0"))
  for (case in invalid) expect_error(Neoantimon2:::event_import_project(
    event_test_project(write_annovar(case[[1L]]), "somatic_annovar")), case[[2L]])
})

test_that("MAF standard coordinates, tumor ALTs, and Variant_Type are checked", {
  write_maf <- function(row) {
    path <- tempfile(fileext = ".maf")
    writeLines(c(paste(c("Chromosome", "Start_Position", "End_Position", "Reference_Allele",
      "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Variant_Type"), collapse = "\t"), row), path)
    path
  }
  valid <- list(c("chr1\t1\t1\tA\tA\tC\tSNP", "SNV"),
    c("chr1\t5\t6\tAC\tAC\t-\tDEL", "DEL"),
    c("chr1\t5\t6\t-\t-\tT\tINS", "INS"))
  for (case in valid) expect_identical(Neoantimon2:::event_import_project(
    event_test_project(write_maf(case[[1L]]), "somatic_maf"))[[1L]]$event_type, case[[2L]])
  two <- Neoantimon2:::event_import_project(event_test_project(
    write_maf("chr1\t1\t1\tA\tC\tG\tSNP"), "somatic_maf"))
  expect_length(two, 2L)
  expect_setequal(vapply(two, `[[`, "", "alt"), c("C", "G"))
  for (row in c("chr1\t5\t6\t-\t-\tT\tDEL", "chr1\t5\t6\tAC\tAC\t-\tINS"))
    expect_error(Neoantimon2:::event_import_project(event_test_project(
      write_maf(row), "somatic_maf")), "Variant_Type contradicts")
  for (row in c("chr1\t6\t5\t-\t-\tT\tINS", "chr1\t5\t5\t-\t-\tT\tINS"))
    expect_error(Neoantimon2:::event_import_project(event_test_project(
      write_maf(row), "somatic_maf")), "MAF")
})

test_that("valid zero-event VCF runs complete and remain distinct from pre-run", {
  write_vcf <- function(records) {
    path <- tempfile(fileext = ".vcf")
    writeLines(c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample-1",
      records), path)
    path
  }
  for (records in list(character(), "chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/0")) {
    p <- event_test_project(write_vcf(records), "germline_vcf")
    expect_error(neo_events(p), "not available")
    p <- neo_run(p)
    expect_identical(Neoantimon2:::read_manifest(p$path)$status, "COMPLETED")
    event_cp <- Filter(function(x) identical(x$step_id, "event_import"),
      Neoantimon2:::read_manifest(p$path)$checkpoints)[[1L]]
    expect_identical(event_cp$status, "COMPLETED")
    expect_equal(nrow(neo_events(p)), 0L)
    expect_equal(nrow(neo_origins(p)), 0L)
    expect_equal(nrow(neo_results(p)), 0L)
    expect_equal(nrow(neo_multisample(p)), 0L)
    expect_silent(Neoantimon2:::validate_table(
      Neoantimon2:::read_runtime_snapshot(p)$events, "events"))
  }
  expect_error(neo_run(event_test_project(character(), character())), "requires event input")
  bad <- Neoantimon2:::initialize_run(event_test_project("bad_ref.vcf", "somatic_vcf"))
  expect_error(neo_run(bad), "REF mismatch")
  expect_identical(Neoantimon2:::read_manifest(bad$path)$status, "FAILED")
})

test_that("separate phase sources retain raw GT and PS without interpretation", {
  p <- event_test_project(c("phase_somatic.vcf", "phase_germline.vcf"),
                          c("somatic_vcf", "germline_vcf"))
  rows <- Neoantimon2:::event_import_project(p)
  expect_length(rows, 2L)
  tokens <- vapply(rows, function(x) x$event_details_json$source_records[[1L]]$raw_format$sample, "")
  expect_setequal(tokens, c("0|1:44", "1|0:44"))
  expect_true(all(vapply(rows, function(x) length(x$phase_links_json) == 0L, TRUE)))
  expect_true(all(vapply(rows, function(x) is.null(x$phase_set), TRUE)))
  expect_setequal(vapply(rows, `[[`, "", "event_role"), c("somatic", "germline"))
})

test_that("symbolic and BND records defer; star and invalid input fail explicitly", {
  rows <- Neoantimon2:::event_import_project(event_test_project("symbolic.vcf", "somatic_vcf"))
  expect_setequal(vapply(rows, `[[`, "", "event_type"), c("SV", "BND"))
  expect_identical(Filter(function(x) x$event_type == "SV", rows)[[1L]]$sv_type, "DEL")
  expect_true(all(vapply(rows, function(x) x$reconstruction_status == "UNSUPPORTED", TRUE)))
  cases <- c(star.vcf = "spanning-deletion", bad_ref.vcf = "REF mismatch",
             alias.vcf = "Contig absent", range.vcf = "outside",
             bad_allele.vcf = "Unsupported VCF ALT")
  for (name in names(cases)) {
    p <- event_test_project(name, "somatic_vcf")
    expect_error(Neoantimon2:::event_import_project(p), cases[[name]], info = name)
  }
  p <- event_test_project("bad_ref.vcf", "somatic_vcf")
  p <- Neoantimon2:::initialize_run(p)
  expect_error(neo_run(p), "REF mismatch")
  manifest <- Neoantimon2:::read_manifest(p$path)
  expect_identical(manifest$status, "FAILED")
  expect_identical(Filter(function(x) x$step_id == "event_import", manifest$checkpoints)[[1L]]$status, "FAILED")
})

test_that("MAF distinct nonreference tumor alleles split without guessing", {
  rows <- Neoantimon2:::event_import_project(event_test_project("ambiguous.maf", "somatic_maf"))
  expect_length(rows, 2L)
  expect_setequal(vapply(rows, `[[`, "", "alt"), c("C", "G"))
  expect_setequal(vapply(rows, function(x)
    x$event_details_json$source_records[[1L]]$selected_alt_index, 0L), c(1L, 2L))
  expect_error(Neoantimon2:::event_import_project(event_test_project("contradictory.maf", "somatic_maf")),
               "Malformed MAF tumor allele")
})

test_that("canonical T002 TSV validates source identity and rebinds with provenance", {
  source <- Neoantimon2:::event_import_project(event_test_project("duplicate.vcf", "somatic_vcf"))[[1L]]
  contract <- Neoantimon2:::schema_contract()
  fields <- vapply(contract$tables$events$fields, `[[`, "", "name")
  encode <- function(x) if (is.null(x)) "\\N" else if (is.list(x))
    as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null")) else as.character(x)
  path <- tempfile(fileext = ".tsv")
  writeLines(c(paste(fields, collapse = "\t"),
               paste(vapply(source[fields], encode, ""), collapse = "\t")), path)
  p <- event_test_project(path, "canonical_events")
  imported <- Neoantimon2:::event_import_project(p)[[1L]]
  expect_identical(imported$event_key, source$event_key)
  expect_identical(imported$event_id, source$event_id)
  expect_false(identical(imported$analysis_id, source$analysis_id))
  expect_identical(imported$event_details_json$source_records[[1L]]$original$source_event_id,
                   source$event_id)
  expect_identical(imported$event_details_json$source_records[[1L]]$original$source_analysis_context,
                   source$analysis_context)
  expect_identical(imported$event_details_json$source_records[[1L]]$original$source_patient_id,
                   source$patient_id)
  expect_identical(imported$event_details_json$source_records[[1L]]$original$source_canonical_event,
                   source$canonical_event)
  write_source <- function(row) {
    path <- tempfile(fileext = ".tsv")
    writeLines(c(paste(fields, collapse = "\t"),
      paste(vapply(row[fields], encode, ""), collapse = "\t")), path)
    path
  }
  wrong_patient <- source
  wrong_patient$patient_id <- "other-patient"
  wrong_patient$event_id <- Neoantimon2:::biological_id("event_id", list(
    context_id = wrong_patient$analysis_context, patient_id = wrong_patient$patient_id,
    sample_id = wrong_patient$sample_id, event_key = wrong_patient$event_key), contract)
  expect_error(Neoantimon2:::event_import_project(event_test_project(write_source(wrong_patient),
    "canonical_events")), "binding mismatch")
  other_allele <- source
  other_allele$alt <- if (source$alt == "C") "G" else "C"
  expect_error(Neoantimon2:::event_import_project(event_test_project(write_source(other_allele),
    "canonical_events")), "canonical_event disagrees")
  wrong_namespace <- source
  wrong_namespace$canonical_event$namespace <- "other-namespace"
  wrong_namespace$event_key <- Neoantimon2:::biological_id("event_key", list(
    genome_assembly = wrong_namespace$genome_assembly,
    canonical_event = wrong_namespace$canonical_event), contract)
  wrong_namespace$event_id <- Neoantimon2:::biological_id("event_id", list(
    context_id = wrong_namespace$analysis_context, patient_id = wrong_namespace$patient_id,
    sample_id = wrong_namespace$sample_id, event_key = wrong_namespace$event_key), contract)
  expect_error(Neoantimon2:::event_import_project(event_test_project(write_source(wrong_namespace),
    "canonical_events")), "canonical_event disagrees")
})

test_that("BND mate coordinates require exact primary reference bounds", {
  make_bnd <- function(mate) {
    path <- tempfile(fileext = ".vcf")
    writeLines(c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
      paste("chr1", 5, "bnd", "A", paste0("A]", mate, "]"), ".", "PASS", ".", sep = "\t")), path)
    path
  }
  for (mate in c("chr1:20", "chrM:5")) {
    rows <- Neoantimon2:::event_import_project(event_test_project(make_bnd(mate), "somatic_vcf"))
    expect_length(rows, 1L)
    expect_identical(rows[[1L]]$alt, paste0("A]", mate, "]"))
  }
  for (mate in c("absent:1", "1:5", "chr1:0", "chr1:121"))
    expect_error(Neoantimon2:::event_import_project(event_test_project(make_bnd(mate),
      "somatic_vcf")), "BND mate", info = mate)
})

test_that("germline GT selects asserted alleles and no GT leaves genotype unknown", {
  make_gt <- function(gt = NULL) {
    path <- tempfile(fileext = ".vcf")
    record <- "chr1\t1\tg\tA\tC,G\t.\tPASS\t."
    header <- "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    if (!is.null(gt)) { header <- paste0(header, "\tFORMAT\tsample-1"); record <- paste0(record, "\tGT\t", gt) }
    writeLines(c("##fileformat=VCFv4.2", header, record), path)
    path
  }
  for (case in list(list("0/1", 1L), list("1/2", 2L), list("0/0", 0L), list(NULL, 2L))) {
    p <- event_test_project(make_gt(case[[1L]]), "germline_vcf")
    rows <- if (identical(case[[2L]], 0L)) {
      reference <- Neoantimon2:::event_reference_open(Neoantimon2:::reference_bundle(p))
      value <- Neoantimon2:::event_parse_vcf(p, p$inputs[[1L]], reference)
      Neoantimon2:::event_reference_close(reference)
      value
    } else Neoantimon2:::event_import_project(p)
    expect_length(rows, case[[2L]])
    if (length(rows)) {
      src <- rows[[1L]]$event_details_json$source_records[[1L]]
      expect_identical(src$genotype_status, if (is.null(case[[1L]])) "UNKNOWN_NO_GT" else "RAW_GT_PRESENT")
      if (!is.null(case[[1L]])) expect_identical(src$raw_format$sample, case[[1L]])
    }
  }
  for (gt in c("1/3", "0/.", "1/", "a/b"))
    expect_error(Neoantimon2:::event_import_project(event_test_project(make_gt(gt),
      "germline_vcf")), "GT", info = gt)
})

test_that("legacy source adapters retain source mode without inventing coordinates", {
  write_source <- function(lines) { path <- tempfile(fileext = ".tsv"); writeLines(lines, path); path }
  fragment <- write_source(c("label\tsequence\tgroup_id\treference_nm_id\treference_gene_symbol\tsource_annotation",
    "frag-A\tMKTA\tgroup-1\tNM_001\tGENE\tRNA source"))
  entire <- write_source(c("label\tinput_nm_id\treading_frame\tgroup_id\tsource_annotation",
    "region-A\tNM_002\t1\tgroup-2\treference only"))
  for (case in list(list(fragment, "sequence_fragment", "fragment"),
                    list(entire, "sequence_entire_region", "entire-region"))) {
    row <- Neoantimon2:::event_import_project(event_test_project(case[[1L]], case[[2L]]))[[1L]]
    expect_identical(row$event_type, "SEQUENCE")
    expect_null(row$contig); expect_null(row$pos); expect_null(row$ref); expect_null(row$alt)
    expect_identical(row$event_details_json$source_records[[1L]]$original$source_mode, case[[3L]])
    expect_silent(Neoantimon2:::validate_table(list(row), "events"))
  }
  annovar <- write_source(c("Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene",
    "chr1\t1\t1\tA\tC\texonic"))
  row <- Neoantimon2:::event_import_project(event_test_project(annovar, "somatic_annovar"))[[1L]]
  expect_identical(row$alt, "C")
  expect_identical(row$event_details_json$source_records[[1L]]$original$Func.refGene, "exonic")
  vep <- write_source(c("Uploaded_variation\tLocation\tAllele\tConsequence",
    "variant-A\tchr1:1\tC\tmissense_variant"))
  row <- Neoantimon2:::event_import_project(event_test_project(vep, "somatic_vep"))[[1L]]
  expect_identical(row$event_details_json$source_records[[1L]]$reference_allele_origin,
                   "selected_T005_reference_not_source_column")
  vep_deletion <- write_source(c("Uploaded_variation\tLocation\tAllele\tConsequence",
    "variant-del\tchr1:2-3\t-\tframeshift_variant"))
  deletion <- Neoantimon2:::event_import_project(event_test_project(vep_deletion, "somatic_vep"))[[1L]]
  expect_identical(deletion$event_type, "DEL")
  expect_identical(deletion$event_details_json$source_records[[1L]]$original$Location, "chr1:2-3")
  vep_insertion <- write_source(c("Uploaded_variation\tLocation\tAllele\tConsequence",
    "variant-ins\tchr1:3-2\tA\tframeshift_variant"))
  insertion <- Neoantimon2:::event_import_project(event_test_project(vep_insertion, "somatic_vep"))[[1L]]
  expect_identical(insertion$event_type, "INS")
  expect_identical(insertion$event_details_json$source_records[[1L]]$original$Allele, "A")
  p <- neo_run(event_test_project(fragment, "sequence_fragment"))
  expect_equal(nrow(neo_events(p)), 1L)
  expect_equal(nrow(neo_origins(p)), 0L)
  snps <- write_source(c("CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
    "chr1\t1\trs1\tA\tC\t.\tPASS\t.\tGT:AD"))
  snp <- Neoantimon2:::event_import_project(event_test_project(snps, "germline_snp_table"))[[1L]]
  expect_identical(snp$event_role, "germline")
  expect_identical(snp$event_details_json$source_records[[1L]]$genotype_status,
                   "UNKNOWN_NO_SAMPLE_GT")
  bnd <- write_source(c("Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tmateID",
    "chr1\t5\t5\tA\tA]chrM:5]\texonic\tGENE\tmate-1"))
  bnd_row <- Neoantimon2:::event_import_project(event_test_project(bnd, "somatic_bnd_table"))[[1L]]
  expect_identical(bnd_row$event_type, "BND")
  expect_identical(bnd_row$alt, "A]chrM:5]")
  expect_identical(bnd_row$event_details_json$source_records[[1L]]$original$mateID, "mate-1")
  bad_bnd <- write_source(c("Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tmateID",
    "chr1\t5\t5\tA\tA]1:5]\texonic\tGENE\tmate-1"))
  expect_error(Neoantimon2:::event_import_project(event_test_project(bad_bnd, "somatic_bnd_table")),
               "BND mate contig absent")
})

test_that("INFO order, source row order, and caller order do not alter biological IDs", {
  original <- readLines(event_fixture("somatic.vcf"))
  reordered <- original
  record <- grep("^chr1\\t1\\t", reordered)[1L]
  reordered[[record]] <- sub("AF=0.2;DP=20", "DP=20;AF=0.2", reordered[[record]], fixed = TRUE)
  data_rows <- grep("^chr1\\t", reordered)
  reordered[data_rows] <- rev(reordered[data_rows])
  path <- tempfile(fileext = ".vcf"); writeLines(reordered, path)
  a <- Neoantimon2:::event_import_project(event_test_project("somatic.vcf", "somatic_vcf"))
  b <- Neoantimon2:::event_import_project(event_test_project(path, "somatic_vcf"))
  expect_identical(sort(vapply(a, `[[`, "", "event_key")),
                   sort(vapply(b, `[[`, "", "event_key")))
  expect_identical(sort(vapply(a, `[[`, "", "event_type")),
                   sort(vapply(b, `[[`, "", "event_type")))
})

test_that("source or selected reference changes cannot reuse event import", {
  p <- event_test_project("duplicate.vcf", "somatic_vcf")
  source <- p$inputs[[1L]]$resolved_path
  changed <- tempfile(fileext = ".vcf")
  writeLines(readLines(source), changed)
  p2 <- event_test_project(changed, "somatic_vcf")
  write("#changed", file = changed, append = TRUE)
  expect_error(Neoantimon2:::event_import_project(p2), "changed since project binding")
  p3 <- event_test_project("duplicate.vcf", "somatic_vcf")
  genome <- Filter(function(x) x$logical_name == "reference.genome_fasta", p3$resources)[[1L]]$resolved_path
  write("A", file = genome, append = TRUE)
  expect_error(Neoantimon2:::event_import_project(p3), "Reference changed")
})

test_that("completed event checkpoint corruption is detected on immutable run", {
  p <- neo_run(event_test_project("duplicate.vcf", "somatic_vcf"))
  path <- Neoantimon2:::event_stage_path(p)
  write("corrupt", path, append = TRUE)
  expect_error(neo_events(p), "Declared checksum mismatch")
  expect_error(neo_run(p), "Declared checksum mismatch")
  expect_identical(Neoantimon2:::read_manifest(p$path)$status, "COMPLETED")
})

test_that("VCF declared reference conflicts fail and source strand is not applied", {
  source <- readLines(event_fixture("somatic.vcf"))
  bad_length <- tempfile(fileext = ".vcf")
  writeLines(sub("length=120", "length=119", source, fixed = TRUE), bad_length)
  expect_error(Neoantimon2:::event_import_project(event_test_project(bad_length, "somatic_vcf")),
               "contig length mismatch")
  bad_build <- tempfile(fileext = ".vcf")
  writeLines(sub("##reference=GRCh38", "##reference=GRCh37", source, fixed = TRUE), bad_build)
  expect_error(Neoantimon2:::event_import_project(event_test_project(bad_build, "somatic_vcf")),
               "contradicts selected reference")
  strand <- tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
               "chr1\t1\t.\tA\tC\t.\tPASS\tSTRAND=-"), strand)
  row <- Neoantimon2:::event_import_project(event_test_project(strand, "somatic_vcf"))[[1L]]
  expect_identical(row$ref, "A")
  expect_identical(row$alt, "C")
  expect_identical(row$event_details_json$source_records[[1L]]$raw_info, "STRAND=-")
})

test_that("deferred structural assertions distinguish END and ignore row order", {
  lines <- readLines(event_fixture("symbolic.vcf"))
  changed <- tempfile(fileext = ".vcf")
  writeLines(sub("END=10", "END=11", lines, fixed = TRUE), changed)
  a <- Neoantimon2:::event_import_project(event_test_project("symbolic.vcf", "somatic_vcf", "symbolic"))
  b <- Neoantimon2:::event_import_project(event_test_project(changed, "somatic_vcf", "symbolic"))
  del <- function(rows) Filter(function(x) x$event_type == "SV", rows)[[1L]]
  expect_false(identical(del(a)$event_key, del(b)$event_key))
  expect_identical(del(a)$end, 10L)
  reordered <- tempfile(fileext = ".vcf")
  idx <- grep("^chr1\\t", lines)
  lines[idx] <- rev(lines[idx]); writeLines(lines, reordered)
  c <- Neoantimon2:::event_import_project(event_test_project(reordered, "somatic_vcf", "symbolic"))
  expect_identical(sort(vapply(a, `[[`, "", "event_key")),
                   sort(vapply(c, `[[`, "", "event_key")))
})
