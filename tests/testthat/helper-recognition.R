recognition_fixture_path <- function(file) normalizePath(testthat::test_path("..", "fixtures", "tcr_recognition", file))
recognition_functionality_cases <- function() jsonlite::read_json(recognition_fixture_path("functionality.json"))
recognition_functionality_input <- function(flags, locus = "TRB", sample = "sample-1", id = "functionality") {
  rows <- Neoantimon2:::source_table(recognition_fixture_path("unique-paired.tsv"),
    c("sequence_id", "locus", "junction_aa"), "AIRR")
  for (i in seq_along(rows)) if (rows[[i]]$locus == locus) {
    rows[[i]][c("productive", "stop_codon", "vj_in_frame")] <- NULL
    for (key in names(flags)) rows[[i]][[key]] <- flags[[key]]
  }
  path <- tempfile(fileext = ".json")
  x <- jsonlite::read_json(recognition_fixture_path("unique-cell.json"))
  x$Rearrangement <- rows
  jsonlite::write_json(x, path, auto_unbox = TRUE, null = "null")
  list(input_id = id, input_role = "tcr_airr", sample_id = sample, path = path, version = "synthetic-2.0")
}
recognition_functionality_repertoire <- function(flags, locus = "TRB") Neoantimon2:::tcr_import(
  recognition_fixture_project(list(recognition_functionality_input(flags, locus))))
recognition_fixture_input <- function(file, role = "tcr_rearrangement", id = file, sample = "sample-1")
  list(input_id = id, input_role = role, sample_id = sample, path = recognition_fixture_path(file), version = "synthetic-2.0")
recognition_fixture_project <- function(inputs = list(recognition_fixture_input("bulk.tsv"))) {
  cfg <- neo_config(env = character(), paths = list(run_root = tempfile("recognition-runs-")), checksum = "sha256")
  neo_project(samples = data.frame(sample_id = c("sample-1", "normal"), patient_id = "patient-1"),
    inputs = inputs, config = cfg, hla = list("patient-1" = "HLA-A*02:01"), fixture_kind = "synthetic")
}
recognition_fixture_repertoire <- function(paired = FALSE) Neoantimon2:::tcr_import(recognition_fixture_project(
  if (paired) list(recognition_fixture_input("paired.tsv"), recognition_fixture_input("airr.json", "tcr_airr")) else
    list(recognition_fixture_input("bulk.tsv"))))
recognition_test_candidate <- function(peptide = "GILGFVFTL", hla = "HLA-A*02:01", sample = "sample-1") {
  p <- Neoantimon2:::peptide_presenters(list(hla))[[1L]]
  list(candidate_id = paste0("candidate:", digest::digest(list(peptide, p$mhc_key, sample))), ligand_sequence = peptide,
    peptide_sequence = peptide, sample_id = sample, hla_class = p$hla_class, mhc_key = p$mhc_key,
    hla_allele = p$hla_allele, hla_alpha = p$hla_alpha, hla_beta = p$hla_beta,
    final_grade = "NOT_GRADED", priority_origin_id = NULL)
}
recognition_test_cap <- function(provider = "fixture", model = "beta_only", ...) {
  x <- c(list(provider = provider, model = model, fixture = TRUE), list(...))
  if (provider == "prime") x$model <- "prime-2.1"
  Neoantimon2:::recognition_capability(x)
}
recognition_test_settings <- function(...) {
  x <- list(candidate_ids = NULL, profile_ids = NULL, max_queries = 10000, max_pairing_hypotheses = 1000, batch_size = 32)
  for (n in names(list(...))) x[[n]] <- list(...)[[n]]
  x
}
recognition_paired_profile <- function(r = recognition_fixture_repertoire(TRUE))
  Filter(function(p) p$completeness == "paired_cdr3_profile" && p$features$alpha_cdr3_aa == "CAVRDSNYQLIW" && is.null(p$features$t_cell_type), r$profiles)[[1L]]
recognition_test_predict <- function(provider, file, model = "beta_only", use_mhc = TRUE) {
  cap <- recognition_test_cap(provider, model, use_mhc = use_mhc)
  r <- recognition_fixture_repertoire(TRUE)
  prepared <- Neoantimon2:::recognition_prepare(list(recognition_test_candidate()), r, cap,
    recognition_test_settings(profile_ids = list(recognition_paired_profile(r)$receptor_feature_profile_id)))
  input <- Neoantimon2:::recognition_inputs(recognition_fixture_project(list(
    recognition_fixture_input(file, "recognition_output"))))[[1L]]
  Neoantimon2:::recognition_predict(prepared, cap, input)
}
recognition_temp_settings <- function(providers = list(list(provider = "fixture", fixture = TRUE, model = "paired")), ...) {
  path <- tempfile(fileext = ".json")
  jsonlite::write_json(c(list(version = "recognition-run-settings-v1", providers = providers), list(...)), path,
    auto_unbox = TRUE, pretty = TRUE)
  list(input_id = "recognition-settings", input_role = "recognition_settings", sample_id = "sample-1", path = path, version = "1")
}
recognition_mock_presentation_executable <- function(provider) {
  path <- tempfile("recognition-presentation-mock-")
  version <- if (provider == "netmhcpan") "NetMHCpan version 4.2" else "NetMHCIIpan version 4.3"
  writeLines(c("#!/bin/sh", paste0("if [ \"$1\" = \"-h\" ]; then echo '# ", version, "'; exit 0; fi"),
    "input=''", "allele=''", "while [ $# -gt 0 ]; do",
    "  case \"$1\" in -p|-f) input=$2; shift 2 ;; -a) allele=$2; shift 2 ;; *) shift ;; esac", "done",
    "if [ -z \"$input\" ]; then exit 2; fi",
    if (provider == "netmhcpan") "echo 'Pos HLA Peptide Score_EL %Rank_EL'" else "echo 'Pos MHC Peptide Score_EL %Rank_EL'",
    "n=0", "while IFS= read -r peptide; do", "  n=$((n+1))",
    "  case $((n%4)) in 1) rank=0.5 ;; 2) rank=0.5 ;; 3) rank=2.0 ;; 0) rank=6.0 ;; esac",
    "  printf '%s %s %s 0.5 %s\\n' \"$n\" \"$allele\" \"$peptide\" \"$rank\"",
    "done < \"$input\""), path)
  Sys.chmod(path, "0755"); path
}
recognition_gold_project <- function() {
  base <- reconstruction_test_project()
  cfg <- base$config; cfg$peptides <- peptide_test_policy(8L, 9L)
  cfg$predictors <- list(list(provider = "netmhcpan", task = "mhc1_presentation",
    executable = recognition_mock_presentation_executable("netmhcpan"), options = list()),
    list(provider = "netmhciipan", task = "mhc2_presentation",
      executable = recognition_mock_presentation_executable("netmhciipan"), options = list()))
  path <- tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.2", "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample-1",
    "chrW\t4\t.\tA\tG\t.\tPASS\t.\tGT\t0/1", "chrW\t13\t.\tA\tG\t.\tPASS\t.\tGT\t0/1"), path)
  inputs <- list(list(input_id = "recognition-gold-vcf", input_role = "somatic_vcf", sample_id = "sample-1", path = path))
  hla <- list("patient-1" = c("HLA-A*02:01", "HLA-DRA*01:01~HLA-DRB1*04:01"))
  pre <- neo_project(samples = base$samples, config = cfg, hla = hla, inputs = inputs)
  events <- Neoantimon2:::event_import_project(pre)
  dna <- tempfile(fileext = ".tsv"); rna <- tempfile(fileext = ".tsv"); cn <- tempfile(fileext = ".tsv")
  writeLines(c("event_id\ttumor_ref_reads\ttumor_alt_reads\ttumor_depth\tnormal_ref_reads\tnormal_alt_reads\tnormal_depth",
    vapply(events, function(e) paste(e$event_id, 20, 5, 25, 19, 1, 20, sep = "\t"), "")), dna)
  writeLines(c("event_id\trna_ref_reads\trna_alt_reads\trna_depth",
    paste(events[[1L]]$event_id, 9, 1, 10, sep = "\t"), paste(events[[2L]]$event_id, 10, 0, 10, sep = "\t")), rna)
  writeLines(c("chrom\tloc.start\tloc.end\ttcn.em\tlcn.em\tnormal_cn\tpurity", "W\t1\t42\t2\t0\t2\t0.8"), cn)
  cfg$somatic_evidence <- list(sources = list(list(kind = "dna_reads", sample_id = "sample-1", path = dna),
    list(kind = "rna_reads", sample_id = "sample-1", path = rna), list(kind = "facets_cn", sample_id = "sample-1",
      path = cn, reference_build = cfg$scientific$genome_assembly)))
  neo_project(samples = base$samples, config = cfg, hla = hla, inputs = inputs)
}
recognition_clone_inputs <- function(project) lapply(project$inputs, function(i)
  list(input_id = i$input_id, input_role = i$input_role, sample_id = i$sample_id, path = i$resolved_path, version = i$version))
recognition_clone_with <- function(project, extra) neo_project(samples = project$samples, config = project$config,
  hla = lapply(project$hla, unlist), inputs = c(recognition_clone_inputs(project), extra), fixture_kind = project$fixture_kind)
recognition_invariance_fields <- function(rows) {
  fields <- intersect(c("candidate_id", "final_grade", "grading_status", "technical_eligibility", "rank_within_sample_class",
    "rank_status", "priority_origin_id", "primary_rank", "primary_predictor", "primary_metric", "primary_prediction_status",
    "prediction_status", "grade_reason", "policy_version"), names(rows))
  rows[order(rows$candidate_id), fields, drop = FALSE]
}
recognition_e2e_pair <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    before <- neo_run(recognition_gold_project()); rows <- neo_results(before)
    first <- rows[which(rows$hla_class == "I")[[1L]], ]
    r <- list(reactivity_id = "e2e-measurement", cell_id = "observed", repertoire_id = "cells", ligand_type = "MHC:peptide",
      antigen_type = "peptide", antigen = "synthetic:antigen", peptide_sequence_aa = first$ligand_sequence,
      mhc_allele_1 = first$mhc_key, mhc_class = first$hla_class, reactivity_method = "ELISpot",
      reactivity_readout = "spots", reactivity_value = 120, reactivity_unit = "SFU", qualitative_annotation = "positive",
      threshold = 10, threshold_reference = "synthetic:study#cutoff", reactivity_refs = list("synthetic:study"))
    path <- tempfile(fileext = ".json"); jsonlite::write_json(list(Info = list(version = "2.0"), Reactivity = list(r)), path, auto_unbox = TRUE)
    # User-supplied PRIME-format output for every supported core Class I peptide.
    prime_path <- tempfile(fileext = ".tsv")
    peptides <- sort(unique(rows$ligand_sequence[rows$hla_class == "I"]))
    writeLines(c("Peptide\t%Rank_A0201\tScore_A0201\t%RankBinding_A0201",
      paste(peptides, "0.25", "0.8", "0.5", sep = "\t")), prime_path)
    extra <- list(recognition_fixture_input("unique-paired.tsv"), recognition_fixture_input("unique-cell.json", "tcr_airr"),
      list(input_id = "e2e-reactivity", input_role = "tcr_reactivity", sample_id = "sample-1", path = path, version = "2.0"),
      list(input_id = "e2e-prime", input_role = "recognition_output", sample_id = "sample-1", path = prime_path, version = "synthetic-2.1"),
      recognition_temp_settings(list(list(provider = "fixture", fixture = TRUE, model = "paired"),
        list(provider = "prime", provider_version = "2.1", model = "prime-2.1", output_input_id = "e2e-prime", fixture = TRUE))))
    after <- neo_run(recognition_clone_with(before, extra))
    cache <<- list(before = before, after = after, measured_candidate = first$candidate_id)
    cache
  }
})
