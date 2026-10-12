longitudinal_fixture <- function(file) normalizePath(testthat::test_path("..", "fixtures", "longitudinal", file))
longitudinal_settings_fixture <- function() jsonlite::read_json(longitudinal_fixture("contexts.json"))
longitudinal_project <- function(settings, extra = list(), fixture_kind = "real") {
  path <- tempfile(fileext = ".json"); jsonlite::write_json(settings, path, auto_unbox = TRUE, null = "null")
  samples <- lapply(settings$contexts, function(s) list(sample_id = s$sample_id, patient_id = s$subject_id))
  cfg <- neo_config(env = character(), checksum = "sha256", paths = list(run_root = tempfile("longitudinal-runs-")))
  neo_project(samples = samples, config = cfg, inputs = c(list(list(input_id = "longitudinal-settings",
    input_role = "longitudinal_settings", sample_id = samples[[1L]]$sample_id, path = path, version = "1")), extra), fixture_kind = fixture_kind)
}
longitudinal_assay <- function(sample, key, domain = "candidate", state = "evaluated_negative") list(
  sample_id = sample, feature_key = key, domain = domain, state = state, adequate = TRUE,
  assay = "synthetic source-aware coverage and generation audit", evidence_ref = "synthetic:evaluability",
  policy_version = "fixture-opportunity-v1")
longitudinal_snv <- function(sample = "S1", pos = 4L, alt = "G", ref = "A") list(
  event_id = paste0("event:", sample, ":", pos, alt), sample_id = sample, genome_assembly = "GRCh38.p14",
  contig = "chrW", pos = pos, ref = ref, alt = alt, canonical_event = list(namespace = "t006-primary-small-variant-v1",
    identity = list(event_role = "somatic", contig = "chrW", pos = as.integer(pos), ref = ref, alt = alt)))
longitudinal_cells <- function(states, ordered = TRUE, regions = TRUE) {
  contexts <- lapply(seq_along(states), function(i) {
    s <- list(sample_id = paste0("S", i), subject_id = "patient-1")
    if (ordered) s$ordinal_time <- as.numeric(i)
    if (regions) s$region_id <- paste0("R", i)
    s
  })
  cells <- lapply(seq_along(states), function(i) list(sample_id = paste0("S", i), state = states[[i]], observations = list()))
  Neoantimon2:::longitudinal_summary(cells, contexts, list())
}

# A ligand/presenter has origin X in S1, origin Y in S2, and no candidate in
# S3. Both events and their imported clonal evidence exist in every sample.
longitudinal_origin_fixture <- function(mapping = "mapping-identity.tsv") {
  settings <- longitudinal_settings_fixture()
  settings$providers <- list(list(input_id = "pc", provider = "PyClone-VI",
    provider_version = "synthetic-0.1.6", run_id = "external-origin-fixture",
    subject_id = "patient-1", mapping_input_id = "map", options = list(seed = 1)))
  extra <- list(list(input_id = "pc", input_role = "pyclone_vi", sample_id = "S1",
    path = longitudinal_fixture("candidate-origin-clonal.tsv"), version = "1"),
    list(input_id = "map", input_role = "clonal_mapping", sample_id = "S1",
      path = longitudinal_fixture(mapping), version = "mapping-v1"))
  project <- longitudinal_project(settings, extra)
  key <- Neoantimon2:::biological_id("candidate_key", list(species = "NCBITaxon:9606",
    hla_class = "I", ligand_sequence = "GILGFVFTL", mhc_key = "HLA-A*02:01"),
    Neoantimon2:::peptide_contract())
  graphs <- lapply(seq_len(3L), function(i) {
    sample <- paste0("S", i); id <- paste0("candidate:", sample)
    events <- lapply(seq_len(2L), function(j) {
      e <- longitudinal_snv(sample, pos = c(4L, 13L)[[j]])
      e$event_id <- paste0(c("X", "Y")[[j]], "_", sample); e
    })
    row <- data.frame(sample_id = sample, candidate_id = id, candidate_key = key,
      ligand_sequence = "GILGFVFTL", mhc_key = "HLA-A*02:01", rna_support = FALSE)
    candidate <- list(candidate_id = id, candidate_key = key, ligand_sequence = "GILGFVFTL",
      mhc_key = "HLA-A*02:01", hla_class = "I", origin_ids = list(paste0("origin:", sample)))
    origins <- if (i < 3L) list(list(sample_id = sample,
      source_event_ids = list(events[[i]]$event_id))) else list()
    list(sample_id = sample, subject_id = "patient-1", events = events,
      candidates = if (i < 3L) list(candidate) else list(), results = row,
      project = list(sample_id = sample, hla = list("patient-1" = "HLA-A*02:01")),
      explanation = list(result = as.list(row[1L, ]), origins = origins))
  })
  testthat::with_mocked_bindings(
    Neoantimon2:::longitudinal_construct(project, settings, graphs),
    longitudinal_base_explain = function(project, id) list(candidates = list(
      Filter(function(g) g$sample_id == project$sample_id, graphs)[[1L]]$explanation)),
    .package = "Neoantimon2")
}

longitudinal_without_physical <- function(x) {
  if (!is.list(x)) return(x)
  x[c("artifact", "mapping_artifact")[c("artifact", "mapping_artifact") %in% names(x)]] <- NULL
  lapply(x, longitudinal_without_physical)
}
longitudinal_real_sample <- function(id, later = FALSE, extra_variant = FALSE) {
  base <- reconstruction_test_project(); cfg <- base$config
  cfg$peptides <- peptide_test_policy(8L, 9L)
  cfg$predictors <- list(list(provider = "netmhcpan", task = "mhc1_presentation",
    executable = recognition_mock_presentation_executable("netmhcpan"), options = list()),
    list(provider = "netmhciipan", task = "mhc2_presentation",
      executable = recognition_mock_presentation_executable("netmhciipan"), options = list()))
  path <- tempfile(fileext = ".vcf")
  records <- c("chrW\t4\t.\tA\tG\t.\tPASS\t.\tGT\t0/1",
    paste0("chrW\t13\t.\tA\t", if (later) "C" else "G", "\t.\tPASS\t.\tGT\t0/1"))
  if (extra_variant) records <- c(records, "chrF\t4\t.\tA\tC\t.\tPASS\t.\tGT\t0/1")
  writeLines(c("##fileformat=VCFv4.2", "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
    paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", id), records), path)
  samples <- list(list(sample_id = id, patient_id = "patient-1"))
  inputs <- list(list(input_id = paste0(id, "-vcf"), input_role = "somatic_vcf", sample_id = id, path = path, version = "synthetic-1"))
  hla <- list("patient-1" = c("HLA-A*02:01", "HLA-A*24:02", "HLA-DRA*01:01~HLA-DRB1*04:01"))
  pre <- neo_project(samples = samples, config = cfg, hla = hla, inputs = inputs)
  events <- Neoantimon2:::event_import_project(pre)
  dna <- tempfile(fileext = ".tsv"); cn <- tempfile(fileext = ".tsv"); rna <- tempfile(fileext = ".tsv")
  writeLines(c("event_id\ttumor_ref_reads\ttumor_alt_reads\ttumor_depth\tnormal_ref_reads\tnormal_alt_reads\tnormal_depth",
    vapply(events, function(e) paste(e$event_id, 20, 5, 25, 20, 0, 20, sep = "\t"), "")), dna)
  writeLines(c("event_id\trna_ref_reads\trna_alt_reads\trna_depth", vapply(events, function(e) paste(e$event_id, 9, 1, 10, sep = "\t"), "")), rna)
  writeLines(c("chrom\tloc.start\tloc.end\ttcn.em\tlcn.em\tnormal_cn\tpurity", "W\t1\t42\t2\t1\t2\t0.8", "F\t1\t51\t2\t1\t2\t0.8"), cn)
  cfg$somatic_evidence <- list(sources = list(list(kind = "dna_reads", sample_id = id, path = dna),
    list(kind = "rna_reads", sample_id = id, path = rna), list(kind = "facets_cn", sample_id = id, path = cn, reference_build = "GRCh38.p14")))
  loh <- tempfile(fileext = ".tsv")
  tab <- readLines(normalizePath(testthat::test_path("..", "fixtures", "hla_loh", "lohhla-standard.HLAlossPrediction_CI.xls")))
  if (!later) tab[[2L]] <- sub("0.0099", "0.2", sub("0.4999", "1.0", tab[[2L]], fixed = TRUE), fixed = TRUE)
  writeLines(tab, loh)
  cfg$hla_loh <- list(sources = list(list(sample_id = id, path = loh, provider_version = "synthetic-1")))
  inputs <- c(inputs, list(list(input_id = "apm", input_role = "apm_evidence", sample_id = id,
    path = longitudinal_fixture(if (later) "apm-after.tsv" else "apm-before.tsv"), version = "1")))
  if (!later) inputs <- c(inputs, list(recognition_fixture_input("bulk.tsv", sample = id)))
  p <- neo_run(neo_project(samples = samples, config = cfg, hla = hla, inputs = inputs))
  if (id == "S1") {
    ms <- presentation_write_psms(unique(neo_results(p)$ligand_sequence)); ms$sample_id <- id
    p <- neo_run(recognition_clone_with(p, list(ms)))
  }
  p
}
longitudinal_e2e <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    runs <- list(longitudinal_real_sample("S1", extra_variant = TRUE), longitudinal_real_sample("S2"), longitudinal_real_sample("S3", later = TRUE))
    before <- lapply(runs, neo_results)
    keys <- lapply(runs, function(p) vapply(Neoantimon2:::recognition_candidates(p, Neoantimon2:::read_runtime_snapshot(p)), `[[`, "", "candidate_key"))
    s <- longitudinal_settings_fixture(); s$sample_runs <- lapply(seq_along(runs), function(i) list(sample_id = paste0("S", i), path = runs[[i]]$path))
    for (i in seq_along(s$contexts)) s$contexts[[i]]$hla_genotype <- unname(unlist(runs[[i]]$hla))
    groups <- list(A = Reduce(intersect, keys), B = setdiff(keys[[2L]], keys[[3L]]),
      C = setdiff(keys[[3L]], union(keys[[1L]], keys[[2L]])), D = setdiff(keys[[1L]], union(keys[[2L]], keys[[3L]])))
    if (any(lengths(groups) == 0L)) {
      saveRDS(list(runs=runs,keys=keys,groups=groups),tempfile("neo-longitudinal-distribution-debug-", fileext=".rds"))
      stop("Real E2E requires all four biological distributions: ",paste(lengths(groups),collapse=","))
    }
    s$evaluability <- c(lapply(groups$B, function(k) longitudinal_assay("S3", k)),
      unlist(lapply(c("S1", "S2"), function(id) lapply(groups$C, function(k) longitudinal_assay(id, k))), recursive = FALSE))
    event <- lapply(runs, function(p) Neoantimon2:::read_runtime_snapshot(p)$events)
    for (i in c(1L, 2L)) for (e in event[[3L]]) if (!Neoantimon2:::longitudinal_event_key(e) %in% vapply(event[[i]], Neoantimon2:::longitudinal_event_key, ""))
      s$evaluability <- c(s$evaluability, list(longitudinal_assay(paste0("S", i), Neoantimon2:::longitudinal_event_key(e), "event")))
    for (e in event[[2L]]) if (!Neoantimon2:::longitudinal_event_key(e) %in% vapply(event[[3L]], Neoantimon2:::longitudinal_event_key, ""))
      s$evaluability <- c(s$evaluability, list(longitudinal_assay("S3", Neoantimon2:::longitudinal_event_key(e), "event")))
    path <- tempfile(fileext = ".tsv")
    writeLines(c("mutation_id\tsample_id\tcluster_id\tcellular_prevalence\tcellular_prevalence_std\tcluster_assignment_prob",
      vapply(seq_along(event), function(i) paste(event[[i]][[1L]]$event_id, paste0("S", i), 1, .72 - i * .1, .08, .65, sep = "\t"), "")), path)
    s$providers <- list(list(input_id = "pyclone", provider = "PyClone-VI", provider_version = "synthetic-0.1.6", run_id = "external-fit-1", subject_id = "patient-1", options = list(seed = 4), reference_build = "GRCh38.p14"))
    extra <- list(list(input_id = "pyclone", input_role = "pyclone_vi", sample_id = "S1", path = path, version = "synthetic-1"))
    p <- neo_run(longitudinal_project(s, extra))
    cache <<- list(project = p, overlay = Neoantimon2:::longitudinal_read(p), runs = runs, before = before, groups = groups, settings = s, extra = extra)
    cache
  }
})
