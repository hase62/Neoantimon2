presentation_fixture_path <- function(file) normalizePath(testthat::test_path("..", "fixtures", "presentation_evidence", file))
presentation_fixture_input <- function(file, role = "processing_output", id = file, sample = "sample-1")
  list(input_id = id, input_role = role, sample_id = sample, path = presentation_fixture_path(file), version = "synthetic-1")
presentation_bound_input <- function(file, role = "processing_output", sample = "sample-1") {
  p <- recognition_fixture_project(list(presentation_fixture_input(file, role, sample = sample)))
  p$inputs[[1L]]
}
presentation_candidate <- function(peptide = "GILGFVFTL", hla = "HLA-A*02:01", sample = "sample-1") {
  c <- recognition_test_candidate(peptide, hla, sample)
  c$peptide_key <- Neoantimon2:::biological_id("peptide_key", list(ligand_sequence = peptide), Neoantimon2:::peptide_contract())
  c$origin_ids <- list("origin:1", "origin:2")
  c
}
presentation_origin <- function(peptide = "GILGFVFTL", n = "AAAAA", c = "CCCCC", id = "origin:1")
  list(origin_id = id, sequence_context_id = paste0("protein:", id), source_protein_id = paste0("source:", id),
    peptide_key = presentation_candidate(peptide)$peptide_key, ligand_sequence = peptide,
    peptide_start_aa = nchar(n) + 1L, peptide_end_aa = nchar(n) + nchar(peptide),
    mutant_sequence_ref = paste0("authoritative#", id))
presentation_query <- function(provider = "fixture", peptide = "GILGFVFTL", n = "AAAAA", c = "CCCCC",
                               id = "origin:1", hla = "HLA-A*02:01", ...) {
  cap <- Neoantimon2:::processing_capability(c(list(provider = provider), list(...)))
  Neoantimon2:::processing_query(presentation_candidate(peptide, hla), presentation_origin(peptide, n, c, id),
    paste0(n, peptide, c), cap)
}
presentation_experiment <- function(ip = "class_i_pan", presenter = NULL, basis = NULL, class = "I")
  list(input_id = "psms.tsv", ms_run = "run-1", ip_type = ip, hla_class = class, presenter = presenter,
    assignment_basis = basis, assignment_reference = if (is.null(basis)) NULL else "synthetic:experiment",
    context = list(sample_id = "sample-1", region = NULL, timepoint = NULL, tissue = NULL, cell_context = NULL))
presentation_ms <- function(experiment = presentation_experiment(), file = "psms.tsv", role = "immunopeptidomics_psm") {
  input <- presentation_bound_input(file, role)
  psms <- if (role == "immunopeptidomics_psm") Neoantimon2:::ms_generic_import(input, experiment) else
    Neoantimon2:::ms_mzidentml_import(input, experiment)
  Neoantimon2:::ms_observations(psms)
}
presentation_temp_settings <- function(providers = list(list(provider = "fixture")), experiments = list(), max_queries = 10000L) {
  path <- tempfile(fileext = ".json")
  jsonlite::write_json(list(version = "presentation-run-settings-v1", providers = providers,
    experiments = experiments, max_queries = max_queries), path, auto_unbox = TRUE, pretty = TRUE)
  list(input_id = "presentation-settings", input_role = "presentation_settings", sample_id = "sample-1", path = path, version = "1")
}
presentation_write_psms <- function(peptides, id = "e2e-ms") {
  path <- tempfile(fileext = ".tsv")
  writeLines(c("psm_id\tspectrum_id\tpeptide_sequence\tdecoy\tpassThreshold\tq_value",
    vapply(seq_along(peptides), function(i) paste(paste0("psm", i), paste0("scan=", i), peptides[[i]], "FALSE", "TRUE", "0.001", sep = "\t"), "")), path)
  list(input_id = id, input_role = "immunopeptidomics_psm", sample_id = "sample-1", path = path, version = "synthetic-1")
}
presentation_e2e_pair <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    before <- neo_run(recognition_gold_project()); rows <- neo_results(before)
    class_i <- rows[rows$hla_class == "I", ]
    stab <- tempfile(fileext = ".txt")
    writeLines(c("pos HLA peptide Identity Pred Thalf(h) %Rank_Stab BindLevel",
      vapply(seq_len(nrow(class_i)), function(i) paste(0, class_i$mhc_key[[i]], class_i$ligand_sequence[[i]], "PEPLIST", "0.95 12 0.01"), "")), stab)
    extra <- list(presentation_write_psms(unique(rows$ligand_sequence)),
      presentation_fixture_input("apm.tsv", "apm_evidence"), presentation_fixture_input("hla-expression.tsv", "hla_expression"),
      list(input_id = "e2e-stability", input_role = "processing_output", sample_id = "sample-1", path = stab, version = "1.0"),
      presentation_temp_settings(list(list(provider = "fixture"), list(provider = "netmhcstabpan", output_input_id = "e2e-stability")),
        list(list(input_id = "e2e-ms", ip_type = "pan_hla", ms_run = "e2e-run"))))
    after <- neo_run(recognition_clone_with(before, extra))
    cache <<- list(before = before, after = after)
    cache
  }
})
presentation_source_pair <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    inputs <- list(long_read_fixture_input("novel-exon.gtf"),
      long_read_fixture_input("settings-denovo.json", "long_read_settings", "long-read-settings"),
      list(input_id="fusion",input_role="star_fusion",sample_id="sample-1",
        path=normalizePath(testthat::test_path("..","fixtures","source","star-fusion.tsv")),version="synthetic-v1"),
      list(input_id="splice",input_role="splice_structure",sample_id="sample-1",
        path=normalizePath(testthat::test_path("..","fixtures","source","splice-structure.tsv")),version="synthetic-v1"))
    base <- long_read_fixture_project(inputs,hla=list("patient-1"="HLA-A*02:01"))
    cfg <- base$config; cfg$peptides <- Neoantimon2:::peptide_policy(list(class_i_lengths=as.list(6:11)))
    before <- neo_run(neo_project(samples=base$samples,config=cfg,hla=lapply(base$hla,unlist),inputs=inputs))
    rows <- neo_results(before)
    after <- neo_run(recognition_clone_with(before,list(presentation_write_psms(unique(rows$ligand_sequence)),
      presentation_temp_settings(experiments=list(list(input_id="e2e-ms",ip_type="pan_hla"))))))
    cache <<- list(before=before,after=after);cache
  }
})
presentation_loh_pair <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    p <- recognition_gold_project(); cfg <- p$config
    cfg$hla_loh <- Neoantimon2:::hla_loh_config(list(sources=list(list(sample_id="sample-1",
      path=normalizePath(testthat::test_path("..","fixtures","hla_loh","lohhla-standard.HLAlossPrediction_CI.xls")),provider_version="synthetic-1"))))
    hla <- list("patient-1"=c("HLA-A*02:01","HLA-A*24:02","HLA-DRA*01:01~HLA-DRB1*04:01"))
    before <- neo_run(neo_project(samples=p$samples,config=cfg,hla=hla,inputs=recognition_clone_inputs(p)))
    after <- neo_run(recognition_clone_with(before,list(presentation_fixture_input("hla-expression.tsv","hla_expression"))))
    cache <<- list(before=before,after=after);cache
  }
})
presentation_multi_origin_gate <- function(x) {
  eligible <- Filter(function(q) q$input_status == "supported" && q$task == "antigen_processing", x$queries)
  groups <- lapply(sort(unlist(x$candidate_ids)), function(id) Filter(function(q) q$processing_query_id %in%
    vapply(Filter(function(l) l$candidate_id == id && l$applicability == "applicable", x$processing_links),
      `[[`, "", "processing_query_id"), eligible))

  found <- Filter(function(g) length(unique(vapply(g,function(q) q$context$origin_id,"")))>1L &&
    length(unique(vapply(g,function(q) paste(q$n_flank,q$c_flank,sep="|"),"")))>1L,groups)
  if (!length(found)) stop("No real pipeline candidate with distinct biological origin flanks")
  g <- found[[1L]]
  links <- Filter(function(l) l$processing_query_id %in% vapply(g, `[[`, "", "processing_query_id"), x$processing_links)
  list(candidate_id=links[[1L]]$candidate_id,queries=g,links=links,evidence=Filter(function(e)
    e$processing_query_id %in% vapply(g,`[[`,"","processing_query_id"),x$processing_evidence))
}
presentation_mhcflurry_pair <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    base <- recognition_gold_project(); cfg <- base$config
    exe <- tempfile("synthetic-mhcflurry-"); model <- tempfile("synthetic-model-");dir.create(model)
    writeLines("synthetic model identity",file.path(model,"manifest.txt"))
    writeLines(c("#!/bin/sh",'if [ "$1" = "--version" ]; then echo "mhcflurry 2.0.1"; exit 0; fi',
      'input=$1', 'output=$3',
      'echo "allele,peptide,mhcflurry_affinity,mhcflurry_affinity_percentile,mhcflurry_presentation_score,mhcflurry_presentation_percentile" > "$output"',
      'while IFS=, read -r allele peptide; do', '  if [ "$allele" = "allele" ]; then continue; fi',
      '  printf "%s,%s,80,0.9,0.77,0.7\\n" "$allele" "$peptide" >> "$output"', 'done < "$input"'),exe)
    Sys.chmod(exe,"0755")
    cfg$predictors <- c(cfg$predictors,list(list(provider="mhcflurry",task="mhc1_presentation",
      executable=exe,model_path=model,options=list())))
    before <- neo_run(neo_project(samples=base$samples,config=cfg,hla=lapply(base$hla,unlist),inputs=recognition_clone_inputs(base)))
    snapshot <- Neoantimon2:::read_runtime_snapshot(before)
    catalog <- Neoantimon2:::presentation_source_catalog(before,snapshot)
    candidates <- Neoantimon2:::presentation_candidates(before,snapshot,catalog)
    cap <- Neoantimon2:::processing_capability(list(provider="mhcflurry_processing",provider_version="2.0.1"))
    queries <- list()
    for(c in candidates) if(c$hla_class=="I") for(o in Filter(function(o)o$origin_id %in% unlist(c$origin_ids),catalog$origins))
      queries[[length(queries)+1L]] <- Neoantimon2:::processing_query(c,o,catalog$proteins[[o$sequence_context_id]],cap,catalog$boundaries[[o$sequence_context_id]])
    output <- tempfile(fileext=".tsv")
    writeLines(c("processing_query_id\tpeptide\tprocessing_score\tprocessing_percentile",
      vapply(Neoantimon2:::recognition_unique(queries,"processing_query_id"),function(q)paste(q$processing_query_id,q$peptide,.99,.01,sep="\t"),"")),output)
    after <- neo_run(recognition_clone_with(before,list(list(input_id="mhcflurry-processing-import",input_role="processing_output",
      sample_id="sample-1",path=output,version="synthetic-2.0.1"),presentation_temp_settings(list(list(provider="mhcflurry_processing",
        provider_version="2.0.1",output_input_id="mhcflurry-processing-import"))))))
    cache <<- list(before=before,after=after);cache
  }
})

# Explicit valid origin catalog shared by all candidate-granularity release gates.
presentation_test_catalog <- function(two_origins = FALSE) {
  origins <- list(presentation_origin())
  proteins <- list("protein:origin:1" = "AAAAAGILGFVFTLCCCCC")
  if (two_origins) {
    origins <- c(origins, list(presentation_origin(n = "DDDDD", c = "EEEEE", id = "origin:2")))
    proteins[["protein:origin:2"]] <- "DDDDDGILGFVFTLEEEEE"
  }
  origins <- lapply(origins, function(o) { o$origin_validity <- "VALID"; o })
  list(origins = origins, proteins = proteins, boundaries = list())
}
presentation_cross_hla <- function(provider, candidates = list(presentation_candidate(),
  presentation_candidate(hla = "HLA-B*07:02")), two_origins = FALSE, limit = Inf) {
  aliases <- if (provider %in% c("netchop", "netctlpan")) list(protein1 = "protein:origin:1") else NULL
  cap <- Neoantimon2:::processing_capability(list(provider = provider, protein_aliases = aliases))
  prepared <- Neoantimon2:::processing_prepare(candidates, presentation_test_catalog(two_origins), cap, limit)
  output <- presentation_bound_input(switch(provider, mhcflurry_processing = "mhcflurry.tsv",
    netchop = "netchop.txt", netctlpan = "netctlpan.txt", netmhcstabpan = "netmhcstabpan.txt"))
  result <- Neoantimon2:::processing_predict(prepared$queries, cap, output)
  c(prepared, result)
}
