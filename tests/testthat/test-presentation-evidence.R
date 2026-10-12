test_that("presentation overlay has independent versioned IDs and approved public exports", {
  expect_identical(Neoantimon2:::presentation_version(), "presentation-evidence-1.0.0")
  expect_identical(Neoantimon2:::presentation_id_versions()$query, "processing-query-id-v1")
  expect_identical(Neoantimon2:::presentation_id_versions()$evidence, "processing-evidence-id-v1")
  expect_setequal(getNamespaceExports("Neoantimon2"), current_public_api)
  expect_identical(Neoantimon2:::apm_component_version(), "apm-component-map-1.0.0")
})
test_that("one frozen candidate retains two origin-specific processing contexts/results", {
  a <- presentation_query(); b <- presentation_query(n = "DDDDD", c = "EEEEE", id = "origin:2")
  expect_false("candidate_id" %in% names(a))
  expect_identical(a$peptide_id, b$peptide_id)
  expect_false(identical(a$processing_query_id, b$processing_query_id))
  expect_identical(a$context$origin_id, "origin:1")
  expect_identical(b$context$origin_id, "origin:2")
  r <- Neoantimon2:::processing_predict(list(a,b), Neoantimon2:::processing_capability(list(provider = "fixture")))
  expect_length(r$evidence, 2L)
  expect_false(identical(r$evidence[[1L]]$raw_value, r$evidence[[2L]]$raw_value))
  expect_setequal(vapply(r$evidence, `[[`, "", "processing_query_id"), c(a$processing_query_id, b$processing_query_id))
})
for (mode in c("both_flanks", "n_flank_only", "c_flank_only", "no_flanks")) local({
  selected <- mode
  test_that(paste("MHCflurry processing retains biological context", selected), {
    n <- if (selected %in% c("both_flanks", "n_flank_only")) "AAAAA" else ""
    c <- if (selected %in% c("both_flanks", "c_flank_only")) "CCCCC" else ""
    q <- presentation_query("mhcflurry_processing", n = n, c = c)
    expect_identical(q$context_completeness, selected)
    expect_identical(q$n_flank, n); expect_identical(q$c_flank, c)
    r <- Neoantimon2:::processing_predict(list(q), Neoantimon2:::processing_capability(list(provider = "mhcflurry_processing")),
      presentation_bound_input("mhcflurry.tsv"))
    expect_length(r$evidence, 2L)
    expect_setequal(vapply(r$evidence, `[[`, "", "metric"), c("processing_score", "processing_percentile"))
    expect_true(all(vapply(r$evidence, function(e) is.null(e$presenter), TRUE)))
    expect_identical(q$context$n_terminus, n == "")
    expect_identical(q$context$c_terminus, c == "")
  })
})
test_that("MHCflurry never assigns a peptide-only output to origin flanks", {
  row <- list(peptide = "GILGFVFTL", processing_score = "0.9")
  q <- presentation_query("mhcflurry_processing")
  expect_false(Neoantimon2:::processing_row_matches(row, q, list(provider = "mhcflurry_processing")))
  a <- presentation_query("mhcflurry_processing", options = list(use_n_flank = FALSE, use_c_flank = FALSE))
  expect_false(identical(q$processing_query_id, a$processing_query_id))
  expect_identical(a$context_completeness, "no_flanks")
  expect_identical(a$context$n_flank, "AAAAA")
  expect_null(a$n_flank)
})
test_that("MHCflurry retains separate scores for alternative origins", {
  qs <- list(presentation_query("mhcflurry_processing"), presentation_query("mhcflurry_processing", n = "DDDDD", c = "EEEEE", id = "origin:2"))
  r <- Neoantimon2:::processing_predict(qs, Neoantimon2:::processing_capability(list(provider = "mhcflurry_processing")), presentation_bound_input("mhcflurry.tsv"))
  expect_length(r$evidence, 4L)
  scores <- Filter(function(e) e$metric == "processing_score", r$evidence)
  expect_setequal(vapply(scores, `[[`, 0, "raw_value"), c(.8,.7))
  expect_length(unique(vapply(scores, `[[`, "", "processing_query_id")), 2L)
})
for (provider in c("mhcflurry_processing", "netchop", "netctlpan", "netmhcstabpan")) local({
  p <- provider
  test_that(paste(p, "rejects Class II processing/stability without relabeling frozen EL evidence"), {
    q <- presentation_query(p, hla = "HLA-DRA*01:01~HLA-DRB1*04:01")
    candidate <- presentation_candidate(hla = "HLA-DRA*01:01~HLA-DRB1*04:01")
    prepared <- Neoantimon2:::processing_prepare(list(candidate), presentation_test_catalog(),
      Neoantimon2:::processing_capability(list(provider = p)))
    expect_identical(prepared$links[[1L]]$status, "not_applicable_class_ii")
    expect_identical(prepared$queries[[1L]]$input_status, "not_applicable_class_ii")
    expect_length(Neoantimon2:::processing_predict(prepared$queries, Neoantimon2:::processing_capability(list(provider = p)))$evidence, 0L)
  })
  test_that(paste(p, "does not convert unsupported U or nonstandard residues"), {
    for (seq in c("GILGFVFUL", "GILGFVFXL", "GILGFVFBZ")) {
      q <- presentation_query(p, peptide = seq)
      expect_identical(q$input_status, "unsupported_residue")
      expect_identical(q$peptide, seq)
      expect_length(Neoantimon2:::processing_predict(list(q), Neoantimon2:::processing_capability(list(provider = p)))$evidence, 0L)
    }
    if (p != "netmhcstabpan") expect_identical(presentation_query(p, n = "AAUAA")$input_status, "unsupported_residue")
  })
})
test_that("NetChop maps cleavage after residue to exact N and C boundaries", {
  cfg <- list(provider = "netchop", protein_aliases = list(protein1 = "protein:origin:1"), options = list(threshold = .5))
  cap <- Neoantimon2:::processing_capability(cfg)
  q <- presentation_query("netchop", protein_aliases = cfg$protein_aliases, options = cfg$options)
  r <- Neoantimon2:::processing_predict(list(q), cap, presentation_bound_input("netchop.txt"))
  expect_length(r$evidence, 2L)
  c <- Filter(function(e) e$boundary == "c_terminal", r$evidence)[[1L]]
  n <- Filter(function(e) e$boundary == "n_terminal", r$evidence)[[1L]]
  expect_equal(c$position, 14); expect_equal(n$position, 5)
  expect_equal(c$raw_value, .95); expect_equal(n$raw_value, .25)
  expect_identical(c$residue, "L"); expect_identical(c$provider_binary_call, "S")
  expect_identical(n$provider_binary_call, "."); expect_equal(c$provider_threshold, .5)
  expect_false(any(grepl("generated", unlist(lapply(r$evidence, names)))))
})
test_that("NetChop protein termini remain explicit and positions are never conflated", {
  cap <- Neoantimon2:::processing_capability(list(provider = "netchop", protein_aliases = list(terminal = "protein:origin:1")))
  q <- presentation_query("netchop", n = "", c = "", protein_aliases = list(terminal = "protein:origin:1"))
  r <- Neoantimon2:::processing_predict(list(q), cap, presentation_bound_input("netchop.txt"))
  expect_true(q$context$n_terminus); expect_true(q$context$c_terminus)
  expect_length(r$evidence, 1L); expect_identical(r$evidence[[1L]]$boundary, "c_terminal")
  expect_equal(r$evidence[[1L]]$position, 9)
  moved <- presentation_query("netchop", n = "AAAAAA")
  expect_false(identical(q$processing_query_id, moved$processing_query_id))
})
test_that("NetCTLpan retains native MHC TAP cleavage integrated and rank components", {
  cfg <- list(provider = "netctlpan", protein_aliases = list(protein1 = "protein:origin:1"))
  cap <- Neoantimon2:::processing_capability(cfg)
  q <- presentation_query("netctlpan", protein_aliases = cfg$protein_aliases)
  r <- Neoantimon2:::processing_predict(list(q), cap, presentation_bound_input("netctlpan.txt"))
  expect_length(r$evidence, 5L)
  expect_setequal(vapply(r$evidence, `[[`, "", "metric"), c("mhc_binding_score", "tap_transport_score", "c_terminal_cleavage_score", "provider_integrated_score", "integrated_percent_rank"))
  expect_equal(Filter(function(e) e$metric == "provider_integrated_score", r$evidence)[[1L]]$raw_value, .65)
  expect_equal(Filter(function(e) e$metric == "tap_transport_score", r$evidence)[[1L]]$raw_value, -.2)
  expect_identical(q$presenter, "HLA-A*02:01")
})
for (provider in c("netctlpan", "netmhcstabpan")) local({
  p <- provider
  test_that(paste(p, "enforces length bounds and preserves non-9mer caveat"), {
    lengths <- if (p == "netctlpan") 8:11 else c(8,9,14)
    for (n in lengths) expect_identical(presentation_query(p, peptide = paste(rep("A",n),collapse=""))$input_status, "supported")
    for (n in if (p == "netctlpan") c(7,12) else c(7,15))
      expect_identical(presentation_query(p, peptide = paste(rep("A",n),collapse=""))$input_status, "unsupported_length")
    expect_match(presentation_query(p)$scope_note, "approximations")
  })
})
test_that("NetMHCstabpan is pMHC-scoped and retains hours score and rank distinctly", {
  cap <- Neoantimon2:::processing_capability(list(provider = "netmhcstabpan"))
  q <- presentation_query("netmhcstabpan")
  expect_null(q$context)
  r <- Neoantimon2:::processing_predict(list(q), cap, presentation_bound_input("netmhcstabpan.txt"))
  expect_length(r$evidence, 8L)
  half <- Filter(function(e) e$metric == "half_life", r$evidence)[[1L]]
  expect_equal(half$raw_value, 3.5); expect_identical(half$unit, "hours"); expect_null(half$rank)
  rank <- Filter(function(e) e$metric == "stability_percent_rank", r$evidence)[[1L]]
  expect_equal(rank$rank, 1); expect_identical(rank$score_direction, "lower")
  expect_null(half$protein_origin_id)
  wrong <- presentation_query("netmhcstabpan", hla = "HLA-B*07:02")
  expect_length(Neoantimon2:::processing_predict(list(wrong), cap, presentation_bound_input("netmhcstabpan.txt"))$evidence, 0L)
})
test_that("source flanks are sliced only from authoritative matching coordinates", {
  o <- presentation_origin(); p <- "AAAAAGILGFVFTLCCCCC"
  c <- Neoantimon2:::processing_origin_context(o, p, 3L)
  expect_identical(c$n_flank, "AAA"); expect_identical(c$c_flank, "CCC")
  expect_equal(c$available_n_flank, 3); expect_equal(c$available_c_flank, 3)
  expect_error(Neoantimon2:::processing_origin_context(o, "CCCCCGILGFVFTAEEEEE"), "mismatch")
  expect_identical(Neoantimon2:::processing_origin_context(o, NULL)$status, "missing_context")
})
test_that("capability detection and command planning never execute providers", {
  exe <- tempfile(); writeLines("#!/bin/sh\nexit 99", exe); Sys.chmod(exe,"0755")
  for (p in c("netchop", "netctlpan", "netmhcstabpan")) {
    cap <- Neoantimon2:::processing_capability(list(provider = p, executable = exe))
    expect_identical(cap$execution_status, "local_executable_detected_not_executed")
    cmd <- Neoantimon2:::processing_command(cap, "local-input", "HLA-A02:01", 9L)
    expect_false(cmd$execute); expect_identical(cmd$executable, normalizePath(exe))
    expect_true("local-input" %in% unlist(cmd$args))
  }
  expect_error(Neoantimon2:::processing_config(list(provider="netchop", model="unknown")), "Unsupported")
  expect_error(Neoantimon2:::processing_config(list(provider="netctlpan", options=list(weight=2))), "Unknown")
})
test_that("generic PSMs retain independent quality and all provenance", {
  ms <- presentation_ms(); expect_length(ms$psms, 8L)
  p <- Filter(function(p) p$psm_source_id == "p1", ms$psms)[[1L]]
  expect_equal(p$q_value, .001); expect_null(p$pep); expect_true(p$passThreshold)
  expect_identical(p$rank, 1); expect_identical(p$decoy, FALSE)
  expect_match(p$raw_artifact$sha256, "^[a-f0-9]{64}$")
  no_q <- Filter(function(p) p$psm_source_id == "p2", ms$psms)[[1L]]
  expect_null(no_q$q_value)
  observations <- Filter(function(o) o$assignment_status == "unique_assignment" && o$modification_status == "unmodified", ms$observations)
  expect_length(observations, 1L); expect_length(observations[[1L]]$psm_ids, 2L)
  expect_identical(observations[[1L]]$origin_assignment, "unresolved_shared_peptide")
})
test_that("ambiguous spectra decoys threshold failures and unknown decoys are not unique observations", {
  ms <- presentation_ms()
  for (id in c("p3","p4")) expect_identical(Filter(function(p) p$psm_source_id == id,ms$psms)[[1L]]$assignment_status,"ambiguous_assignment")
  expect_identical(Filter(function(p) p$psm_source_id == "p6",ms$psms)[[1L]]$assignment_status,"decoy")
  expect_identical(Filter(function(p) p$psm_source_id == "p7",ms$psms)[[1L]]$assignment_status,"threshold_failed")
  expect_identical(Filter(function(p) p$psm_source_id == "p8",ms$psms)[[1L]]$assignment_status,"unresolved")
  links <- Neoantimon2:::ms_candidate_links(list(presentation_candidate()),ms$observations)
  expect_equal(sum(vapply(links,function(l) l$status == "exact_unmodified_ligand_match",TRUE)),1)
})
test_that("modified MS ligand is a backbone match and never an exact unmodified match", {
  ms <- presentation_ms()
  o <- Filter(function(o) o$modification_status == "modified",ms$observations)
  expect_length(o,1L); expect_identical(o[[1L]]$modifications[[1L]]$accession,"UNIMOD:21")
  links <- Neoantimon2:::ms_candidate_links(list(presentation_candidate()),o)
  expect_identical(links[[1L]]$status,"backbone_match_with_modification")
  expect_identical(o[[1L]]$modified_sequence,"GILGFVFT[Phospho]L")
})
test_that("pan-HLA observed ligand supports peptide without assigning either predicted allele", {
  ms <- presentation_ms()
  cs <- list(presentation_candidate(),presentation_candidate(hla="HLA-B*07:02"))
  links <- Neoantimon2:::ms_candidate_links(cs,ms$observations)
  direct <- Filter(function(l) l$status == "exact_unmodified_ligand_match",links)
  expect_length(direct,2L)
  expect_true(all(vapply(direct,function(l) l$direct_hla_assignment == "hla_unresolved",TRUE)))
  expect_true(all(vapply(ms$observations,function(o) is.null(o$restriction$presenter),TRUE)))
})
test_that("monoallelic and allele-specific controls require explicit provenance", {
  for (basis in c("monoallelic_system","allele_specific_isolation","explicit_source_assignment")) {
    exp <- presentation_experiment("monoallelic","HLA-A*02:01",basis)
    ms <- presentation_ms(exp)
    links <- Neoantimon2:::ms_candidate_links(list(presentation_candidate(),presentation_candidate(hla="HLA-B*07:02")),ms$observations)
    direct <- Filter(function(l) l$status == "exact_unmodified_ligand_match",links)
    expect_length(direct,1L); expect_identical(direct[[1L]]$direct_hla_assignment,"direct_allele_resolved")
    exp$assignment_reference <- NULL
    expect_identical(Neoantimon2:::ms_restriction(exp)$status,"hla_unresolved")
  }
  exp <- presentation_experiment(presenter="HLA-A*02:01",basis="predicted_binding")
  expect_identical(Neoantimon2:::ms_restriction(exp)$status,"hla_unresolved")
})
test_that("Class II observed ligands preserve explicit heterodimers and unresolved pan-Class-II restriction", {
  hla <- "HLA-DRA*01:01~HLA-DRB1*04:01"
  ms <- presentation_ms(presentation_experiment("class_ii_pan", class="II"))
  links <- Neoantimon2:::ms_candidate_links(list(presentation_candidate(hla=hla)),ms$observations)
  direct <- Filter(function(l) l$status == "exact_unmodified_ligand_match",links)
  expect_length(direct,1L); expect_identical(direct[[1L]]$direct_hla_assignment,"hla_unresolved")
  exp <- presentation_experiment("allele_specific",hla,"explicit_source_assignment",class="II")
  expect_identical(Neoantimon2:::ms_restriction(exp)$presenter,hla)
  exp$presenter <- "HLA-DPB1*04:01"
  expect_identical(Neoantimon2:::ms_restriction(exp)$status,"hla_unresolved")
})
test_that("MS context and sequence matching do not transfer evidence or use I/L equivalence", {
  ms <- presentation_ms()
  for (c in list(presentation_candidate(sample="normal"),presentation_candidate(peptide="GLLGFVFTL")))
    expect_length(Neoantimon2:::ms_candidate_links(list(c),ms$observations),0L)
  os <- lapply(ms$observations,function(o) { o$context$timepoint <- "post-treatment";o })
  expect_length(Neoantimon2:::ms_candidate_links(list(presentation_candidate()),os),0L)
})
test_that("multiple MS runs remain distinct observations without counting probability", {
  a <- presentation_ms(); b <- presentation_ms(within(presentation_experiment(),ms_run <- "run-2"))
  expect_false(identical(a$psms[[1L]]$psm_id,b$psms[[1L]]$psm_id))
  expect_false(identical(a$observations[[1L]]$immunopeptidomics_observation_id,b$observations[[1L]]$immunopeptidomics_observation_id))
  pick <- function(ms) Filter(function(o) o$peptide_sequence == "GILGFVFTL" && o$assignment_status == "unique_assignment" &&
    o$modification_status == "unmodified",ms$observations)[[1L]]
  expect_identical(pick(a)$biological_peptide_id,pick(b)$biological_peptide_id)
})
test_that("mzIdentML semantic extraction preserves software databases spectra peptide evidence and CVs", {
  ms <- presentation_ms(file="identifications.mzid",role="immunopeptidomics_mzidentml")
  expect_length(ms$psms,3L)
  p <- Filter(function(p) p$psm_source_id == "sii1",ms$psms)[[1L]]
  expect_equal(p$q_value,.001); expect_equal(p$charge_state,2)
  expect_equal(p$experimental_mz,500.3); expect_equal(p$calculated_mz,500.2)
  expect_identical(p$analysis_software[[1L]]$attributes$version,"1.0")
  expect_identical(p$search_database[[1L]]$attributes$location,"synthetic:database.fasta")
  expect_identical(p$search_database[[1L]]$cv_params[[1L]]$value,
    digest::digest(file=presentation_fixture_path("search-database.fasta"),algo="sha256",serialize=FALSE))
  expect_identical(p$protein_accessions[[1L]],"synthetic:protein1")
  expect_identical(p$peptide_source_id,"peptide1")
  expect_identical(p$spectrum_file$spectra_data_ref,"spectra1")
  expect_setequal(vapply(p$cv_params,`[[`,"","accession"),c("MS:1002354","MS:1001171","MS:1002351"))
  expect_equal(p$fdr_local,.003)
  expect_null(p$fdr)
  modified <- Filter(function(p) p$psm_source_id == "sii2",ms$psms)[[1L]]
  expect_null(modified$q_value); expect_equal(modified$pep,.002)
  expect_identical(modified$modifications[[1L]]$cv_params[[1L]]$accession,"UNIMOD:21")
  decoy <- Filter(function(p) p$psm_source_id == "sii3",ms$psms)[[1L]]
  expect_true(decoy$decoy); expect_false(decoy$passThreshold)
})
test_that("HLA expression preserves imbalance zero missing and per-chain measurements", {
  x <- Neoantimon2:::apm_import(presentation_bound_input("hla-expression.tsv","hla_expression"))
  expect_length(x,7L)
  pick <- function(a) Filter(function(e) e$allele == a,x)[[1L]]
  expect_equal(pick("HLA-A*02:01")$allele_fraction,.999)
  expect_equal(pick("HLA-A*01:01")$raw_count,1)
  expect_equal(pick("HLA-B*07:02")$abundance,0)
  expect_null(pick("HLA-B*08:01")$abundance)
  expect_true(all(vapply(x,function(e) e$state == "evidence_present_unresolved",TRUE)))
  expect_true(all(vapply(x,function(e) is.null(e$provider_loss_call),TRUE)))
  ii <- Neoantimon2:::apm_candidate_links(list(presentation_candidate(hla="HLA-DRA*01:01~HLA-DRB1*04:01")),x)
  applicable <- Filter(function(l) l$applicability == "applicable",ii)
  expect_length(applicable,2L); expect_setequal(vapply(applicable,`[[`,"","gene"),c("HLA-DRA","HLA-DRB1"))
})
test_that("provider ASE interpretation remains provider supplied", {
  x <- Neoantimon2:::apm_import(presentation_bound_input("hla-provider-call.tsv","hla_expression"))[[1L]]
  expect_identical(x$state,"provider_loss_call")
  expect_identical(x$provider_loss_call,"ASE_loss")
  expect_identical(x$interpretation_reference,"synthetic:provider-model")
})
for (file in c("arcas-alleles.tsv","arcas-homozygous.tsv","arcas-genes.tsv")) local({
  f <- file
  test_that(paste("arcasHLA quant uses documented named headers",f), {
    x <- Neoantimon2:::apm_import(presentation_bound_input(f,"arcas_hla_quant"))
    expect_true(length(x)>0)
    expect_true(all(vapply(x,function(e) e$provider == "arcasHLA",TRUE)))
    expect_true(all(vapply(x,function(e) e$unit == "TPM",TRUE)))
    expect_true(all(vapply(x,function(e) e$state == "evidence_present_unresolved",TRUE)))
    expect_true(all(vapply(x,function(e) !is.null(e$raw_fields$provider_fields),TRUE)))
    expect_true(all(vapply(x,function(e) is.null(e$provider_loss_call),TRUE)))
    if (f == "arcas-alleles.tsv") {
      expect_equal(x[[1L]]$raw_count,999);expect_equal(x[[2L]]$abundance,.1)
      expect_identical(x[[1L]]$raw_fields$provider_fields$baf,"0.0")
      expect_null(x[[1L]]$allele_fraction)
    }
  })
})
test_that("APM map preserves independent roles without broad IFN loss inference", {
  x <- Neoantimon2:::apm_import(presentation_bound_input("apm.tsv","apm_evidence"))
  genes <- unique(vapply(x,`[[`,"","gene"))
  expect_true(all(c("B2M","TAP1","TAP2","TAPBP","CALR","PDIA3","ERAP1","ERAP2","PSMB8","PSMB9","PSMB10","NLRC5","CD74","HLA-DMA","CIITA") %in% genes))
  b2m <- Filter(function(e) e$gene == "B2M",x)
  expect_length(b2m,4L)
  missense <- Filter(function(e) e$measurement_type == "somatic_alteration",b2m)[[1L]]
  expect_identical(missense$state,"evidence_present_unresolved")
  expect_identical(missense$consequence,"missense_variant")
  expr <- Filter(function(e) e$measurement_type == "gene_expression",b2m)
  expect_true(all(vapply(expr,function(e) e$state == "evidence_present_unresolved",TRUE)))
  strong <- Filter(function(e) e$measurement_type == "copy_number",b2m)[[1L]]
  expect_identical(strong$state,"direct_loss_supported")
  expect_equal(strong$copy_number,0)
  expect_identical(Filter(function(e) e$gene == "JAK1",x)[[1L]]$role,"outside_core_apm_map")
})
test_that("B2M and TAP evidence apply only to Class I; Class II loading is separate", {
  x <- Neoantimon2:::apm_import(presentation_bound_input("apm.tsv","apm_evidence"))
  i <- presentation_candidate(); ii <- presentation_candidate(hla="HLA-DRA*01:01~HLA-DRB1*04:01")
  links <- Neoantimon2:::apm_candidate_links(list(i,ii),x)
  for (id in c(i$candidate_id,ii$candidate_id)) {
    b <- Filter(function(l) l$candidate_id == id && l$gene %in% c("B2M","TAP1","TAP2"),links)
    expect_true(all(vapply(b,function(l) l$applicability == if (id == i$candidate_id) "applicable" else "not_applicable",TRUE)))
  }
  classii <- Filter(function(l) l$candidate_id == ii$candidate_id && l$gene %in% c("CD74","HLA-DMA","CIITA"),links)
  expect_true(all(vapply(classii,function(l) l$applicability == "applicable",TRUE)))
})
test_that("absence of supplied machinery evidence never reports pathway intact", {
  view <- Neoantimon2:::presentation_candidate_view(NULL,presentation_candidate()$candidate_id)
  expect_identical(view$summary$apm_evidence_status,"not_evaluated")
  expect_identical(view$summary$presenter_expression_status,"not_evaluated")
  expect_identical(view$summary$immunopeptidomics_status,"no_evidence_supplied")
  expect_false(any(unlist(view$summary) %in% c("intact","normal","no_escape")))
})
test_that("presentation E2E preserves frozen candidate IDs eligibility grades ranks and priority origins", {
  p <- presentation_e2e_pair()
  a <- neo_results(p$before); b <- neo_results(p$after)
  expect_true(nrow(a)>0)
  expect_identical(recognition_invariance_fields(a),recognition_invariance_fields(b))
  expect_identical(nrow(a),nrow(b))
  expect_true(all(b$exact_observed_ligand_count>0))
  expect_true(all(b$ms_psm_count>0))
  expect_true(all(b$direct_hla_assignment_status == "hla_unresolved"))
  expect_true(all(b$processing_evidence_count[b$hla_class == "I"]>0))
  expect_true(all(b$pmhc_stability_evidence_count[b$hla_class == "I"]>0))
  expect_true(all(b$processing_evidence_count[b$hla_class == "II"] == 0))
  expect_true(all(b$apm_evidence_status[b$hla_class == "I"] == "direct_loss_supported"))
  expect_true(all(b$apm_evidence_status[b$hla_class == "II"] != "direct_loss_supported"))
  x <- Neoantimon2:::presentation_read(p$after)
  expect_true(Neoantimon2:::presentation_validate(x))
  expect_true(all(vapply(x$provider_runs,function(r) r$external_runs == 0L,TRUE)))
  # Strong additive sidecars leave all existing pMHC prediction rows unchanged.
  pred <- function(p) unlist(lapply(Neoantimon2:::somatic_prediction_runs(p), `[[`, "evidence"), recursive=FALSE)
  expect_identical(pred(p$before),pred(p$after))
})
test_that("results filters and structured/text explanations expose separate compact families", {
  p <- presentation_e2e_pair(); rows <- neo_results(p$after)
  id <- rows$candidate_id[[1L]]
  expect_equal(nrow(neo_results(p$after,candidate_id=id)),1)
  e <- neo_explain(p$after,candidate_id=id)$candidates[[1L]]$presentation_evidence
  expect_true(all(c("antigen_processing","pmhc_stability","immunopeptidomics","hla_expression","antigen_presentation_machinery") %in% names(e)))
  text <- neo_explain(p$after,candidate_id=id,format="text")
  for (section in c("ANTIGEN PROCESSING","PMHC BINDING / PRESENTATION","PMHC STABILITY","IMMUNOPEPTIDOMICS","HLA AVAILABILITY / EXPRESSION","ANTIGEN-PRESENTATION MACHINERY","HLA-SPECIFIC LOH"))
    expect_match(text,section,fixed=TRUE)
  expect_false(grepl("definitely presented|immunogenicity_score|APM_score",text))
})
test_that("presentation checkpoints verify raw artifact and sidecar integrity", {
  p <- presentation_e2e_pair(); x <- Neoantimon2:::presentation_read(p$after)
  mutated <- x; mutated$processing_evidence[[1L]]$raw_value <- 500
  expect_error(Neoantimon2:::presentation_validate(mutated),"identity mismatch|fingerprint mismatch")
  mutated <- x; mutated$queries[[1L]]$n_flank <- "wrong"
  expect_error(Neoantimon2:::presentation_validate(mutated),"identity mismatch|fingerprint mismatch")
  expect_identical(x$run_id,p$after$run_id)
})
test_that("bundled source candidates stay NOT_GRADED despite direct MS and processing evidence", {
  p <- presentation_source_pair(); a <- neo_results(p$before); b <- neo_results(p$after)
  expect_identical(recognition_invariance_fields(a),recognition_invariance_fields(b))
  expect_true(all(b$final_grade == "NOT_GRADED"))
  expect_true(all(b$exact_observed_ligand_count>0))
  expect_true(any(grepl("fusion",b$source_type)))
  expect_true(any(grepl("novel_junction",b$source_type)))
  expect_true(any(grepl("long_read_isoform",b$source_type)))
  x <- Neoantimon2:::presentation_read(p$after)
  expect_true(all(vapply(x$ms_links,function(l) l$origin_assignment == "unresolved_shared_peptide",TRUE)))
  expect_true(all(vapply(Filter(function(q) q$input_status == "supported",x$queries),function(q)
    !is.null(q$context$protein_sequence),TRUE)))
  expect_true(any(b$processing_origin_count>1))
})
test_that("HLA DNA LOH and RNA expression remain parallel nonconcordance-tolerant evidence", {
  p <- presentation_loh_pair(); a <- neo_results(p$before); b <- neo_results(p$after)
  expect_identical(recognition_invariance_fields(a),recognition_invariance_fields(b))
  expect_identical(a$presenter_hla_loh_status,b$presenter_hla_loh_status)
  rows <- b[b$mhc_key == "HLA-A*02:01",]
  expect_true(nrow(rows)>0)
  expect_true(all(rows$presenter_expression_status == "evidence_present"))
  expect_true(all(rows$presenter_hla_loh_status == "evaluated"))
  expect_true(all(rows$presenter_hla_loh_call == "loss_supported"))
  e <- neo_explain(p$after,candidate_id=rows$candidate_id[[1L]])$candidates[[1L]]
  expect_equal(Filter(function(x) x$allele == "HLA-A*02:01",e$presentation_evidence$hla_expression)[[1L]]$abundance,99.9)
  expect_identical(e$hla_loh$allele_evidence$loh_call,"loss_supported")
})
test_that("query-limit status is explicit with no partial execution or silent truncation", {
  base <- presentation_e2e_pair()$before
  p <- neo_run(recognition_clone_with(base,list(presentation_temp_settings(max_queries=1L))))
  x <- Neoantimon2:::presentation_read(p)
  expect_true(length(x$queries)>1)
  expect_length(x$processing_evidence,0L)
  expect_identical(x$provider_runs[[1L]]$status,"not_run/query_limit")
  expect_true(all(vapply(x$queries,function(q) q$input_status %in%
    c("not_run/query_limit", "not_applicable_class_ii"),TRUE)))
  expect_true(all(vapply(Filter(function(l) l$applicability == "not_applicable", x$processing_links),
    function(l) l$status == "not_applicable_class_ii", TRUE)))
})
test_that("APM requires documented support for a direct loss state", {
  input <- presentation_bound_input("apm.tsv","apm_evidence")
  rows <- Neoantimon2:::presentation_table(input$resolved_path,c("gene","measurement_type"))
  row <- Filter(function(r) identical(r$loss_support,"documented_homozygous_deletion"),rows)[[1L]]
  row$interpretation_reference <- NULL
  path <- tempfile(fileext=".tsv")
  fields <- names(row);writeLines(c(paste(fields,collapse="\t"),paste(vapply(row,function(x)if(is.null(x))"" else x,""),collapse="\t")),path)
  input$resolved_path <- path
  x <- Neoantimon2:::apm_import(input)[[1L]]
  expect_identical(x$state,"evidence_present_unresolved")
  expect_equal(x$copy_number,0)
})
test_that("PSM named sample/context metadata is preserved and contradictions rejected", {
  input <- presentation_bound_input("psms.tsv","immunopeptidomics_psm")
  row <- list(psm_id="context-psm",spectrum_id="scan=123",peptide_sequence="GILGFVFTL",decoy=FALSE,
    sample_id="sample-1",timepoint="post",tissue="tumor")
  p <- Neoantimon2:::ms_psm(row,input,presentation_experiment())
  expect_identical(p$context$timepoint,"post");expect_identical(p$context$tissue,"tumor")
  row$sample_id <- "normal"
  expect_error(Neoantimon2:::ms_psm(row,input,presentation_experiment()),"sample binding conflict")
})
test_that("mzIdentML ignores definition ordering and rejects dangling refs and entities", {
  input <- presentation_bound_input("identifications.mzid","immunopeptidomics_mzidentml")
  baseline <- Neoantimon2:::ms_mzidentml_import(input,presentation_experiment())
  path <- tempfile(fileext=".mzid")
  text <- paste(readLines(input$resolved_path),collapse="\n")
  writeLines(sub('peptideEvidence_ref="pe1"','peptideEvidence_ref="missing"',text,fixed=TRUE),path)
  input$resolved_path <- path
  expect_error(Neoantimon2:::ms_mzidentml_import(input,presentation_experiment()),"Unresolved peptide evidence")
  writeLines('<!DOCTYPE MzIdentML [<!ENTITY test SYSTEM "file:///etc/passwd">]>',path)
  expect_error(Neoantimon2:::ms_mzidentml_import(input,presentation_experiment()),"DTD/entity")
  expect_equal(length(baseline),3)
})
test_that("processing fingerprints change with source position HLA and provider options", {
  q <- presentation_query("netctlpan")
  for (other in list(presentation_query("netctlpan",n="AAAAAA"),presentation_query("netctlpan",hla="HLA-B*07:02"),
    presentation_query("netctlpan",provider_version="1.1b"))) {
    expect_false(identical(q$processing_query_id,other$processing_query_id))
    expect_false(identical(q$fingerprint,other$fingerprint))
  }
  cap <- Neoantimon2:::processing_capability(list(provider="fixture"))
  a <- Neoantimon2:::processing_predict(list(presentation_query(),presentation_query(id="origin:2")),cap)
  b <- Neoantimon2:::processing_predict(rev(list(presentation_query(),presentation_query(id="origin:2"))),cap)
  expect_identical(a$evidence,b$evidence)
})
test_that("no supplied APM evidence in a real MS-only run stays not evaluated", {
  p <- neo_run(recognition_clone_with(presentation_e2e_pair()$before,
    list(presentation_write_psms("AAAAAAAAA"))))
  rows <- neo_results(p)
  expect_true(all(rows$apm_evidence_status == "not_evaluated"))
  expect_true(all(rows$presenter_expression_status == "not_evaluated"))
  expect_true(all(rows$exact_observed_ligand_count == 0))
})
test_that("a bounded protein end is never claimed as a deterministic terminus", {
  o <- presentation_origin(n="",c="")
  meta <- list(protein_origin_id="real-protein-origin",n_terminus_known=TRUE,c_terminus_known=FALSE,
    translation_boundary="TRANSCRIPT_END_NO_STOP")
  context <- Neoantimon2:::processing_origin_context(o,"GILGFVFTL",source_metadata=meta)
  expect_true(context$n_terminus);expect_null(context$c_terminus)
  expect_identical(context$protein_origin_id,"real-protein-origin")
  cap <- Neoantimon2:::processing_capability(list(provider="netchop"))
  q <- Neoantimon2:::processing_query(presentation_candidate(),o,"GILGFVFTL",cap,meta)
  expect_identical(q$input_status,"missing_context")
})
test_that("shared candidate overlay retains all cross-source origins", {
  p <- presentation_source_pair()$after
  snap <- Neoantimon2:::read_runtime_snapshot(p)
  catalog <- Neoantimon2:::presentation_source_catalog(p,snap)
  cs <- Neoantimon2:::presentation_candidates(p,snap,catalog)
  for (c in cs) {
    own <- Filter(function(o) identical(o$sample_peptide_id,c$sample_peptide_id),catalog$origins)
    expect_true(all(vapply(own,`[[`,"","origin_id") %in% unlist(c$origin_ids)))
  }
})
test_that("real frozen candidate has distinct valid biological origins and flanked processing queries", {
  p <- presentation_source_pair()$after
  x <- Neoantimon2:::presentation_read(p)
  gate <- presentation_multi_origin_gate(x)
  expect_length(gate$queries, 6L)
  expect_length(gate$evidence, 6L)
  expect_length(unique(vapply(gate$links,`[[`,"","candidate_id")),1L)
  expect_length(unique(vapply(gate$queries,`[[`,"","peptide_id")),1L)
  expect_true(length(unique(vapply(gate$queries,function(q) q$context$origin_id,"")))>1L)
  expect_true(length(unique(vapply(gate$queries,function(q) paste(q$n_flank,q$c_flank,sep="|"),"")))>1L)
  catalog <- Neoantimon2:::presentation_source_catalog(p,Neoantimon2:::read_runtime_snapshot(p))
  origins <- Filter(function(o) o$origin_id %in% vapply(gate$queries,function(q)q$context$origin_id,""),catalog$origins)
  expect_true(all(vapply(origins,function(o) o$origin_validity == "VALID",TRUE)))
  expect_true(length(unique(vapply(gate$evidence,`[[`,0,"raw_value")))>1L)
})
test_that("new presentation evidence leaves frozen recognition records unchanged", {
  p <- recognition_e2e_pair()$after
  inputs <- lapply(p$inputs,function(i) list(input_id=i$input_id,input_role=i$input_role,sample_id=i$sample_id,
    path=i$declared_path,version=i$version))
  after <- neo_run(neo_project(samples=p$samples,config=p$config,hla=lapply(p$hla,unlist),
    inputs=c(inputs,list(presentation_fixture_input("apm.tsv","apm_evidence"))),fixture_kind=p$fixture_kind))
  a <- Neoantimon2:::recognition_read(p); b <- Neoantimon2:::recognition_read(after)
  scrub <- function(rows) lapply(rows,function(e) {e$run_id<-NULL;e})
  expect_identical(scrub(a$evidence),scrub(b$evidence))
  expect_identical(a$queries,b$queries)
  expect_identical(a$repertoire,b$repertoire)
  expect_identical(recognition_invariance_fields(neo_results(p)),recognition_invariance_fields(neo_results(after)))
})
test_that("strong MHCflurry processing leaves affinity and presentation rows exactly unchanged", {
  p <- presentation_mhcflurry_pair()
  extract <- function(p) Filter(function(e)e$predictor=="mhcflurry",
    unlist(lapply(Neoantimon2:::somatic_prediction_runs(p),`[[`,"evidence"),recursive=FALSE))
  a <- extract(p$before); b <- extract(p$after)
  expect_true(length(a)>0)
  expect_setequal(unique(vapply(a,`[[`,"","metric")),c("BA_nM","BA_rank","presentation_score","presentation_percentile"))
  expect_identical(a,b)
  expect_identical(recognition_invariance_fields(neo_results(p$before)),recognition_invariance_fields(neo_results(p$after)))
  evidence <- Neoantimon2:::presentation_read(p$after)$processing_evidence
  expect_true(length(evidence)>0)
  expect_true(all(vapply(Filter(function(e)e$metric=="processing_score",evidence),function(e)e$raw_value==.99,TRUE)))
  expect_false(any(vapply(evidence,function(e)e$task=="mhc1_binding",TRUE)))
})
test_that("native NetCTLpan non-9-mer component outputs retain provider scope and coordinates", {
  cap <- Neoantimon2:::processing_capability(list(provider="netctlpan",protein_aliases=list(protein1="protein:origin:1")))
  for(peptide in c("GILGFVFT","GILGFVFTLA","GILGFVFTLAA")) {
    q <- presentation_query("netctlpan",peptide=peptide,protein_aliases=list(protein1="protein:origin:1"))
    r <- Neoantimon2:::processing_predict(list(q),cap,presentation_bound_input("netctlpan.txt"))
    expect_length(r$evidence,5L)
    expect_true(all(vapply(r$evidence,function(e) grepl("approximations",e$scope_note),TRUE)))
    expect_equal(q$context$peptide_start_aa,6)
  }
})
test_that("native NetMHCstabpan 8-mer and 14-mer outputs keep stability and rank metrics", {
  cap <- Neoantimon2:::processing_capability(list(provider="netmhcstabpan"))
  for(peptide in c("GILGFVFT","GILGFVFTLAAAAA")) {
    q <- presentation_query("netmhcstabpan",peptide=peptide)
    r <- Neoantimon2:::processing_predict(list(q),cap,presentation_bound_input("netmhcstabpan.txt"))
    expect_length(r$evidence,8L)
    half <- Filter(function(e)e$metric=="half_life",r$evidence)[[1L]]
    expect_identical(half$unit,"hours");expect_null(half$rank)
    expect_true(any(vapply(r$evidence,function(e)e$metric=="stability_percent_rank",TRUE)))
  }
})
test_that("analysis artifact references change provenance fingerprints but not biological processing queries", {
  c <- presentation_candidate();o <- presentation_origin()
  cap <- Neoantimon2:::processing_capability(list(provider="mhcflurry_processing"))
  meta <- list(protein_origin_id="protein-sha256:stable",source_record_id="translation:run-one",n_terminus_known=TRUE,c_terminus_known=TRUE)
  a <- Neoantimon2:::processing_query(c,o,"AAAAAGILGFVFTLCCCCC",cap,meta)
  meta$source_record_id <- "translation:run-two";o$mutant_sequence_ref <- "another-analysis#protein"
  b <- Neoantimon2:::processing_query(c,o,"AAAAAGILGFVFTLCCCCC",cap,meta)
  expect_identical(a$processing_query_id,b$processing_query_id)
  expect_false(identical(a$fingerprint,b$fingerprint))
  row <- list(processing_query_id=a$processing_query_id,peptide="WRONGSEQ")
  expect_false(Neoantimon2:::processing_row_matches(row,a,cap$config))
})

for (provider in c("mhcflurry_processing", "netchop")) local({
  p <- provider
  test_that(paste(p, "shares one cross-HLA computation with two distinct frozen candidate links"), {
    x <- presentation_cross_hla(p)
    expect_length(x$queries, 1L); expect_length(x$links, 2L)
    expect_length(unique(vapply(x$links, `[[`, "", "candidate_id")), 2L)
    expect_true(all(vapply(x$links, function(l) l$applicability == "applicable", TRUE)))
    expect_length(x$evidence, 2L)
    expect_length(unique(vapply(x$evidence, `[[`, "", "processing_query_id")), 1L)
    expect_false(any(vapply(x$evidence, function(e) "candidate_id" %in% names(e), TRUE)))
    for (candidate in list(presentation_candidate(), presentation_candidate(hla = "HLA-B*07:02"))) {
      sidecar <- list(processing_links = x$links, queries = x$queries, processing_evidence = x$evidence, statuses = x$statuses)
      view <- Neoantimon2:::presentation_candidate_view(sidecar, candidate$candidate_id)
      expect_equal(view$summary$processing_evidence_count, 2L)
      expect_equal(view$summary$processing_origin_count, 1L)
    }
  })
})
test_that("NetCTLpan separates presenters without candidate IDs in computation identity", {
  x <- presentation_cross_hla("netctlpan")
  expect_length(x$queries, 2L); expect_length(x$links, 2L)
  expect_setequal(vapply(x$queries, `[[`, "", "presenter"), c("HLA-A*02:01", "HLA-B*07:02"))
  expect_false(any(vapply(x$queries, function(q) "candidate_id" %in% names(q), TRUE)))
  candidate <- presentation_candidate(); candidate$candidate_id <- "another-wrapper"; candidate$sample_id <- "another-sample"
  cap <- Neoantimon2:::processing_capability(list(provider = "netctlpan"))
  q <- Neoantimon2:::processing_query(candidate, presentation_origin(), "AAAAAGILGFVFTLCCCCC", cap)
  expect_identical(q$processing_query_id, presentation_query("netctlpan")$processing_query_id)
})
test_that("NetMHCstabpan shares one peptide-presenter query across distinct valid origin paths", {
  x <- presentation_cross_hla("netmhcstabpan", candidates = list(presentation_candidate()), two_origins = TRUE)
  expect_length(x$queries, 1L); expect_length(x$links, 1L)
  expect_length(x$links[[1L]]$origin_contexts, 2L)
  expect_length(x$evidence, 8L)
  expect_null(x$queries[[1L]]$context)
  expect_identical(presentation_query("netmhcstabpan")$processing_query_id,
    presentation_query("netmhcstabpan", n = "DDDDD", c = "EEEEE", id = "origin:2")$processing_query_id)
})
test_that("NetMHCstabpan keeps different canonical presenters as two computations", {
  x <- presentation_cross_hla("netmhcstabpan", two_origins = TRUE)
  expect_length(x$queries, 2L); expect_length(x$links, 2L)
  expect_setequal(vapply(x$queries, `[[`, "", "presenter"), c("HLA-A*02:01", "HLA-B*07:02"))
  expect_true(all(vapply(x$links, function(l) length(l$origin_contexts) == 2L, TRUE)))
})
test_that("candidate HLA origin and link construction order cannot change IDs or link sets", {
  candidates <- list(presentation_candidate(), presentation_candidate(hla = "HLA-B*07:02"))
  for (provider in c("mhcflurry_processing", "netchop", "netctlpan", "netmhcstabpan")) {
    cap <- Neoantimon2:::processing_capability(list(provider = provider,
      protein_aliases = list(protein1 = "protein:origin:1")))
    catalog <- presentation_test_catalog(TRUE)
    a <- Neoantimon2:::processing_prepare(candidates, catalog, cap)
    catalog$origins <- rev(catalog$origins)
    b <- Neoantimon2:::processing_prepare(rev(candidates), catalog, cap)
    expect_identical(a, b)
    expect_identical(Neoantimon2:::sorted_records(rev(a$links)), b$links)
    # Synthetic metric output exercises evidence identity for every prepared query.
    output <- tempfile(fileext = if (provider == "mhcflurry_processing") ".tsv" else ".txt")
    if (provider == "mhcflurry_processing") {
      writeLines(c("processing_query_id\tpeptide\tprocessing_score", vapply(a$queries,
        function(q) paste(q$processing_query_id, q$peptide, .8, sep = "\t"), "")), output)
    } else file.copy(presentation_fixture_path(paste0(provider, ".txt")), output)
    bound <- list(resolved_path = output, input_id = "order-fixture", version = "synthetic")
    ea <- Neoantimon2:::processing_predict(a$queries, cap, bound)
    eb <- Neoantimon2:::processing_predict(rev(b$queries), cap, bound)
    expect_identical(ea, eb)
    expect_true(length(ea$evidence) > 0L)
  }
})
test_that("adding a candidate changes only links for HLA-independent queries and evidence", {
  for (provider in c("mhcflurry_processing", "netchop")) {
    a <- presentation_cross_hla(provider, candidates = list(presentation_candidate()))
    b <- presentation_cross_hla(provider)
    expect_identical(a$queries, b$queries)
    expect_identical(a$evidence, b$evidence)
    expect_length(a$links, 1L); expect_length(b$links, 2L)
    candidate_c <- presentation_candidate(); candidate_c$candidate_id <- "third-wrapper"
    c <- presentation_cross_hla(provider, candidates = list(presentation_candidate(),
      presentation_candidate(hla = "HLA-B*07:02"), candidate_c))
    expect_identical(b$queries, c$queries); expect_identical(b$evidence, c$evidence)
    expect_length(c$links, 3L)
  }
})
test_that("query limit counts unique computations after candidate deduplication", {
  a <- presentation_cross_hla("mhcflurry_processing", limit = 1)
  expect_equal(a$requested_queries, 1L); expect_false(a$limited)
  expect_length(a$queries, 1L); expect_length(a$evidence, 2L)
  b <- presentation_cross_hla("netctlpan", limit = 1)
  expect_equal(b$requested_queries, 2L); expect_true(b$limited)
  expect_length(b$queries, 2L); expect_length(b$evidence, 0L)
  expect_true(all(vapply(b$queries, function(q) q$input_status == "not_run/query_limit", TRUE)))
  expect_true(all(vapply(b$statuses, function(s) s$status == "not_run/query_limit", TRUE)))
  for (provider in c("mhcflurry_processing", "netctlpan")) {
    project <- recognition_fixture_project(list(presentation_temp_settings(
      providers = list(list(provider = provider)), max_queries = 1L)))
    x <- Neoantimon2:::presentation_construct(project,
      list(presentation_candidate(), presentation_candidate(hla = "HLA-B*07:02")), presentation_test_catalog())
    Neoantimon2:::presentation_validate(x)
    expect_equal(x$provider_runs[[1L]]$requested_queries, if (provider == "netctlpan") 2L else 1L)
    expect_identical(x$provider_runs[[1L]]$status, if (provider == "netctlpan") "not_run/query_limit" else "prepared")
  }
})
test_that("shared Class I queries cannot attach processing evidence to Class II links", {
  ii <- presentation_candidate(hla = "HLA-DRA*01:01~HLA-DRB1*04:01")
  for (provider in c("mhcflurry_processing", "netchop")) {
    x <- presentation_cross_hla(provider, candidates = list(ii, presentation_candidate()), limit = 1)
    expect_length(x$queries, 1L); expect_length(x$links, 2L)
    expect_equal(x$requested_queries, 1L); expect_length(x$evidence, 2L)
    sidecar <- list(processing_links = x$links, queries = x$queries, processing_evidence = x$evidence, statuses = x$statuses)
    view <- Neoantimon2:::presentation_candidate_view(sidecar, ii$candidate_id)
    expect_equal(view$summary$processing_evidence_count, 0L)
    expect_identical(view$processing_statuses[[1L]]$status, "not_applicable_class_ii")
  }
})
test_that("unsupported biological inputs share queries while retaining separate candidate links", {
  for (peptide in c("GILGFVFUL", "GILGFVFTLAAAAAAAA")) {
    candidates <- list(presentation_candidate(peptide), presentation_candidate(peptide, "HLA-B*07:02"))
    origin <- presentation_origin(peptide)
    catalog <- list(origins = list(origin), proteins = setNames(list(paste0("AAAAA", peptide, "CCCCC")), origin$sequence_context_id), boundaries = list())
    cap <- Neoantimon2:::processing_capability(list(provider = "mhcflurry_processing"))
    x <- Neoantimon2:::processing_prepare(candidates, catalog, cap)
    expect_length(x$queries, 1L); expect_length(x$links, 2L)
    expect_identical(x$queries[[1L]]$input_status, if (grepl("U", peptide)) "unsupported_residue" else "unsupported_length")
  }
  cap <- Neoantimon2:::processing_capability(list(provider = "mhcflurry_processing"))
  catalog <- presentation_test_catalog(); catalog$proteins <- list()
  x <- Neoantimon2:::processing_prepare(list(presentation_candidate(), presentation_candidate(hla = "HLA-B*07:02")), catalog, cap)
  expect_length(x$queries, 1L); expect_length(x$links, 2L)
  expect_identical(x$queries[[1L]]$input_status, "missing_context")
})

test_that("computational identity excludes unused options and result provenance wrappers", {
  a <- presentation_query("netchop", protein_aliases = list(protein1 = "protein:origin:1"))
  b <- presentation_query("netchop", protein_aliases = list(protein1 = "protein:origin:1"),
    options = list(use_n_flank = FALSE, use_c_flank = FALSE, flank_length = 1))
  expect_identical(a$processing_query_id, b$processing_query_id)
  cap <- Neoantimon2:::processing_capability(list(provider = "netchop", protein_aliases = list(protein1 = "protein:origin:1")))
  output <- presentation_bound_input("netchop.txt")
  ea <- Neoantimon2:::processing_predict(list(a), cap, output)$evidence
  eb <- Neoantimon2:::processing_predict(list(b), cap, output)$evidence
  expect_identical(vapply(ea, `[[`, "", "processing_evidence_id"), vapply(eb, `[[`, "", "processing_evidence_id"))
  expect_false(identical(ea[[1L]]$fingerprint, eb[[1L]]$fingerprint))
  for (provider in c("mhcflurry_processing", "netctlpan", "netmhcstabpan"))
    expect_identical(presentation_query(provider)$processing_query_id,
      presentation_query(provider, options = list(threshold = .9))$processing_query_id)
  output$input_id <- "another-artifact-wrapper"; output$sample_id <- "another-sample"
  ec <- Neoantimon2:::processing_predict(list(a), cap, output)$evidence
  expect_identical(vapply(ea, `[[`, "", "processing_evidence_id"), vapply(ec, `[[`, "", "processing_evidence_id"))
  expect_false(identical(ea[[1L]]$fingerprint, ec[[1L]]$fingerprint))
  expect_false(identical(a$processing_query_id,
    presentation_query("netchop", options = list(threshold = .9))$processing_query_id))
})
