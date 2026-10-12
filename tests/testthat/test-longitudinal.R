test_that("longitudinal contract remains separate and public API is explicitly approved", {
  expect_identical(Neoantimon2:::longitudinal_version(), "longitudinal-integration-1.0.0")
  expect_setequal(unlist(Neoantimon2:::longitudinal_namespaces()), c("subject-context-id-v1", "cross-sample-event-key-v1",
    "cross-sample-candidate-key-v1", "clonal-evidence-id-v1", "trajectory-evidence-id-v1"))
  expect_setequal(getNamespaceExports("Neoantimon2"), current_public_api)
})

test_that("nine explicit context fixtures and five trajectory release fixtures validate", {
  cases<-jsonlite::read_json(longitudinal_fixture("context-cases.json"))
  expect_length(cases,9L)
  for(s in cases) expect_silent(Neoantimon2:::longitudinal_settings(longitudinal_project(s)))
  for(c in jsonlite::read_json(longitudinal_fixture("trajectory-cases.json")))
    expect_identical(longitudinal_cells(unlist(c$states),ordered=c$ordered)$trajectory_state,c$expected)
})

test_that("explicit subject assignments and relations enforce subject boundaries", {
  s <- longitudinal_settings_fixture(); p <- longitudinal_project(s)
  expect_identical(Neoantimon2:::longitudinal_settings(p)$subjects[[1L]]$subject_id, "patient-1")
  s$contexts[[1L]]$subject_id <- NULL
  expect_error(Neoantimon2:::longitudinal_settings(longitudinal_project(s)), "subject_id|nonempty")
  s <- longitudinal_settings_fixture(); s$subjects <- c(s$subjects, list(list(subject_id = "other")))
  s$contexts[[3L]]$subject_id <- "other"
  expect_error(Neoantimon2:::longitudinal_settings(longitudinal_project(s)), "Cross-subject")
})

test_that("time labels do not define chronology and ISO/relative/ordinal axes are explicit", {
  s <- list(sample_id = "S1", timepoint = "T1")
  expect_null(Neoantimon2:::longitudinal_time(s)$value)
  for (label in c("baseline", "T2", "followup")) {
    s$timepoint <- label; expect_null(Neoantimon2:::longitudinal_time(s)$value)
  }
  expect_equal(Neoantimon2:::longitudinal_time(list(relative_time = 3, time_unit = "days", time_origin = "surgery"))$value, 3)
  expect_equal(Neoantimon2:::longitudinal_time(list(ordinal_time = 2))$value, 2)
  expect_equal(Neoantimon2:::longitudinal_time(list(collection_time = "2026-10-01T00:00:00Z"))$axis, "absolute_utc")
  expect_error(Neoantimon2:::longitudinal_time(list(collection_time = "2026-02-30")), "Invalid")
  expect_error(Neoantimon2:::longitudinal_time(list(collection_time = "2026-10-01T00:00:00")), "UTC")
  expect_identical(Neoantimon2:::longitudinal_order(list(list(sample_id="T2"),list(sample_id="T1")))$ordering, "unordered")
  expect_identical(Neoantimon2:::longitudinal_order(list(list(sample_id="B",ordinal_time=0),list(sample_id="A",ordinal_time=1)))$contexts[[1]]$sample_id, "B")
})

test_that("matched normals and technical replicates retain explicit context relations", {
  s <- longitudinal_settings_fixture()
  s$contexts[[2]]$role <- "normal"; s$relations <- list(list(from="S1",to="S2",type="matched_normal"))
  expect_identical(Neoantimon2:::longitudinal_settings(longitudinal_project(s))$relations[[1]]$type, "matched_normal")
  s$relations <- list(list(from="S1",to="S2",type="technical_replicate"))
  s$contexts[[2]]$ordinal_time <- 2
  a <- list(sample_id="S1",state="observed"); b <- list(sample_id="S2",state="not_observed_evaluable")
  expect_null(Neoantimon2:::longitudinal_transition(a,b,s$contexts,s$relations))
  s$relations<-list(list(from="S1",to="S2",type="technical_replicate"),list(from="S2",to="S3",type="technical_replicate"))
  b$sample_id<-"S3"
  expect_null(Neoantimon2:::longitudinal_transition(a,b,s$contexts,s$relations))
})

test_that("conflicting temporal axes are rejected and single axes or unordered remain valid", {
  axes <- list(collection_time = "2026-10-01", relative_time = 3, ordinal_time = 2)
  combinations <- c(combn(names(axes), 2L, simplify = FALSE), list(names(axes)))
  for (fields in combinations) {
    s <- longitudinal_settings_fixture()
    s$contexts[[1L]][names(axes)] <- NULL
    s$contexts[[1L]][fields] <- axes[fields]
    s$contexts[[1L]]$time_unit <- "days"; s$contexts[[1L]]$time_origin <- "collection"
    expect_error(Neoantimon2:::longitudinal_settings(longitudinal_project(s)), "Conflicting temporal axes: at most one")
    expect_error(Neoantimon2:::longitudinal_time(s$contexts[[1L]]), "Conflicting temporal axes: at most one")
  }
  for (fields in c(lapply(names(axes), identity), list(character()))) {
    s <- longitudinal_settings_fixture()
    s$contexts[[1L]][names(axes)] <- NULL
    s$contexts[[1L]][fields] <- axes[fields]
    s$contexts[[1L]]$time_unit <- "days"; s$contexts[[1L]]$time_origin <- "collection"
    expect_silent(Neoantimon2:::longitudinal_settings(longitudinal_project(s)))
  }
})

test_that("small variant equivalence uses normalized allele and build not consequence", {
  key <- Neoantimon2:::longitudinal_event_key
  expect_identical(key(longitudinal_snv("S1")),key(longitudinal_snv("S2")))
  expect_false(identical(key(longitudinal_snv()),key(longitudinal_snv(pos=13L))))
  expect_identical(key(longitudinal_snv("S1",ref="AA",alt="A")),key(longitudinal_snv("S2",ref="AA",alt="A")))
  e <- longitudinal_snv(); e$genome_assembly <- "GRCh37"
  expect_false(identical(key(e),key(longitudinal_snv())))
  e$canonical_event <- NULL; expect_null(key(e))
})

test_that("fusion and splice keys require canonical exact structure", {
  key <- Neoantimon2:::longitudinal_event_key
  s <- list(sample_id="S1",reference_build="GRCh38",source_type="fusion",left=list(contig="1",position=4L,strand="+"),
    right=list(contig="2",position=8L,strand="-"),genes=list(list(id="g1",name="G1"),list(id="g2",name="G2")))
  e <- list(source_type="fusion",structure=s); f <- e; f$structure$sample_id <- "S2"
  expect_identical(key(e),key(f)); f$structure$right$position <- 9L; expect_false(identical(key(e),key(f)))
  e <- list(source_type="novel_junction",structure=list(sample_id="S1",reference_build="GRCh38",source_type="novel_junction",
    gene_id="g1",contig="1",strand="+",donor=10L,acceptor=20L,transcript_id="t1"))
  f <- e; f$structure$sample_id <- "S2"; expect_identical(key(e),key(f))
  f$structure$acceptor <- 21L; expect_false(identical(key(e),key(f)))
})

test_that("long-read equivalence references frozen isoform structure identity", {
  s <- Neoantimon2:::long_read_structure("GRCh38","1","+",data.frame(start=c(1L,20L),end=c(10L,30L)))
  id <- Neoantimon2:::long_read_structure_id(s)
  expect_identical(Neoantimon2:::longitudinal_event_key(list(isoform_structure_id=id,sample_id="S1")),id)
  expect_identical(Neoantimon2:::longitudinal_event_key(list(isoform_structure_id=id,sample_id="S2")),id)
})

test_that("candidate equivalence reuses frozen ligand and canonical presenter identity", {
  key <- function(peptide,hla) Neoantimon2:::biological_id("candidate_key",list(species="NCBITaxon:9606",hla_class="I",
    ligand_sequence=peptide,mhc_key=hla),Neoantimon2:::peptide_contract())
  a <- list(candidate_key=key("GILGFVFTL","HLA-A*02:01"),sample_id="S1",source_type="somatic",origin_ids=list("o1","o2"))
  b <- a; b$sample_id <- "S2"; b$source_type <- "fusion"
  expect_identical(Neoantimon2:::longitudinal_candidate_key(a),Neoantimon2:::longitudinal_candidate_key(b))
  b$candidate_key <- key("GILGFVFTL","HLA-B*07:02"); expect_false(identical(a$candidate_key,b$candidate_key))
  b$candidate_key <- key("GILGFVFTI","HLA-A*02:01"); expect_false(identical(a$candidate_key,b$candidate_key))
  expect_identical(a$origin_ids,list("o1","o2"))
})

test_that("compound origin equivalence preserves full biological event requirements", {
  links <- lapply(c(4L,13L),function(pos) {e<-longitudinal_snv(pos=pos);list(sample_id="S1",event_id=e$event_id,feature_key=Neoantimon2:::longitudinal_event_key(e))})
  o <- list(sample_id="S1",source_event_ids=as.list(vapply(links,`[[`,"","event_id")),transcript_id="t1")
  full <- Neoantimon2:::longitudinal_origin_basis(o,links)
  expect_length(full$required_event_keys,2L)
  o$source_event_ids <- o$source_event_ids[1L]
  expect_false(identical(full,Neoantimon2:::longitudinal_origin_basis(o,links)))
  o$source_event_ids <- list("unmapped"); expect_identical(Neoantimon2:::longitudinal_origin_basis(o,links)$status,"unresolved")
})

test_that("evaluability assertions require domain opportunity and explicit evidence", {
  s <- longitudinal_settings_fixture(); s$evaluability <- list(longitudinal_assay("S2","key"))
  expect_identical(Neoantimon2:::longitudinal_cell(s,"S2","candidate","key")$state,"not_observed_evaluable")
  for (domain in c("rna","ms","repertoire","event")) expect_identical(Neoantimon2:::longitudinal_cell(s,"S2",domain,"key")$state,"not_evaluated")
  s$evaluability[[1]]$adequate <- FALSE
  expect_error(Neoantimon2:::longitudinal_settings(longitudinal_project(s)),"adequate")
  s$evaluability[[1]]$adequate <- TRUE; s$evaluability[[1]]$evidence_ref <- NULL
  expect_error(Neoantimon2:::longitudinal_settings(longitudinal_project(s)),"evidence_ref")
})

test_that("all six typed states survive and contradictory observations remain ambiguous", {
  states <- c(evaluated_positive="observed",evaluated_negative="not_observed_evaluable",insufficient_evidence="not_observed_insufficient",
    not_evaluated="not_evaluated",not_applicable="not_applicable",ambiguous="ambiguous")
  for (state in names(states)) {
    s<-list(evaluability=list(longitudinal_assay("S2","k",state=state)))
    expect_identical(Neoantimon2:::longitudinal_cell(s,"S2","candidate","k")$state,unname(states[[state]]))
  }
  s<-list(evaluability=list(longitudinal_assay("S2","k")))
  expect_identical(Neoantimon2:::longitudinal_cell(s,"S2","candidate","k",list(list(candidate_id="c")))$state,"ambiguous")
})

test_that("gain and loss require both opportunity and chronology", {
  expect_identical(longitudinal_cells(c("observed","not_observed_evaluable"))$trajectory_state,"loss_with_supported_absence")
  expect_identical(longitudinal_cells(c("not_observed_evaluable","observed"))$trajectory_state,"gain_with_prior_supported_absence")
  for (state in c("not_evaluated","not_observed_insufficient","ambiguous"))
    expect_identical(longitudinal_cells(c("observed",state))$trajectory_state,"insufficient_to_classify")
  expect_length(longitudinal_cells(c("observed","not_observed_evaluable"),ordered=FALSE)$transitions,0L)
  expect_identical(longitudinal_cells(c("observed","not_observed_evaluable"),ordered=FALSE)$spatial_state,"region_specific_observation")
})

test_that("private/shared terminology does not convert incomplete comparison into private", {
  x<-longitudinal_cells(c("observed","not_evaluated"),ordered=FALSE)
  expect_identical(x$distribution,"observed_in_one_sample_only")
  expect_identical(x$trajectory_state,"insufficient_to_classify")
  expect_identical(x$spatial_state,"insufficient_evaluability")
  y<-longitudinal_cells(c("observed","not_observed_evaluable","not_observed_evaluable"),ordered=FALSE)
  expect_identical(y$spatial_state,"region_specific_observation")
  expect_identical(longitudinal_cells(c("observed","observed"),ordered=FALSE)$spatial_state,"observed_all_regions")
})

test_that("HLA genotype conflicts are reported without replacing sample declarations", {
  s<-list(list(subject_id="patient-1"))
  ct<-list(list(subject_id="patient-1",hla_genotype=list("A*02:01")),list(subject_id="patient-1",hla_genotype=list("HLA-A*02:01")))
  expect_identical(Neoantimon2:::longitudinal_hla_consistency(ct,s)[[1]]$status,"hla_genotype_consistent")
  ct[[2]]$hla_genotype<-list("HLA-B*07:02")
  expect_identical(Neoantimon2:::longitudinal_hla_consistency(ct,s)[[1]]$status,"hla_genotype_conflict")
  expect_identical(ct[[2]]$hla_genotype,list("HLA-B*07:02"))
  expect_identical(Neoantimon2:::longitudinal_hla_consistency(ct[1],s)[[1]]$status,"not_evaluated")
})

test_that("PyClone-VI mapping is explicit with ambiguous/unmapped states and run scope", {
  s<-longitudinal_settings_fixture()
  s$providers<-list(list(input_id="pc",provider="PyClone-VI",provider_version="0.1.6",run_id="run-1",subject_id="patient-1",mapping_input_id="map",options=list(seed=1)))
  extra<-list(list(input_id="pc",input_role="pyclone_vi",sample_id="S1",path=longitudinal_fixture("pyclone-vi.tsv"),version="1"),
    list(input_id="map",input_role="clonal_mapping",sample_id="S1",path=longitudinal_fixture("mapping.tsv"),version="1"))
  p<-longitudinal_project(s,extra)
  links<-list(list(sample_id="S1",subject_id="patient-1",event_id="e1"),list(sample_id="S2",subject_id="patient-1",event_id="e2"),list(sample_id="S1",subject_id="patient-1",event_id="e3"))
  a<-Neoantimon2:::longitudinal_clonal_import(p,s,links)
  expect_setequal(vapply(a,`[[`,"","mapping_status"),c("mapped","ambiguous","unmapped"))
  expect_equal(Filter(function(e)e$provider_mutation_id=="m1"&&e$sample_id=="S1",a)[[1]]$cluster_assignment_prob,.65)
  s$providers[[1]]$run_id<-"run-2"; b<-Neoantimon2:::longitudinal_clonal_import(p,s,links)
  expect_false(identical(vapply(a,`[[`,"","clonal_evidence_id"),vapply(b,`[[`,"","clonal_evidence_id")))
  expect_identical(a[[1]]$cluster_scope$run_id,"run-1")
  expect_false(identical(a[[1]]$cluster_scope,b[[1]]$cluster_scope))
  s$providers[[1]]$mapping_input_id<-NULL
  expect_true(all(vapply(Neoantimon2:::longitudinal_clonal_import(p,s,links),function(e)e$mapping_status=="unmapped",TRUE)))
})

test_that("provider row order does not change clonal identity", {
  s<-longitudinal_settings_fixture(); s$providers<-list(list(input_id="pc",provider="PyClone-VI",provider_version="1",run_id="run",subject_id="patient-1"))
  path<-tempfile(fileext=".tsv"); rows<-readLines(longitudinal_fixture("pyclone-vi.tsv"));writeLines(rows,path)
  extra<-list(list(input_id="pc",input_role="pyclone_vi",sample_id="S1",path=path,version="1"))
  a<-Neoantimon2:::longitudinal_clonal_import(longitudinal_project(s,extra),s,list())
  writeLines(c(rows[1],rev(rows[-1])),path)
  b<-Neoantimon2:::longitudinal_clonal_import(longitudinal_project(s,extra),s,list())
  expect_setequal(vapply(a,`[[`,"","clonal_evidence_id"),vapply(b,`[[`,"","clonal_evidence_id"))
})

test_that("provider file row order does not change trajectory evidence identities", {
  x<-longitudinal_e2e();extra<-x$extra
  lines<-readLines(extra[[1]]$path);path<-tempfile(fileext=".tsv")
  writeLines(c(lines[[1]],rev(lines[-1L])),path);extra[[1]]$path<-path
  o<-Neoantimon2:::longitudinal_read(neo_run(longitudinal_project(x$settings,extra)))
  expect_setequal(vapply(o$candidate_trajectories,`[[`,"","trajectory_evidence_id"),
    vapply(x$overlay$candidate_trajectories,`[[`,"","trajectory_evidence_id"))
  expect_setequal(vapply(o$event_trajectories,`[[`,"","trajectory_evidence_id"),
    vapply(x$overlay$event_trajectories,`[[`,"","trajectory_evidence_id"))
})

test_that("clonal mapping row order and path preserve all longitudinal semantics", {
  a <- longitudinal_origin_fixture()
  b <- longitudinal_origin_fixture("mapping-identity-reversed.tsv")
  expect_identical(longitudinal_without_physical(a), longitudinal_without_physical(b))
  expect_identical(vapply(a$event_links, `[[`, "", "feature_key"), vapply(b$event_links, `[[`, "", "feature_key"))
  expect_identical(vapply(a$clonal_evidence, `[[`, "", "clonal_evidence_id"), vapply(b$clonal_evidence, `[[`, "", "clonal_evidence_id"))
  for (field in c("event_trajectories", "candidate_trajectories", "receptor_trajectories")) {
    expect_identical(vapply(a[[field]], `[[`, "", "trajectory_evidence_id"), vapply(b[[field]], `[[`, "", "trajectory_evidence_id"))
    expect_identical(lapply(a[[field]], `[[`, "summary"), lapply(b[[field]], `[[`, "summary"))
  }
  expect_false(identical(a$clonal_evidence[[1L]]$mapping_artifact$sha256, b$clonal_evidence[[1L]]$mapping_artifact$sha256))
  expect_false(identical(a$clonal_evidence[[1L]]$mapping_artifact$path, b$clonal_evidence[[1L]]$mapping_artifact$path))
  expect_identical(a$clonal_evidence[[1L]]$mapping_semantics$role, "clonal_mapping")
  expect_identical(a$clonal_evidence[[1L]]$mapping_semantics$version, "mapping-v1")
  record <- a$clonal_evidence[[1L]]
  changed <- record; changed$mapping_semantics$rows[[1L]]$event_id <- "different-event"
  expect_false(identical(Neoantimon2:::longitudinal_id("clonal", Neoantimon2:::longitudinal_clonal_semantics(record)),
    Neoantimon2:::longitudinal_id("clonal", Neoantimon2:::longitudinal_clonal_semantics(changed))))
  expect_silent(Neoantimon2:::longitudinal_validate(a))
  expect_silent(Neoantimon2:::longitudinal_validate(b))
})

test_that("grouped candidate cells isolate clonal evidence to sample-local origins", {
  x <- longitudinal_origin_fixture()
  expect_length(x$candidate_trajectories, 1L)
  cells <- x$candidate_trajectories[[1L]]$cells
  for (i in seq_len(2L)) {
    cell <- Filter(function(c) c$sample_id == paste0("S", i), cells)[[1L]]
    expected <- paste0(c("X", "Y")[[i]], "_S", i)
    expect_length(cell$clonal_evidence, 1L)
    expect_identical(cell$clonal_evidence[[1L]]$event_ids, list(expected))
    expect_identical(cell$source_event_states[[1L]]$observations[[1L]]$event_id, expected)
    expect_length(cell$source_event_states, 1L)
    expect_length(Filter(function(e) e$sample_id == cell$sample_id, x$clonal_evidence), 2L)
  }
  missing <- Filter(function(c) c$sample_id == "S3", cells)[[1L]]
  expect_length(missing$clonal_evidence, 0L)
  expect_length(missing$source_event_states, 0L)
  expect_false(missing$candidate_generated)
})

test_that("real three-sample E2E retains A/B/C/D distributions and typed matrices", {
  x<-longitudinal_e2e(); o<-x$overlay
  expect_true(all(lengths(x$groups)>0))
  find<-function(k) Filter(function(t)t$feature_key==k,o$candidate_trajectories)[[1]]
  expect_identical(find(x$groups$A[[1]])$summary$trajectory_state,"observed_all_evaluable_samples")
  expect_identical(find(x$groups$B[[1]])$summary$trajectory_state,"loss_with_supported_absence")
  expect_identical(find(x$groups$C[[1]])$summary$trajectory_state,"gain_with_prior_supported_absence")
  expect_identical(find(x$groups$D[[1]])$summary$trajectory_state,"insufficient_to_classify")
  expect_length(find(x$groups$A[[1]])$summary$transitions,2L)
  expect_true(length(o$sample_matrix)>0)
  expect_identical(o$external_clonal_runs,0L)
  expect_true(length(o$event_trajectories)>0)
  expect_true(all(vapply(o$event_links,function(l)l$subject_id=="patient-1",TRUE)))
})

test_that("sample evidence and grade/rank remain frozen across integration", {
  x<-longitudinal_e2e(); rows<-neo_results(x$project)
  for(i in seq_along(x$runs)) {
    expect_identical(neo_results(x$runs[[i]]),x$before[[i]])
    grouped<-rows[rows$sample_id==paste0("S",i),,drop=FALSE]
    expect_equal(recognition_invariance_fields(grouped),recognition_invariance_fields(x$before[[i]]),ignore_attr=TRUE)
    expect_identical(grouped$candidate_id,x$before[[i]]$candidate_id)
  }
  expect_true(all(c("subject_id","cross_sample_candidate_key","trajectory_state","longitudinal_escape_status")%in%names(rows)))
  empty<-neo_results(x$project,sample_id="not-a-sample")
  expect_equal(nrow(empty),0L)
  expect_identical(names(empty),names(rows))
  expect_true(all(vapply(x$overlay$candidate_links,function(l)l$sample_result$sample_id==l$sample_id,TRUE)))
})

test_that("event CCF and imported cellular prevalence both survive with uncertainty", {
  x<-longitudinal_e2e(); e<-x$overlay$clonal_evidence[[1]]
  l<-Filter(function(l)l$sample_id==e$sample_id&&any(unlist(l$sample_result$source_event_ids)%in%unlist(e$event_ids)),x$overlay$candidate_links)[[1]]
  expect_false(identical(l$sample_result$ccf_value,e$cellular_prevalence))
  expect_true(is.numeric(l$sample_result$ccf_value))
  expect_equal(e$cellular_prevalence_std,.08);expect_equal(e$cluster_assignment_prob,.65)
  p<-x$runs[[match(e$sample_id,c("S1","S2","S3"))]]
  expect_identical(l$evidence$priority$ccf_solutions,neo_explain(p,candidate_id=l$candidate_id)$candidates[[1]]$priority$ccf_solutions)
})

test_that("HLA LOH/APM trajectories retain class specificity and no causality", {
  x<-longitudinal_e2e();t<-Filter(function(t)t$feature_key%in%x$groups$A&&any(vapply(t$escape$hla_loh_transitions,function(z)
    identical(z$description,"HLA loss newly supported in later sample"),TRUE)),x$overlay$candidate_trajectories)
  expect_true(length(t)>0)
  expect_true(any(vapply(t,function(t)length(t$escape$apm_transitions)>0,TRUE)))
  for(l in Filter(function(l)l$sample_id=="S3",x$overlay$candidate_links)) {
    status<-l$evidence$presentation_evidence$summary$apm_evidence_status
    if(l$hla_class=="I") expect_identical(status,"direct_loss_supported") else expect_identical(status,"not_evaluated")
  }
  text<-neo_explain(x$project,candidate_id=x$overlay$candidate_links[[1]]$candidate_id,format="text")
  for(label in c("BIOLOGICAL CANDIDATE IDENTITY","SAMPLE OBSERVATIONS","SOURCE EVENT TRAJECTORY","CLONAL / CELLULAR","RNA / EXPRESSION","TREATMENT/TIME","EVALUABILITY / UNCERTAINTY")) expect_match(text,label,fixed=TRUE)
  expect_false(grepl("selected by treatment|caused by therapy|therapy induced",text))
})

test_that("MS evidence stays in its sample and receptor sampling never implies disappearance", {
  x<-longitudinal_e2e()
  for(l in x$overlay$candidate_links) {
    status<-l$evidence$presentation_evidence$summary$immunopeptidomics_status
    if(l$sample_id=="S1") expect_identical(status,"exact_unmodified_ligand_observed") else expect_identical(status,"no_evidence_supplied")
  }
  expect_true(length(x$overlay$receptor_trajectories)>0)
  for(t in x$overlay$receptor_trajectories) {
    expect_identical(t$summary$trajectory_state,"sample_specific_receptor_observations")
    expect_length(t$summary$transitions,0L)
    expect_identical(Filter(function(c)c$sample_id=="S3",t$cells)[[1]]$state,"not_evaluated")
  }
})

test_that("input sample order does not alter identity or trajectory output", {
  x<-longitudinal_e2e();s<-x$settings
  for(k in c("contexts","subjects","relations","evaluability","sample_runs")) s[[k]]<-rev(s[[k]])
  p<-neo_run(longitudinal_project(s,x$extra));o<-Neoantimon2:::longitudinal_read(p)
  expect_identical(o$candidate_trajectories,x$overlay$candidate_trajectories)
  expect_identical(o$event_trajectories,x$overlay$event_trajectories)
})

test_that("adding a sample leaves prior biological keys unchanged", {
  x<-longitudinal_e2e();s<-x$settings
  s$contexts<-s$contexts[1:2];s$sample_runs<-s$sample_runs[1:2]
  s$relations<-Filter(function(r)all(c(r$from,r$to)%in%c("S1","S2")),s$relations)
  s$treatments<-list();s$evaluability<-Filter(function(e)e$sample_id%in%c("S1","S2"),s$evaluability)
  s$providers<-list()
  o<-Neoantimon2:::longitudinal_read(neo_run(longitudinal_project(s)))
  expect_true(all(vapply(o$candidate_links,`[[`,"","feature_key")%in%vapply(x$overlay$candidate_links,`[[`,"","feature_key")))
  expect_true(all(vapply(o$event_links,`[[`,"","feature_key")%in%vapply(x$overlay$event_links,`[[`,"","feature_key")))
})

test_that("single-sample results preserve exact ordering and frozen IDs and grades", {
  p<-neo_run(recognition_gold_project());before<-neo_results(p)
  s<-list(version=Neoantimon2:::longitudinal_version(),subjects=list(list(subject_id="patient-1")),
    contexts=list(list(sample_id="sample-1",subject_id="patient-1")))
  f<-tempfile(fileext=".json");jsonlite::write_json(s,f,auto_unbox=TRUE)
  after<-neo_run(recognition_clone_with(p,list(list(input_id="ls",input_role="longitudinal_settings",sample_id="sample-1",path=f,version="1"))))
  expect_equal(recognition_invariance_fields(neo_results(after)),recognition_invariance_fields(before))
  expect_identical(neo_results(after)$candidate_id,before$candidate_id)
  expect_identical(names(neo_results(after)),names(before))
})

test_that("single-sample neo_multisample preserves exact legacy behavior with longitudinal settings", {
  p <- neo_run(recognition_gold_project())
  before <- neo_multisample(p)
  s <- list(version = Neoantimon2:::longitudinal_version(), subjects = list(list(subject_id = "patient-1")),
    contexts = list(list(sample_id = "sample-1", subject_id = "patient-1", ordinal_time = 1)))
  f <- tempfile(fileext = ".json"); jsonlite::write_json(s, f, auto_unbox = TRUE)
  after <- neo_run(recognition_clone_with(p, list(list(input_id = "ls", input_role = "longitudinal_settings",
    sample_id = "sample-1", path = f, version = "1"))))
  expect_identical(neo_multisample(after), before)
  for (patient in list(NULL, "patient-1", "missing")) for (peptide in list(NULL, "GILGFVFTL"))
    expect_identical(neo_multisample(after, patient_id = patient, peptide = peptide),
      neo_multisample(p, patient_id = patient, peptide = peptide))
  expect_error(neo_multisample(after, patient_id = 1), "patient_id")
  # A matched normal does not make this a multi-tumor project.
  overlay <- Neoantimon2:::longitudinal_read(after)
  overlay$context_graph$contexts <- c(overlay$context_graph$contexts,
    list(list(sample_id = "normal-1", subject_id = "patient-1", role = "normal")))
  expect_identical(testthat::with_mocked_bindings(neo_multisample(after),
    longitudinal_read = function(project) overlay, .package = "Neoantimon2"), before)
  expect_identical(testthat::with_mocked_bindings(neo_results(after),
    longitudinal_read = function(project) overlay, .package = "Neoantimon2"), neo_results(after))
  expect_identical(testthat::with_mocked_bindings(neo_explain(after, candidate_id = neo_results(after)$candidate_id[[1L]]),
    longitudinal_read = function(project) overlay, .package = "Neoantimon2"),
    neo_explain(after, candidate_id = neo_results(after)$candidate_id[[1L]]))
})

test_that("source NOT_GRADED policies survive longitudinal summaries", {
  p<-presentation_source_pair()$after;rows<-neo_results(p)
  s<-list(version=Neoantimon2:::longitudinal_version(),subjects=list(list(subject_id="patient-1")),contexts=list(list(sample_id="sample-1",subject_id="patient-1")))
  f<-tempfile(fileext=".json");jsonlite::write_json(s,f,auto_unbox=TRUE)
  after<-neo_run(recognition_clone_with(p,list(list(input_id="ls",input_role="longitudinal_settings",sample_id="sample-1",path=f,version="1"))))
  expect_true(all(neo_results(after)$final_grade=="NOT_GRADED"))
  expect_identical(neo_results(after)$candidate_id,rows$candidate_id)
  expect_identical(neo_results(after)$grade_reason,rows$grade_reason)
  o<-Neoantimon2:::longitudinal_read(after)
  expect_true(any(vapply(o$candidate_links,function(l) any(vapply(l$origin_equivalence,
    function(origin) length(origin$required_event_keys)>0L,TRUE)),TRUE)))
  expect_true(all(vapply(o$candidate_links,function(l) all(vapply(l$origin_equivalence,
    function(origin) origin$status=="mapped",TRUE)),TRUE)))
})

test_that("longitudinal sidecar tampering and changed settings fail closed", {
  x<-longitudinal_e2e();s<-x$settings;p<-neo_run(longitudinal_project(s,x$extra))
  root<-Neoantimon2:::longitudinal_path(p);file<-list.files(root,pattern="longitudinal.rds",recursive=TRUE,full.names=TRUE)[[1]]
  writeLines("tampered",file)
  expect_error(neo_results(p),"integrity")
  p<-longitudinal_project(s,x$extra);writeLines("{}",p$inputs[[1]]$resolved_path)
  expect_error(neo_run(p),"Unsupported|changed")
})

test_that("derived overlay validator rejects cross-subject evidence and unsupported absence", {
  x<-longitudinal_e2e()$overlay
  expect_silent(Neoantimon2:::longitudinal_validate(x))
  bad<-x;bad$candidate_links[[1]]$subject_id<-"other"
  expect_error(Neoantimon2:::longitudinal_validate(bad),"Cross-subject")
  bad<-x;bad$candidate_trajectories[[1]]$cells[[1]]$state<-"not_observed_evaluable"
  bad$candidate_trajectories[[1]]$cells[[1]]$evaluability<-NULL
  expect_error(Neoantimon2:::longitudinal_validate(bad),"Unsupported absence")
  bad<-x
  observed_cell<-which(vapply(bad$candidate_trajectories[[1]]$cells,
    function(cell) length(cell$observations)>0L,TRUE))[[1L]]
  bad$candidate_trajectories[[1]]$cells[[observed_cell]]$observations[[1]]$subject_id<-"other"
  expect_error(Neoantimon2:::longitudinal_validate(bad),"bleeding")
})
