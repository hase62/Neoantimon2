expanded_setup <- function(mode="good",profile="dr-expanded-v1") {
  root<-tempfile("expanded ; literal ");dir.create(root)
  loci<-ie_loci();da<-stats::setNames(lapply(loci,function(l)paste0(l,"*",if(l%in%c("MICA","MICB"))"002:01:01"else"02:01")),loci)
  re<-stats::setNames(lapply(loci,function(l)paste0(l,"*",if(l%in%c("MICA","MICB"))"001:01:01"else"01:01")),loci)
  re$DRA<-"DRA*01:01";re$DQB1<-"DQB1*05:01";re$DPA1<-"DPA1*02:02";re$DPB1<-"DPB1*19:01"
  alleles<-unique(c(unlist(da),unlist(re)))
  base<-jsonlite::read_json(system.file("extdata","indirect-autorun","synthetic-reference.json",package="Neoantimon2"))
  template<-base$records[[1]]
  base$records<-lapply(alleles,function(a){r<-template;r$allele<-ie_allele(a);r$sequence_identifier<-paste0("SYNTHETIC:",a)
    r$protein_sequence<-if(a%in%unlist(da))"ACDEFGHIKLMNPQRSTVWY"else"ACDEFGHIKAAAAAA";n<-nchar(r$protein_sequence)
    r$domains<-r$covered_regions<-list(list(name="full_precursor",start=1L,end=n));r})
  ref<-file.path(root,"reference.json");ia_json(base,ref)
  tool<-file.path(root,"model");dir.create(file.path(tool,"data"),recursive=TRUE)
  file.copy(system.file("extdata","indirect-expanded","mock-class2.py",package="Neoantimon2"),file.path(tool,"netMHCIIpan"));Sys.chmod(file.path(tool,"netMHCIIpan"),"0755");writeLines(mode,file.path(tool,"data/mode.txt"))
  list(donor=da,recipient=re,profile=profile,reference=list(path=ref,sha256=ia_hash(ref)),installation=list(path=tool,authorization=list(status="mock_only",terms="Controlled software fixture only")),output_dir=file.path(root,"run"),
    source_loci=loci,self_loci=loci)
}
expanded_call <- function(a)do.call(neo_indirect_class2,a)

test_that("13 sources retain every Origin while exact queries and cores deduplicate",{
 a<-expanded_setup();x<-expanded_call(a)
 expect_identical(x$status,"COMPLETED");expect_false(x$manifest$real_predictor_executed)
 expect_equal(x$result$donor_origin_count,78L);expect_equal(x$result$self_origin_count,13L)
 expect_length(x$result$candidates,24L);expect_length(x$result$candidate_origin_links,312L)
 expect_true(all(vapply(x$plan$presenters,function(p)p$expression_confirmed&&p$expression_basis=="typing_based_assumption"&&p$experimental_expression_status=="not_verified_by_workflow",TRUE)))
 expect_true(all(vapply(x$result$presenter_summary,function(p)p$expression_confirmed&&p$expression_basis=="typing_based_assumption"&&p$experimental_expression_status=="not_verified_by_workflow",TRUE)))
 saved<-jsonlite::read_json(file.path(a$output_dir,"result.json"))
 expect_identical(saved$presenter_summary,x$result$presenter_summary)
 table<-utils::read.delim(file.path(a$output_dir,"presenter_summary.tsv"),check.names=FALSE)
 expect_true(all(table$expression_basis=="typing_based_assumption"&table$experimental_expression_status=="not_verified_by_workflow"))
 expect_equal(x$result$group_summary[[1]]$complete_count,20L)
 expect_setequal(unique(vapply(x$result$locus_summary,`[[`,"","locus")),ie_loci())
 expect_true(all(vapply(x$result$candidates,function(c)length(c$origin_ids)==13L,TRUE)))
 expect_equal(x$result$unique_query_count,28L)
 expect_true(all(vapply(x$result$candidates,function(c)c$prediction_status=="evaluated",TRUE)))
 expect_true(all(vapply(x$plan$presenters,function(p)p$status=="evaluated",TRUE)))
 expect_null(x$result$clinical_risk_score)
 before<-ia_files(a$output_dir);a$resume<-TRUE;y<-expanded_call(a);expect_identical(x$result,y$result);expect_identical(ia_files(a$output_dir),before)
 a$thresholds<-list(DR=500,DQ=1000,DP=1000);expect_error(expanded_call(a),"mismatch")
})

test_that("DQ DP retain only explicitly confirmed pairs and inversion",{
 a<-expanded_setup(profile="dqdp-exploratory-v1");expect_error(expanded_call(a),"explicit confirmed")
 a$presenters<-list(list(pair="DQA1*01:01~DQB1*05:01",pairing="confirmed_cis",expression_confirmed=TRUE),list(pair="DPA1*02:02~DPB1*19:01",pairing="explicit_expressed",expression_confirmed=TRUE))
 x<-expanded_call(a);expect_identical(x$status,"COMPLETED");expect_identical(x$result$schema_version,"indirect-class2-dqdp-exploratory-1.0.0")
 expect_true(jsonvalidate::json_validate(jsonlite::toJSON(x$result,auto_unbox=TRUE,null="null",digits=NA),system.file("schema",paste0(x$result$schema_version,".json"),package="Neoantimon2"),engine="ajv"))
 expect_length(x$plan$presenters,2L);expect_setequal(vapply(x$result$group_summary,`[[`,"","group"),c("DQ","DP"))
 pred<-jsonlite::read_json(file.path(a$output_dir,"predictions.json"))
 dp<-Filter(function(r)grepl("DPB1",r$presenter),pred)
 expect_true(all(vapply(dp,function(r)r$cores[[1]]$orientation=="reverse",TRUE)))
 expect_true(all(vapply(dp,function(r)r$raw_output$forward_ligand_core_start0==6L,TRUE)))
 expect_false(x$result$official_pirche_reproduction);expect_null(x$result$clinical_risk_score)
 expect_error(ie_vendor("DRA*01:01~DRB1*01:01:01"),"no truncation")
 expect_error(ie_presenters(a$presenters,ie_typing(a$recipient,names(a$recipient),"recipient"),"dr-expanded-v1"),"outside")
})

test_that("absence, missing typing, partial sequence and unsupported model remain different",{
 a<-expanded_setup();a$donor$DRB3<-list(status="absent",copy_number=0);x<-expanded_call(a)
 expect_identical(x$status,"COMPLETED");expect_equal(x$result$donor_origin_count,72L)
 a<-expanded_setup();a$donor$DRB3<-NULL;x<-expanded_call(a);expect_identical(x$status,"INCOMPLETE");expect_null(x$result$group_summary[[1]]$complete_count)
 a<-expanded_setup();a$recipient$MICA<-list(status="uncertain",alleles="MICA*001:01:01");x<-expanded_call(a)
 expect_false(x$result$self_complete);expect_null(x$result$presenter_summary[[1]]$complete_count)
 a<-expanded_setup();m<-jsonlite::read_json(a$reference$path);m$records[[1]]$sequence_status<-"partial";m$records[[1]]$covered_regions<-list();ia_json(m,a$reference$path);a$reference$sha256<-ia_hash(a$reference$path)
 x<-expanded_call(a);expect_identical(x$status,"INCOMPLETE");expect_false(x$result$source_complete)
 a<-expanded_setup();a$recipient$DRB1<-c("DRB1*01:01","DRB1*99:99");m<-jsonlite::read_json(a$reference$path);r<-m$records[[1]];r$allele<-"HLA-DRB1*99:99";m$records<-c(m$records,list(r));ia_json(m,a$reference$path);a$reference$sha256<-ia_hash(a$reference$path)
 x<-expanded_call(a);expect_identical(x$status,"INCOMPLETE");expect_true(any(vapply(x$plan$presenters,function(p)p$status=="unsupported",TRUE)))
 expect_true(any(vapply(x$result$candidates,function(c)c$prediction_status=="unsupported"&&is.null(c$affinity_nM),TRUE)))
})

test_that("BA boundaries, EL-only, tampering and input policies fail closed",{
 a<-expanded_setup("boundary");x<-expanded_call(a);expect_equal(x$result$group_summary[[1]]$complete_count,0L)
 a<-expanded_setup("el-only");expect_error(expanded_call(a),"EL-only");expect_identical(jsonlite::read_json(file.path(a$output_dir,"run.json"))$status,"FAILED")
 a<-expanded_setup("badcore");expect_error(expanded_call(a),"core does not map")
 a<-expanded_setup();a$dry_run<-TRUE;a$installation<-NULL;x<-expanded_call(a);expect_identical(x$status,"DRY_RUN");expect_false(x$manifest$real_predictor_executed)
 writeLines("tampered",file.path(a$output_dir,"queries.tsv"));a$resume<-TRUE;expect_error(expanded_call(a),"checksum")
 a<-expanded_setup();a$batch_size<-1025;expect_error(expanded_call(a),"1024")
 a<-expanded_setup();a$source_loci<-c("A","A");expect_error(expanded_call(a),"unique")
 expect_error(ie_typing(list(DRB3=list(status="absent",alleles="DRB3*01:01")),"DRB3","donor"),"Absent")
 expect_identical(ie_typing(list(),"DRB3","donor")[[1]]$status,"untyped")
 expect_identical(ie_allele(" mica*001:01:01 "),"HLA-MICA*001:01:01")
 expect_error(ie_allele("A*01:01","B"),"conflict")
})

test_that("IPD cached files are versioned, verified and never silently repaired",{
 root<-tempfile();dir.create(file.path(root,"3.66.0"),recursive=TRUE)
 c<-file.path(root,"3.66.0");files<-c("LICENCE.md","release_version.txt","Allele_status.txt","hla_prot.fasta")
 writeLines("IPD-IMGT/HLA 3.66.0",file.path(c,files[[2]]));writeLines("Synthetic licence",file.path(c,files[[1]]))
 writeLines(c("Allele,Start,End,Partial","A*01:01:01,1,63,Full","A*01:02:01,4,60,Partial"),file.path(c,files[[3]]))
 writeLines(c(">HLA:FAKE A*01:01:01 20 bp","ACDEFGHIKLMNPQRSTVWY",">HLA:FAKE2 A*01:02:01 19 bp","CDEFGHIKLMNPQRSTVWY"),file.path(c,files[[4]]))
 source<-list(release="3.66.0",commit=paste(rep("a",40),collapse=""),files=stats::setNames(lapply(files,function(f)list(sha256=ia_hash(file.path(c,f)))),files));ia_json(source,file.path(c,"source.json"))
 cache<-ie_ipd_cache(list(release="3.66",cache=root,download=FALSE));rr<-ie_ipd_records(cache,c("HLA-A*01:01","HLA-A*01:02:01","HLA-A*99:99"))
 expect_identical(rr[[1]]$sequence_status,"complete");expect_identical(rr[[2]]$sequence_status,"partial");expect_null(rr[[2]]$coordinate_map$precursor_start)
 expect_identical(rr[[3]]$sequence_status,"unavailable")
 writeLines("changed",file.path(c,"hla_prot.fasta"));expect_error(ie_ipd_cache(list(release="3.66",cache=root)),"integrity")
 expect_error(ie_ipd_cache(list(release="3.66",cache=tempfile(),download=FALSE)),"unavailable")
})


test_that("explicit homozygous copies retain inputs without duplicating synthesis paths",{
 a<-expanded_setup();a$donor$A<-rep(a$donor$A,2);a$recipient$DRA<-rep(a$recipient$DRA,2)
 x<-expanded_call(a);expect_identical(x$status,"COMPLETED");expect_equal(x$result$donor_origin_count,78L)
 input<-jsonlite::read_json(file.path(a$output_dir,"inputs.json"));expect_length(input$donor$A,2L)
 expect_true(all(vapply(Filter(function(c)c$group=="DR",x$result$candidates),function(c)length(c$origin_ids)==13L,TRUE)))
 f<-system.file("schema",paste0(x$result$schema_version,".json"),package="Neoantimon2")
 expect_true(jsonvalidate::json_validate(jsonlite::toJSON(x$result,auto_unbox=TRUE,null="null",digits=NA),f,engine="ajv"))
})


test_that("empty and incompletely characterized pools never become fabricated negatives",{
 a<-expanded_setup();a$source_loci<-a$self_loci<-"A";a$donor$A<-list(status="absent",copy_number=0)
 x<-expanded_call(a);expect_identical(x$status,"COMPLETED");expect_equal(x$result$donor_origin_count,0L);expect_equal(x$result$group_summary[[1]]$complete_count,0L)
 expect_true(all(vapply(x$result$locus_summary,function(l)l$status=="known_no_eligible_windows",TRUE)))
 a<-expanded_setup();a$source_loci<-a$self_loci<-"A";a$donor$A<-NULL;x<-expanded_call(a)
 expect_identical(x$status,"INCOMPLETE");expect_null(x$result$group_summary[[1]]$complete_count)
 expect_true(all(vapply(x$result$locus_summary,function(l)l$status=="incomplete_source"&&is.null(l$complete_count),TRUE)))
 a<-expanded_setup(profile="dqdp-exploratory-v1");a$presenters<-list(list(pair="DQA1*01:01~DQB1*05:01",pairing="confirmed_trans",expression_confirmed=FALSE))
 x<-expanded_call(a);expect_identical(x$status,"INCOMPLETE");expect_false(x$manifest$real_predictor_executed)
 expect_null(x$result$group_summary[[1]]$complete_count);expect_true(all(vapply(x$result$candidates,function(c)is.null(c$affinity_nM),TRUE)))
 expect_error(ie_presenters(list(list(pair="DQA1*01:01~DQB1*05:01",pairing="unknown",expression_confirmed=TRUE)),list(),"dqdp-exploratory-v1"),"explicitly established")
})


test_that("automatic DR coverage cannot hide an untyped additional beta locus",{
 a<-expanded_setup();a$source_loci<-a$self_loci<-"A";a$recipient$DRB3<-NULL
 x<-expanded_call(a);expect_identical(x$status,"INCOMPLETE");expect_false(x$result$presenter_scope_complete)
 expect_null(x$result$group_summary[[1]]$complete_count)
 expect_true(all(vapply(x$result$presenter_summary,function(s)!is.null(s$complete_count),TRUE)))
})


test_that("fresh runs preserve identities and every link resolves to the exact molecule",{
 a<-expanded_setup();x<-expanded_call(a);origins<-x$explanation$origins
 oi<-stats::setNames(origins,vapply(origins,`[[`,"","origin_id"));ci<-stats::setNames(x$result$candidates,vapply(x$result$candidates,`[[`,"","candidate_id"))
 expect_true(all(vapply(x$result$candidate_origin_links,function(l)l$origin_id%in%names(oi)&&l$candidate_id%in%names(ci),TRUE)))
 expect_true(all(vapply(x$result$candidate_origin_links,function(l)oi[[l$origin_id]]$role=="donor"&&oi[[l$origin_id]]$peptide_sequence==ci[[l$candidate_id]]$peptide_sequence&&l$presenter==ci[[l$candidate_id]]$presenter,TRUE)))
 a$output_dir<-paste0(a$output_dir,"-independent");y<-expanded_call(a);expect_identical(x$result,y$result)
})


test_that("valid three-digit protein-number fields are preserved without truncation",{
 expect_identical(ie_vendor("DRA*01:01~DRB1*04:316"),"DRB1_04316")
 expect_identical(ie_vendor("DQA1*01:123~DQB1*06:201"),"HLA-DQA101123-DQB106201")
 expect_error(ie_vendor("DRA*01:01~DRB1*04:316:01"),"no truncation")
})


test_that("official four-digit protein names stay exact and unsupported is not zero",{
 expect_identical(ie_allele("DPB1*1033:01"),"HLA-DPB1*1033:01")
 expect_identical(ie_vendor("DPA1*02:02~DPB1*1033:01"),"HLA-DPA10202-DPB1103301")
 expect_error(ie_pair("DQA1*01:01~DPB1*01:01"),"incompatible")
 a<-expanded_setup(profile="dqdp-exploratory-v1");a$recipient$DPB1<-"DPB1*1033:01"
 m<-jsonlite::read_json(a$reference$path);r<-m$records[[1]];r$allele<-"HLA-DPB1*1033:01";m$records<-c(m$records,list(r));ia_json(m,a$reference$path);a$reference$sha256<-ia_hash(a$reference$path)
 a$presenters<-list(list(pair="DPA1*02:02~DPB1*1033:01",pairing="explicit_expressed",expression_confirmed=TRUE))
 x<-expanded_call(a);expect_identical(x$status,"INCOMPLETE");expect_identical(x$plan$presenters[[1]]$status,"unsupported")
 expect_null(x$result$group_summary[[1]]$complete_count);expect_false(x$manifest$real_predictor_executed)
})


test_that("typing-based DR eligibility is not experimental expression confirmation",{
 typing<-ie_typing(list(DRA="DRA*01:01",DRB1="DRB1*01:01"),c("DRA","DRB1"),"recipient")
 automatic<-ie_presenters(NULL,typing,"dr-expanded-v1")[[1]]
 expect_identical(automatic$status,"eligible")
 expect_true(automatic$expression_confirmed) # Existing analysis admission stays intact.
 expect_identical(automatic$expression_basis,"typing_based_assumption")
 expect_identical(automatic$experimental_expression_status,"not_verified_by_workflow")
 declaration<-list(pair="DRA*01:01~DRB1*01:01",pairing="explicit_DR_single_alpha",
   expression_confirmed=TRUE,evidence="User-supplied assay claim; not verified here")
 explicit<-ie_presenters(list(declaration),typing,"dr-expanded-v1")[[1]]
 expect_identical(explicit$expression_basis,"user_declaration")
 expect_identical(explicit$experimental_expression_status,"not_verified_by_workflow")
 expect_identical(explicit$evidence,declaration$evidence)
 expect_identical(explicit[c("presenter","vendor","pairing","status","expression_confirmed")],automatic[c("presenter","vendor","pairing","status","expression_confirmed")])
 declaration$expression_confirmed<-FALSE
 unconfirmed<-ie_presenters(list(declaration),typing,"dr-expanded-v1")[[1]]
 expect_identical(unconfirmed$status,"expression_unconfirmed")
 expect_identical(unconfirmed$experimental_expression_status,"not_verified_by_workflow")
})
