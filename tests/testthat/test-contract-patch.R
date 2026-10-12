test_that("fixed independent identity vectors pin the revised 1.0.0 contract", {
  for (v in fixture("identity_vectors")) expect_identical(id(v$kind,v$values),v$expected)
  expect_false(any(c("somatic_candidate_id","allo_candidate_id") %in% names(contract$identities)))
  expect_identical(vapply(contract$identities$epitope_key$fields,`[[`,"","name"),
                   c("species","hla_class","binding_core","mhc_key"))
})

test_that("same sample and molecular candidate coexist across contexts", {
  x <- fixture("multi_context")
  expect_true(valid(x))
  a <- x$neoantigen_results[[1]]; b <- x$neoantigen_results[[2]]
  for (n in c("sample_id","candidate_key","peptide_key","mhc_key")) expect_identical(a[[n]],b[[n]])
  expect_false(identical(a$candidate_id,b$candidate_id))
  x$neoantigen_results[[2]]$candidate_id <- a$candidate_id
  expect_error(valid(x),"duplicate primary key")
})

test_that("all shared vocabularies resolve exclusively from the authority", {
  refs <- character()
  for (d in c(contract$tables,contract$records,contract$identities)) for (f in d$fields) {
    for (pair in list(c("vocabulary","enum"),c("items_vocabulary","items_enum"))) {
      ref <- f[[pair[1]]]
      if (!is.null(ref)) {
        refs <- c(refs,ref)
        expect_true(ref %in% names(contract$vocabularies))
        expect_null(f[[pair[2]]])
      }
      if (!is.null(f[[pair[2]]])) for (v in contract$vocabularies)
        expect_false(identical(f[[pair[2]]],v))
    }
  }
  expect_setequal(unique(refs),names(contract$vocabularies))
  for (ref in names(contract$vocabularies)) {
    cc <- contract; cc$vocabularies[[ref]] <- list("TEST_EXTENSION")
    f <- list(type="string",nullable=FALSE,vocabulary=ref)
    expect_true(Neoantimon2:::validate_value("TEST_EXTENSION",f,cc,"test"))
    expect_error(Neoantimon2:::validate_value("TEST_EXTENSION",f,contract,"test"),"enum")
    f <- list(type="array",nullable=FALSE,items_type="string",items_vocabulary=ref)
    expect_true(Neoantimon2:::validate_value(list("TEST_EXTENSION"),f,cc,"test"))
  }
  expect_error(Neoantimon2:::validate_value(list(),list(type="array",items_vocabulary="unknown"),contract,"test"),"Unknown")
})

test_that("required candidate identity and status reject null at both layers", {
  for (n in c("species","hla_class","prediction_status","grading_status","final_grade","grade_before_warnings")) {
    for (missing in list(NULL,NA_character_)) {
      x <- fixture("simple_somatic"); x$neoantigen_results[[1]][n] <- list(missing)
      expect_error(Neoantimon2:::validate_table(x$neoantigen_results,"neoantigen_results",contract),"null is not allowed")
      expect_error(valid(x),"null is not allowed")
    }
  }
  for (status in c("NOT_RUN","RUNNING","SUCCESS","PARTIAL","FAILED","UNSUPPORTED","NOT_APPLICABLE")) {
    x <- fixture("simple_somatic"); x$neoantigen_results[[1]]$prediction_status <- status
    expect_true(Neoantimon2:::validate_table(x$neoantigen_results,"neoantigen_results",contract))
    expect_true(valid(x))
  }
  for (grade in unlist(contract$vocabularies$grade)) {
    x <- fixture("simple_somatic"); r <- x$neoantigen_results[[1]]
    r$final_grade <- r$grade_before_warnings <- grade
    r$grading_status <- if (grade %in% c("INVALID","NOT_GRADED")) grade else "GRADED"
    r$prediction_status <- "SUCCESS"; x$neoantigen_results[[1]] <- r
    if (grade == "INVALID") x$peptide_origins[[1]]$origin_validity <- "INVALID"
    expect_true(Neoantimon2:::validate_table(x$neoantigen_results,"neoantigen_results",contract))
    expect_true(valid(x))
  }
})

test_that("canonical species rejects human aliases before hashing", {
  for (kind in c("candidate_key","epitope_key")) {
    v <- Filter(function(v) v$kind==kind,fixture("identity_vectors"))[[1]]$values
    for (alias in c("human","Homo sapiens","9606","NCBITaxon:09606")) {
      v$species <- alias
      expect_error(id(kind,v),"pattern")
    }
  }
})

test_that("different ligand coordinates and mapping retain one pMHC identity", {
  x <- fixture("class_ii_shared_core"); a <- x$neoantigen_results[[1]]; b <- x$neoantigen_results[[2]]
  expect_true(valid(x))
  expect_false(identical(a$ligand_sequence,b$ligand_sequence))
  expect_false(identical(a$binding_core_start,b$binding_core_start))
  expect_false(identical(a$epitope_assignments_json[[1]]$position_mapping,b$epitope_assignments_json[[1]]$position_mapping))
  expect_identical(a$epitope_key,b$epitope_key)
  v <- a[c("species","hla_class","binding_core","mhc_key")]
  v$binding_core_length <- 99L
  expect_error(id("epitope_key",v),"unknown fields")
  x$neoantigen_results[[1]]$epitope_assignments_json[[1]]$binding_core_length <- 99L
  expect_error(valid(x),"coordinates")
})

test_that("multisample summaries bind each context independently in mixed runs", {
  x <- fixture("multi_context")
  expect_true(valid(x)); expect_length(x$multisample_summary,2)
  x$multisample_summary[[2]]$analysis_context <- x$multisample_summary[[1]]$analysis_context
  expect_error(valid(x),"duplicate primary key")
  x <- fixture("multi_context"); x$multisample_summary[[1]]$analysis_context <- "missing"
  expect_error(valid(x),"foreign key")
  x <- fixture("multi_context"); x$multisample_summary[[1]]$analysis_context <- "ctx-donor-recipient"
  expect_error(valid(x),"Allo multisample")
  x <- fixture("multi_context"); x$multisample_summary[[1]]$patient_id <- "other"
  expect_error(valid(x),"patient/context")
})
