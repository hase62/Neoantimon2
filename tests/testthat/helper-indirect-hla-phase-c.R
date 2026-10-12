phase_c_case <- local({
  cache <- new.env(parent = emptyenv())
  function(name = 'complete') {
    if (exists(name, cache, inherits = FALSE)) return(cache[[name]])
    input <- indirect_input()
    m <- indirect_fixture('binding.json')
    if (name == 'incomplete_self') input$a$self_search[[1L]]$catalog_complete <- FALSE
    if (name == 'incomplete_scope') {
      input <- indirect_second_presenter()
      input$b$presenters <- input$b$presenters[1L]
    }
    if (name == 'accessory') {
      p <- input$b$presenters[[1L]]; p$presenter <- 'HLA-DRA*01:01~HLA-DRB3*01:01'
      input$b$presenters <- c(input$b$presenters, list(p))
    }
    if (name == 'provider_failure') {
      m$rows[[2L]]$status <- 'failed'; m$rows[[2L]]['affinity_nM'] <- list(NULL)
      m$rows[[2L]]$cores <- list(); m$rows[[2L]]$reason_codes <- list('provider_execution_failed')
      input <- indirect_with_binding(m, input)
    }
    if (name == 'nonbinder') {
      m$rows[[2L]]$affinity_nM <- 1000
      input <- indirect_with_binding(m, input)
    }
    if (name == 'union') {
      m$rows[[1L]]$cores[[1L]]$offset <- 1L
      m$rows[[1L]]$cores[[1L]]$core_sequence <- 'CDEFGHIKL'
      input <- indirect_with_binding(m, input)
    }
    if (name == 'two_presenters') input <- indirect_second_presenter()
    if (name == 'missing_provider') {
      m$rows <- Filter(function(r) r$peptide_sequence != 'ACDEFGHIKAAAAAA', m$rows)
      input <- indirect_with_binding(m, input)
    }
    if (name == 'missing_universe') input$a$self_search <- list()
    if (name == 'missing_source') input$a$typing <- Filter(function(t) !(t$subject_id == 'donor-1' && t$locus == 'C'), input$a$typing)
    if (name == 'reference_failure') {
      r <- indirect_fixture('reference.json'); r$status <- 'failed'
      input <- indirect_with_reference(r, input)
    }
    if (name == 'core_ambiguous') {
      m$rows[[2L]]$cores[[1L]]$primary <- FALSE
      input <- indirect_with_binding(m, input)
    }
    x <- indirect_test(input)
    result <- Neoantimon2:::indirect_public_result(x$result, x$foundation)
    explain <- Neoantimon2:::indirect_public_explanation(x$result, x$foundation)
    cache[[name]] <- c(x, list(input = input, public_result = result, public_explanation = explain))
    cache[[name]]
  }
})
phase_c_runtime_project <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- neo_run(recognition_fixture_project(inputs = list()))
    cache
  }
})
phase_c_project <- function(name = 'complete') {
  x <- phase_c_case(name); p <- phase_c_runtime_project()
  p$config$allo <- Neoantimon2:::allo_config(x$input$a); p$allo_foundation <- x$foundation
  p$config$indirect_hla <- x$config; p$indirect_hla <- x$result
  p
}
phase_c_summary <- function(name = 'complete') phase_c_case(name)$public_result$contexts[[1L]]
phase_c_detail <- function(name = 'complete') phase_c_case(name)$public_explanation$contexts[[1L]]
phase_c_semantic <- function(x) { x$provenance <- NULL; x }
phase_c_schema <- function(name) {
  validator <- jsonvalidate::json_validator(system.file('schema', paste0('indirect-hla-', name, '-1.0.0.json'), package = 'Neoantimon2'), engine = 'ajv')
  function(x) as.logical(validator(jsonlite::toJSON(x, auto_unbox = TRUE, null = 'null', digits = NA)))
}
phase_c_json_fixture <- function(name) jsonlite::read_json(testthat::test_path('..', 'fixtures', 'indirect_hla_phase_c', name))
