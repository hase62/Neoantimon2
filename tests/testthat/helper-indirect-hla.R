indirect_fixture <- function(file) jsonlite::read_json(testthat::test_path('..', 'fixtures', 'indirect_hla', file))
indirect_input <- function() {
  a <- indirect_fixture('foundation.json')
  a$reference$path <- normalizePath(testthat::test_path('..', 'fixtures', 'indirect_hla', 'reference.json'))
  b <- indirect_fixture('profile.json')
  b$provider$path <- normalizePath(testthat::test_path('..', 'fixtures', 'indirect_hla', 'binding.json'))
  list(a = a, b = b)
}
indirect_write <- function(x) {
  f <- tempfile(fileext='.json')
  jsonlite::write_json(x, f, pretty=TRUE, auto_unbox=TRUE, null='null', digits=NA)
  f
}
indirect_test <- function(input = indirect_input()) {
  f <- Neoantimon2:::allo_foundation(Neoantimon2:::allo_config(input$a))
  b <- Neoantimon2:::indirect_config(input$b)
  list(foundation=f, config=b, result=Neoantimon2:::indirect_baseline(f,b))
}
indirect_enable <- function(p, input = indirect_input()) {
  x <- indirect_test(input)
  p$config$allo <- Neoantimon2:::allo_config(input$a)
  p$allo_foundation <- x$foundation; p$config$indirect_hla <- x$config
  p$indirect_hla <- x$result; p
}
indirect_with_reference <- function(m, input = indirect_input()) {
  input$a$reference <- list(path=indirect_write(m)); input
}
indirect_with_binding <- function(m, input = indirect_input()) {
  input$b$provider <- list(path=indirect_write(m)); input
}
indirect_context_test <- function(input = indirect_input()) indirect_test(input)$result$contexts[[1L]]
indirect_second_presenter <- function() {
  x <- indirect_input(); m <- indirect_fixture('reference.json'); p <- indirect_fixture('binding.json')
  r <- Filter(function(r) r$allele == 'HLA-DRB1*01:01', m$records)[[1L]]
  r$allele <- 'HLA-DRB1*03:01'; r$sequence_identifier <- 'synthetic:DRB1*03:01'
  m$records <- c(m$records, list(r))
  t <- Filter(function(t) t$subject_id == 'patient-1' && t$locus == 'DRB1', x$a$typing)[[1L]]
  t$alleles <- list('DRB1*03:01'); x$a$typing <- c(x$a$typing, list(t))
  pair <- 'HLA-DRA*01:01~HLA-DRB1*03:01'
  pp <- x$b$presenters[[1L]]; pp$presenter <- pair
  x$b$presenters <- c(x$b$presenters, list(pp))
  p$provider$supported_presenters <- c(p$provider$supported_presenters,list(pair))
  p$rows <- c(p$rows,lapply(p$rows,function(r){r$presenter<-pair;r}))
  indirect_with_binding(p,indirect_with_reference(m,x))
}
indirect_accessory_typing <- function(locus, input = indirect_second_presenter()) {
  stopifnot(locus %in% c('DRB3','DRB4','DRB5'))
  m <- jsonlite::read_json(input$a$reference$path)
  r <- Filter(function(r) r$allele == 'HLA-DRB1*01:01', m$records)[[1L]]
  r$allele <- paste0('HLA-',locus,'*01:01')
  r$sequence_identifier <- paste0('synthetic:accessory-',locus)
  m$records <- c(m$records,list(r))
  t <- Filter(function(t) t$subject_id == 'patient-1' && t$locus == 'DRB1', input$a$typing)[[1L]]
  t$locus <- locus; t$alleles <- list(paste0(locus,'*01:01'))
  input$a$typing <- c(input$a$typing,list(t))
  indirect_with_reference(m,input)
}
# Contract-validation mock only: no historical binary/weights or real execution.
indirect_precomputed_contract <- function() {
  m <- indirect_fixture('binding.json')
  m$provider$evidence_kind <- 'precomputed_predictor_output'
  m$provider$build <- 'validation-only mock executable'
  m$provider$model <- 'validation-only mock model'
  artifact <- function(identity) list(identity=identity, reference=paste0('mock:',identity),
    sha256=digest::digest(identity,algo='sha256',serialize=FALSE))
  m$provenance <- list(kind='precomputed_predictor_output',
    query_artifact=artifact('validation-only query'), result_artifact=artifact('validation-only result'),
    executable_artifact=artifact(m$provider$build), model_artifact=artifact(m$provider$model),
    license=list(status='authorized_external_precomputed',terms='mock license status for validation only'),
    execution=list(origin='imported_external_execution',source='contract-validation mock; no real execution',
      provider_name='NetMHCIIpan',provider_version='3.2',command='mock BA invocation; never executed'),
    metadata=list(purpose='schema/runtime contract validation; not historical predictor evidence'))
  m
}
indirect_schema_valid <- function(x, name='indirect-binding-import', definition=NULL) {
  schema <- system.file('schema',paste0(name,'-1.0.0.json'),package='Neoantimon2')
  if (!is.null(definition)) {
    s <- jsonlite::read_json(schema); s$`$ref` <- paste0('#/$defs/',definition)
    schema <- jsonlite::toJSON(s,auto_unbox=TRUE,null='null')
  }
  validate <- jsonvalidate::json_validator(schema,engine='ajv')
  as.logical(validate(jsonlite::toJSON(x,auto_unbox=TRUE,null='null',digits=NA)))
}
