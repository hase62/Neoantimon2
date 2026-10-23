setup_shebang <- function(path) {
  con <- file(path,'rb');on.exit(close(con))
  raw <- readBin(con,'raw',512L)
  if(length(raw)<2L || !identical(raw[1:2],charToRaw('#!'))) return(NULL)
  line <- strsplit(rawToChar(raw),'\n',fixed=TRUE)[[1L]][[1L]]
  tokens <- strsplit(trimws(substring(line,3L)),'[[:space:]]+')[[1L]]
  if(basename(tokens[[1L]])=='env') {
    if(length(tokens)<2L || startsWith(tokens[[2L]],'-')) return(NULL)
    setup_executable(tokens[[2L]])
  } else setup_executable(tokens[[1L]])
}
setup_python_probe <- function(python,timeout) {
  path <- setup_executable(python)
  if(is.null(path)) return(list(success=FALSE,error='Python interpreter missing; supply python=',status='PYTHON_MISSING'))
  script <- system.file('setup','check_mhcflurry.py',package='Neoantimon2')
  if(!nzchar(script)) stop('Installed Python probe missing; reinstall this development build')
  probe <- setup_command(path,script,timeout)
  lines <- grep('^NEO_SETUP_JSON:',probe$output,value=TRUE)
  if(!length(lines)) return(list(success=FALSE,status='PYTHON_PROBE_FAILED',
    error=paste(probe$output,collapse='\n'),python=path))
  result <- jsonlite::fromJSON(sub('^NEO_SETUP_JSON:','',utils::tail(lines,1L)))
  result$status <- if(isTRUE(result$success) && probe$status==0L) 'PYTHON_CHECKED' else 'PYTHON_DEPENDENCY_MISSING'
  result$success <- isTRUE(result$success) && probe$status==0L
  result
}
setup_predictor_capability <- function(provider,executable,model_path,timeout) {
  spec <- predictor_registry()[[provider]]
  path <- setup_executable(executable,spec$executable)
  if(is.null(path)) return(list(available=FALSE,availability_status='predictor_unavailable',executable_path=NULL))
  probe <- setup_command(path,if(provider=='mhcflurry') '--version' else '-h',timeout)
  version <- predictor_version_from_text(probe$output,provider)
  identity <- predictor_model_identity(model_path)
  model_ok <- !spec$requires_models || (!is.null(identity) && length(identity$files)>0L)
  # Some vendor help commands return nonzero after printing their version.
  # Retain that diagnostic without treating version discovery as a smoke pass.
  reason <- if(probe$status!=0L && is.null(version)) 'execution_probe_failed' else if(is.null(version)) 'version_unresolved' else
    if(!predictor_version_supported(provider,version)) 'unsupported_version' else if(!model_ok) 'model_unavailable' else 'available'
  c(spec,list(available=reason=='available',availability_status=reason,executable_path=path,
    executable_checksum=setup_sha(path),predictor_version=version,
    version_status=if(is.null(version)) 'unknown' else 'detected',model_identity=identity,
    model_version=NULL,model_version_status='unknown',adapter_version=predictor_adapter_version(),
    probe_exit_status=probe$status,probe_output=probe$output))
}
setup_predictor_smoke <- function(provider,capability,request,root,timeout) {
  ii <- provider=='netmhciipan'
  candidate <- list(candidate_id='setup-candidate',sample_peptide_id='setup-peptide',
    origin_ids=list('setup-origin'),ligand_sequence=if(ii) 'ACDEFGHIKLMNPQR' else 'ACDEFGHIK',
    mhc_key=if(ii) 'HLA-DRA*01:01~HLA-DRB1*01:01' else 'HLA-A*02:01',
    hla_class=if(ii) 'II' else 'I',hla_allele=if(ii) NULL else 'HLA-A*02:01',
    hla_alpha=if(ii) 'HLA-DRA*01:01' else NULL,hla_beta=if(ii) 'HLA-DRB1*01:01' else NULL,
    species='NCBITaxon:9606',analysis_context='somatic:setup',sample_id='setup-sample')
  prepared <- predictor_prepare(list(candidate),provider,request$task,capability)
  options <- request$options;options$timeout_seconds <- timeout
  run <- tempfile('neo-predictor-check-',tmpdir=setup_dest(root,TRUE))
  result <- predictor_run(prepared,run,options)
  good <- length(result$evidence)>0L && !length(result$statuses) &&
    any(vapply(result$evidence,function(e) identical(e$task,request$task),TRUE))
  list(success=good,status=if(good) 'PASSED' else 'FAILED',artifact_root=run,
    statuses=result$statuses,note='Synthetic technical query, not scientific accuracy validation')
}

neo_check_predictors <- function(config=NULL,providers=NULL,python=NULL,smoke_test=FALSE,
                                 run_root=tempdir(),timeout=30) {
  setup_flag(smoke_test,'smoke_test')
  if(!is.numeric(timeout) || length(timeout)!=1L || !is.finite(timeout) || timeout<1) stop('Invalid timeout')
  if(!is.null(config) && !inherits(config,'neo_config')) stop('config must be a neo_config object')
  supported <- setdiff(names(predictor_registry()),'fixture')
  requests <- if(is.null(config)) list() else config$predictors
  if(is.null(providers)) providers <- if(length(requests)) unique(vapply(requests,`[[`,'','provider')) else supported
  if(!is.character(providers) || !length(providers) || anyNA(providers) || any(!providers %in% supported))
    stop('Choose supported native predictors: ',paste(supported,collapse=', '))
  selected <- list()
  for(provider in unique(providers)) {
    matches <- Filter(function(r) identical(r$provider,provider),requests)
    if(!length(matches)) matches <- list(list(provider=provider,
      task=predictor_registry()[[provider]]$tasks[[1L]],options=list()))
    selected <- c(selected,matches)
  }
  details <- lapply(selected,function(request) {
    provider <- request$provider
    cap <- tryCatch(setup_predictor_capability(provider,request$executable,request$model_path,timeout),
      error=function(e)list(available=FALSE,availability_status='probe_failed',message=conditionMessage(e)))
    runtime <- NULL
    if(provider=='mhcflurry' && !is.null(cap$executable_path)) {
      interpreter <- if(!is.null(python)) setup_executable(python) else setup_shebang(cap$executable_path)
      if(!is.null(interpreter) && !grepl('^python[0-9.]*([.]exe)?$',basename(interpreter))) interpreter <- NULL
      runtime <- setup_python_probe(interpreter,timeout)
      if(!isTRUE(runtime$success) && cap$available) {
        cap$available <- FALSE;cap$availability_status <- runtime$status
      }
      if(isTRUE(runtime$success) && !is.null(cap$predictor_version) &&
          !identical(runtime$version,cap$predictor_version)) {
        cap$available <- FALSE;cap$availability_status <- 'python_cli_version_mismatch'
      }
    }
    smoke <- list(status='NOT_RUN',success=NA)
    if(smoke_test && isTRUE(cap$available)) smoke <- tryCatch(
      setup_predictor_smoke(provider,cap,request,run_root,timeout),
      error=function(e)list(status='FAILED',success=FALSE,message=conditionMessage(e)))
    list(provider=provider,task=request$task,capability=cap,python=runtime,smoke_test=smoke,
      status=if(!isTRUE(cap$available)) cap$availability_status else if(smoke$status=='FAILED') 'SMOKE_FAILED' else
        if(smoke$status=='PASSED') 'SMOKE_PASSED' else 'CONFIGURED_UNTESTED')
  })
  summary <- do.call(rbind,lapply(details,function(d) data.frame(provider=d$provider,task=d$task,
    status=d$status,version=if(is.null(d$capability$predictor_version)) NA_character_ else d$capability$predictor_version,
    executable=if(is.null(d$capability$executable_path)) NA_character_ else d$capability$executable_path,
    stringsAsFactors=FALSE)))
  list(predictors=summary,details=details,note='Version/configuration checks do not establish prediction or scientific compatibility; smoke checks are explicit')
}

neo_setup_predictor <- function(provider,executable=NULL,model_path=NULL,task=NULL,options=list(),
                                python=NULL,environment=NULL,install=FALSE,version=NULL,
                                dependencies=character(),smoke_test=FALSE) {
  scalar_text(provider,'provider');setup_flag(install,'install')
  setup_flag(smoke_test,'smoke_test')
  if(!provider %in% setdiff(names(predictor_registry()),'fixture')) stop('Unsupported native predictor')
  if(is.null(task)) task <- predictor_registry()[[provider]]$tasks[[1L]]
  if(!is.character(task) || length(task)!=1L || !task %in% predictor_registry()[[provider]]$tasks)
    stop('Task is outside the existing provider contract')
  if(!is.null(executable)) scalar_text(executable,'executable')
  if(!is.null(model_path)) scalar_text(model_path,'model_path')
  predictor_options_validate(provider,task,options)
  installation <- list(status='NOT_REQUESTED')
  if(install) {
    if(provider!='mhcflurry') stop('Obtain licensed predictors from their vendor; configure an existing executable')
    scalar_text(environment,'new environment');scalar_text(version,'pinned MHCflurry version')
    if(!grepl('^2[.][0-9]+([.][0-9]+)*[A-Za-z0-9]*$',version)) stop('Choose an explicit supported 2.x release')
    if(!is.character(dependencies) || anyNA(dependencies) ||
      any(!grepl('^[A-Za-z0-9][A-Za-z0-9_.-]*==[0-9][A-Za-z0-9._+-]*$',dependencies)))
      stop('Dependencies must be explicit package==version pins')
    if(!is.null(executable)) stop('Choose a new environment or an existing executable, not both')
    env <- resolve_path(environment,getwd())
    if(file.exists(env)) stop('Environment already exists and is retained: ',env)
    base <- setup_python(python)
    dir.create(dirname(env),recursive=TRUE,showWarnings=FALSE)
    if(!dir.create(env,showWarnings=FALSE)) stop('Cannot create new environment')
    marker <- file.path(env,'.neo-setup-incomplete')
    writeLines('Installation not complete; do not use as a ready predictor',marker)
    result <- setup_command(base,c('-m','venv',env))
    if(result$status!=0L) stop('Environment creation failed; retained incomplete at ',env,': ',paste(result$output,collapse='\n'))
    bin <- file.path(env,if(.Platform$OS.type=='windows') 'Scripts' else 'bin')
    python <- file.path(bin,if(.Platform$OS.type=='windows') 'python.exe' else 'python')
    result <- setup_command(python,c('-m','pip','install',paste0('mhcflurry==',version),dependencies))
    writeLines(result$output,file.path(env,'neo-setup-install.log'))
    if(result$status!=0L) stop('Package installation failed; retained incomplete at ',env)
    runtime <- setup_python_probe(python,60)
    if(!isTRUE(runtime$success) || !identical(runtime$version,version))
      stop('Python/runtime validation failed; retained incomplete at ',env,': ',runtime$error)
    executable <- file.path(bin,if(.Platform$OS.type=='windows') 'mhcflurry-predict.exe' else 'mhcflurry-predict')
    if(is.null(setup_executable(executable))) stop('Predictor entry point missing; retained incomplete at ',env)
    cli <- setup_predictor_capability(provider,executable,model_path,60)
    if(!cli$availability_status %in% c('available','model_unavailable'))
      stop('Predictor execution/version check failed; retained incomplete at ',env)
    frozen <- setup_command(python,c('-m','pip','freeze'))
    if(frozen$status!=0L) stop('Cannot record installed environment; retained incomplete at ',env)
    writeLines(frozen$output,file.path(env,'neo-setup-requirements.txt'))
    unlink(marker)
    installation <- list(status='ENVIRONMENT_PREPARED',path=env,version=version,
      note='Models are separate; this is not a verified prediction installation')
  } else if(!is.null(environment) || !is.null(version) || length(dependencies))
    stop('Environment/package changes require explicit install=TRUE')
  request <- list(provider=provider,task=task,executable=executable,model_path=model_path,options=options)
  cfg <- neo_config(env=character(),predictors=list(request))
  list(predictors=cfg$predictors,checks=neo_check_predictors(cfg,python=python,smoke_test=smoke_test),
    installation=installation,python=python)
}
