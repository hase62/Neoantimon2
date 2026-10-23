test_that('missing native tools return diagnostics without changing configuration', {
  cfg <- neo_config(env=character(),predictors=list(list(provider='netmhcpan',task='mhc1_binding',
    executable='/absent/netMHCpan',options=list())))
  before <- cfg;root <- tempfile()
  r <- neo_check_predictors(cfg,run_root=root)
  expect_identical(r$predictors$status,'predictor_unavailable')
  expect_identical(cfg,before);expect_false(dir.exists(root))
  expect_error(neo_check_predictors(providers='unknown'),'supported')
  expect_error(neo_check_predictors(providers='fixture'),'supported')
  expect_error(neo_check_predictors(timeout=0),'timeout')
})

test_that('version probes preserve the existing native version contract', {
  cases <- list(netmhcpan=c('4.1','4.2e','4.2.1','3.0'),netmhciipan=c('4.3','4.3k','3.2'),mhcflurry=c('2.0.6','2.3.16','3.0'))
  for(provider in names(cases)) for(version in cases[[provider]]) {
    prefix <- Neoantimon2:::predictor_registry()[[provider]]$name
    path <- setup_test_executable(paste(prefix,version))
    models <- tempfile();dir.create(models);writeLines('synthetic model identity',file.path(models,'weights'))
    before <- Neoantimon2:::predictor_capability(provider,path,models)
    after <- Neoantimon2:::setup_predictor_capability(provider,path,models,10)
    expect_identical(after$predictor_version,before$predictor_version)
    expect_identical(after$available,before$available)
    expect_identical(after$available,!version %in% c('3.0','3.2'))
    expect_identical(after$availability_status,before$availability_status)
    expect_identical(after$model_identity,before$model_identity)
  }
})

test_that('file existence alone does not establish executable readiness', {
  path <- setup_test_executable('not a predictor version')
  cfg <- neo_config(env=character(),predictors=list(list(provider='netmhcpan',task='mhc1_binding',executable=path,options=list())))
  expect_identical(neo_check_predictors(cfg)$predictors$status,'version_unresolved')
  Sys.chmod(path,'0644')
  expect_identical(neo_check_predictors(cfg)$predictors$status,'predictor_unavailable')
  Sys.chmod(path,'0755')
  path <- setup_test_executable(body="printf '%s\\n' 'NetMHCpan 4.2'; exit 1")
  cfg$predictors[[1L]]$executable <- path
  checked <- neo_check_predictors(cfg)
  expect_identical(checked$predictors$status,'CONFIGURED_UNTESTED')
  expect_identical(checked$details[[1L]]$capability$probe_exit_status,1L)
})

test_that('configuration of an existing licensed installation does not install or mutate it', {
  path <- setup_test_executable('NetMHCpan 4.2e');sha <- Neoantimon2:::setup_sha(path)
  result <- neo_setup_predictor('netmhcpan',executable=path)
  expect_identical(result$installation$status,'NOT_REQUESTED')
  expect_identical(result$checks$predictors$status,'CONFIGURED_UNTESTED')
  expect_identical(result$checks$predictors$version,'4.2e')
  expect_identical(Neoantimon2:::setup_sha(path),sha)
  expect_s3_class(neo_config(env=character(),predictors=result$predictors),'neo_config')
  expect_error(neo_setup_predictor('netmhcpan',install=TRUE),'licensed')
  expect_error(neo_setup_predictor('netmhcpan',executable=path,task='mhc2_binding'),'Task')
})

test_that('a synthetic smoke run exercises the existing adapter and preserves a valid zero', {
  fixture <- readLines(test_path('..','fixtures','predictor','netmhcpan-4.1.txt'))
  raw <- fixture[c(1L,2L,which(grepl('ACDEFGHIK synthetic',fixture,fixed=TRUE)))]
  raw[3L] <- sub('85 NA$','0 NA',raw[3L])
  path <- setup_test_executable(body=c('if [ "$1" = "-h" ]; then',
    "printf '%s\\n' 'NetMHCpan 4.1'; exit 0; fi", "cat <<'NEO_OUTPUT'",raw,'NEO_OUTPUT'))
  cfg <- neo_config(env=character(),predictors=list(list(provider='netmhcpan',task='mhc1_binding',executable=path,options=list())))
  r <- neo_check_predictors(cfg,smoke_test=TRUE,run_root=tempdir())
  expect_identical(r$predictors$status,'SMOKE_PASSED')
  expect_true(file.exists(file.path(r$details[[1L]]$smoke_test$artifact_root,
    list.dirs(r$details[[1L]]$smoke_test$artifact_root,recursive=FALSE,full.names=FALSE)[[1L]],'result.rds')))
  run <- readRDS(list.files(r$details[[1L]]$smoke_test$artifact_root,'result.rds',recursive=TRUE,full.names=TRUE)[[1L]])
  expect_true(any(vapply(run$evidence,function(e)identical(e$raw_score,0),TRUE)))
  bad <- setup_test_executable(body=c('if [ "$1" = "-h" ]; then',
    "printf '%s\\n' 'NetMHCpan 4.1'; exit 0; fi", "printf '%s\\n' 'invalid output'"))
  cfg$predictors[[1L]]$executable <- bad
  expect_identical(neo_check_predictors(cfg,smoke_test=TRUE)$predictors$status,'SMOKE_FAILED')
})

test_that('MHCflurry needs models, a runnable Python and matching dependency/version checks', {
  path <- setup_test_executable('MHCflurry 2.0.6')
  cfg <- neo_config(env=character(),predictors=list(list(provider='mhcflurry',task='mhc1_binding',executable=path,options=list())))
  expect_identical(neo_check_predictors(cfg,python='/absent/python')$predictors$status,'model_unavailable')
  models <- tempfile();dir.create(models);writeLines('synthetic',file.path(models,'model'))
  cfg$predictors[[1L]]$model_path <- models
  expect_identical(neo_check_predictors(cfg,python='/absent/python')$predictors$status,'PYTHON_MISSING')
  testthat::local_mocked_bindings(setup_python_probe=function(...)list(success=FALSE,status='PYTHON_DEPENDENCY_MISSING',error='missing tensorflow'),.package='Neoantimon2')
  expect_identical(neo_check_predictors(cfg)$predictors$status,'PYTHON_DEPENDENCY_MISSING')
})

test_that('MHCflurry runtime version mismatches are explicit', {
  path <- setup_test_executable('MHCflurry 2.0.6');models <- tempfile();dir.create(models);writeLines('synthetic',file.path(models,'model'))
  cfg <- neo_config(env=character(),predictors=list(list(provider='mhcflurry',task='mhc1_binding',executable=path,model_path=models,options=list())))
  testthat::local_mocked_bindings(setup_python_probe=function(...)list(success=TRUE,version='2.3.16',status='PYTHON_CHECKED'),.package='Neoantimon2')
  expect_identical(neo_check_predictors(cfg)$predictors$status,'python_cli_version_mismatch')
})

test_that('environment installation requires explicit consent, pins and a new destination', {
  env <- tempfile();dir.create(env);writeLines('keep',file.path(env,'user-file'))
  expect_error(neo_setup_predictor('mhcflurry',environment=env),'install=TRUE')
  expect_error(neo_setup_predictor('mhcflurry',install=TRUE,environment=env,version='2.0.6'),'exists')
  expect_identical(readLines(file.path(env,'user-file')),'keep')
  expect_error(neo_setup_predictor('mhcflurry',install=TRUE,environment=tempfile(),version='latest'),'explicit')
  expect_error(neo_setup_predictor('mhcflurry',install=TRUE,environment=tempfile(),version='2.0.6',dependencies='--user'),'pins')
})

test_that('failed isolated installations retain an incomplete marker without touching the active environment', {
  env <- tempfile();before <- Sys.getenv(c('PATH','VIRTUAL_ENV','PYTHONPATH'))
  testthat::local_mocked_bindings(setup_python=function(...)'test-python',
    setup_command=function(executable,args=character(),timeout=0) {
      if('venv' %in% args) return(list(status=0L,output=character()))
      list(status=1L,output='simulated pip failure')
    },.package='Neoantimon2')
  expect_error(neo_setup_predictor('mhcflurry',install=TRUE,environment=env,version='2.0.6'),'incomplete')
  expect_true(file.exists(file.path(env,'.neo-setup-incomplete')))
  expect_identical(Sys.getenv(c('PATH','VIRTUAL_ENV','PYTHONPATH')),before)
  expect_true(file.exists(file.path(env,'neo-setup-install.log')))
})

test_that('Python venv interpreter paths preserve the environment leaf symlink', {
  skip_on_os('windows')
  root <- tempfile();dir.create(root);bin <- file.path(root,'bin');dir.create(bin)
  base <- Sys.which('Rscript');path <- file.path(bin,'python')
  expect_true(file.symlink(base,path))
  expect_identical(Neoantimon2:::setup_executable(path),file.path(normalizePath(bin),'python'))
  expect_false(identical(Neoantimon2:::setup_executable(path),normalizePath(path)))
})

test_that('successful isolated environment setup records exact pins but keeps missing models explicit', {
  skip_on_os('windows')
  env <- tempfile();calls <- list();before <- Sys.getenv(c('PATH','VIRTUAL_ENV','PYTHONPATH'))
  testthat::local_mocked_bindings(setup_python=function(...)'test-python',
    setup_python_probe=function(...)list(success=TRUE,version='2.0.6',status='PYTHON_CHECKED'),
    setup_command=function(executable,args=character(),timeout=0) {
      calls[[length(calls)+1L]] <<- args
      if('venv' %in% args) {
        bin <- file.path(env,'bin');dir.create(bin)
        for(name in c('python','mhcflurry-predict')) {
          path <- file.path(bin,name);writeLines('#!/bin/sh',path);Sys.chmod(path,'0755')
        }
      }
      list(status=0L,output=if(identical(args,'--version')) 'MHCflurry 2.0.6' else 'mhcflurry==2.0.6')
    },.package='Neoantimon2')
  r <- neo_setup_predictor('mhcflurry',install=TRUE,environment=env,version='2.0.6',
    dependencies='tensorflow==2.15.0')
  expect_identical(r$installation$status,'ENVIRONMENT_PREPARED')
  expect_false(file.exists(file.path(env,'.neo-setup-incomplete')))
  expect_true(file.exists(file.path(env,'neo-setup-requirements.txt')))
  expect_true(any(vapply(calls,function(a) all(c('mhcflurry==2.0.6','tensorflow==2.15.0') %in% a),TRUE)))
  expect_identical(r$checks$predictors$status,'model_unavailable')
  expect_identical(Sys.getenv(c('PATH','VIRTUAL_ENV','PYTHONPATH')),before)
})
