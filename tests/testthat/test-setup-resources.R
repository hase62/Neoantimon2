test_that('the installed catalog preserves profiles, identities and runtime assets', {
  c <- Neoantimon2:::setup_catalog()
  expect_identical(unlist(c$profiles$minimum,use.names=FALSE),c('genome','gencode','mane'))
  expect_setequal(unlist(c$profiles$full,use.names=FALSE),names(c$resources))
  expect_length(names(c$resources),10L)
  roles <- unlist(lapply(c$resources,function(s)lapply(s$files,`[[`,'role')),use.names=FALSE)
  expect_setequal(roles,Neoantimon2:::reference_roles('gencode'))
  expect_identical(Neoantimon2:::setup_selection(c,'minimum',c('mane','imgthla')),c('genome','gencode','mane','imgthla'))
  expect_error(Neoantimon2:::setup_selection(c,NULL,NULL),'Choose')
  expect_error(Neoantimon2:::setup_selection(c,'unknown',NULL),'Unknown')
  expect_error(Neoantimon2:::setup_selection(c,NULL,'unknown'),'Unknown')
  expect_true(file.exists(system.file('setup','download_public_data.py',package='Neoantimon2')))
  text <- paste(readLines(system.file('setup','download_public_data.py',package='Neoantimon2')),collapse='\n')
  expect_match(text,attr(c,'sha256'),fixed=TRUE)
})

test_that('the twelve released API signatures remain unchanged alongside four explicit setup additions', {
  baseline <- jsonlite::read_json(test_path('..','fixtures','setup','released-api-signatures.json'))$functions
  expect_length(baseline,12L)
  expect_setequal(setdiff(getNamespaceExports('Neoantimon2'),names(baseline)),
    c('neo_download_resources','neo_check_resources','neo_check_predictors','neo_setup_predictor'))
  for(name in names(baseline)) {
    f <- getExportedValue('Neoantimon2',name)
    expect_identical(as.list(names(formals(f))),baseline[[name]]$arguments)
    defaults <- lapply(formals(f),function(x)if(missing(x)) '<missing>' else
      paste(deparse(x,width.cutoff=500L),collapse='\n'))
    expect_identical(defaults,baseline[[name]]$defaults)
  }
})

test_that('native R downloads prepare the minimum roles and bind existing configuration', {
  f <- setup_test_catalog();setup_test_bind(f)
  dest <- tempfile('setup-destination-');setup_test_seed_fai(f,dest)
  report <- neo_download_resources(profile='minimum',dest=dest)
  expect_true(report$ok);expect_length(report$resources,9L)
  expect_true(all(report$files$status=='OK'))
  expect_true(all(report$files$integrity=='VERIFIED_RECEIPT'))
  expect_identical(report$reference_validation$status,'NOT_RUN')
  cfg <- neo_config(env=character(),resources=report$resources,scientific=report$scientific,
    paths=list(data_root=dest,run_root=file.path(dest,'runs')))
  checked <- neo_check_resources(config=cfg,validate_reference=TRUE)
  expect_true(checked$ok);expect_identical(checked$reference_validation$status,'PASSED')
  expect_false(dir.exists(file.path(dest,'.neoantimon2-setup-lock')))
  expect_length(list.files(dest,pattern='^\\.neo-setup-',all.files=TRUE),0L)
  before <- vapply(report$files$path,Neoantimon2:::setup_sha,'')
  again <- neo_download_resources(profile='minimum',dest=dest,python='/missing/python')
  expect_identical(vapply(again$files$path,Neoantimon2:::setup_sha,''),before)
  expect_true(again$ok)
})

test_that('checks distinguish missing, corrupt and changed recorded content', {
  f <- setup_test_catalog();setup_test_bind(f);dest <- tempfile()
  missing <- neo_check_resources(resources='mane',dest=dest)
  expect_false(missing$ok);expect_true(all(missing$files$status=='MISSING'));expect_false(dir.exists(dest))
  r <- neo_download_resources(resources='mane',dest=dest)
  path <- r$files$path[[1L]]
  bytes <- readBin(path,'raw',n=file.info(path)$size)
  writeBin(bytes[-seq.int(length(bytes)-4L,length(bytes))],path)
  bad <- neo_check_resources(resources='mane',dest=dest)
  expect_false(bad$ok);expect_identical(bad$files$status[[1L]],'INVALID')
  expect_error(neo_download_resources(resources='mane',dest=dest),'retained')
  fixed <- neo_download_resources(resources='mane',dest=dest,overwrite=TRUE)
  expect_true(fixed$ok)
  receipt <- file.path(dest,'acquisition/r-setup-mane.json')
  x <- jsonlite::read_json(receipt);x$artifacts[[1L]]$sha256 <- paste(rep('0',64),collapse='')
  jsonlite::write_json(x,receipt,auto_unbox=TRUE)
  mismatch <- neo_check_resources(resources='mane',dest=dest,checksum=FALSE)
  expect_false(mismatch$ok);expect_match(mismatch$files$message[[1L]],'SHA256 mismatch')
})

test_that('failed or interrupted acquisition preserves existing data and cleans only owned staging', {
  f <- setup_test_catalog();setup_test_bind(f);dest <- tempfile();dir.create(dest)
  old <- file.path(dest,f$catalog$resources$mane$files[[1L]]$path)
  dir.create(dirname(old),recursive=TRUE);writeLines('previous content',old)
  unrelated <- file.path(dest,'.neo-setup-unrelated');dir.create(unrelated);writeLines('keep',file.path(unrelated,'keep'))
  testthat::local_mocked_bindings(setup_download_file=function(url,path) {
    writeLines('interrupted partial bytes',path);stop('simulated interrupted transfer')
  },.package='Neoantimon2')
  expect_error(neo_download_resources(resources='mane',dest=dest,overwrite=TRUE),'interrupted')
  expect_identical(readLines(old),'previous content')
  expect_true(file.exists(file.path(unrelated,'keep')))
  expect_false(file.exists(file.path(dest,'acquisition/r-setup-mane.json')))
  expect_false(dir.exists(file.path(dest,'.neoantimon2-setup-lock')))
})

test_that('destination, selection and concurrent-operation errors are actionable', {
  f <- setup_test_catalog();setup_test_bind(f)
  dest <- tempfile();writeLines('user file',dest)
  expect_error(neo_download_resources(resources='mane',dest=dest),'not a directory')
  expect_identical(readLines(dest),'user file')
  expect_error(neo_download_resources(profile='bad',dest=tempfile()),'Unknown')
  expect_error(neo_download_resources(dest=tempfile()),'Choose')
  root <- tempfile();dir.create(root);dir.create(file.path(root,'.neoantimon2-setup-lock'))
  expect_error(neo_download_resources(resources='mane',dest=root),'locked')
  expect_true(dir.exists(file.path(root,'.neoantimon2-setup-lock')))
})

test_that('configured resources preserve the original synthetic reference contract and detect inconsistencies', {
  cfg <- setup_test_config();before <- cfg
  result <- neo_check_resources(config=cfg,validate_reference=TRUE)
  expect_true(result$ok);expect_identical(result$reference_validation$status,'PASSED')
  expect_identical(cfg,before)
  cfg$resources[['reference.genome_fasta']]$checksum$value <- paste(rep('0',64),collapse='')
  expect_false(neo_check_resources(config=cfg,checksum=FALSE)$ok)
  cfg <- before;cfg$resources[['reference.mane_gtf']] <- NULL
  expect_false(neo_check_resources(config=cfg)$ok)
  expect_identical(neo_check_resources(config=cfg)$files$status[9L],'MISSING')
  cfg <- before;cfg$resources[['reference.genome_fai']]$path <- tempfile()
  expect_false(neo_check_resources(config=cfg)$ok)
  cfg <- before;cfg$scientific$genome_assembly <- 'GRCh38.p13'
  expect_false(neo_check_resources(config=cfg)$ok)
  cfg <- before;cfg$resources[['reference.protein_fasta']]$version <- '49'
  expect_false(neo_check_resources(config=cfg)$ok)
  cfg <- before;cfg$resources[['reference.refFlat']] <- list(path=tempfile(),type='file')
  expect_false(neo_check_resources(config=cfg)$ok)
  expect_false(dir.exists(cfg$paths$run_root))
  expect_error(neo_check_resources(config=before,profile='minimum'),'Choose config')
})

test_that('FAI preparation delegates only the index operation and verifies its geometry', {
  f <- setup_test_catalog();setup_test_bind(f);dest <- tempfile()
  calls <- list()
  testthat::local_mocked_bindings(setup_executable=function(path,default=NULL)'synthetic-samtools',
    setup_command=function(executable,args=character(),timeout=0) {
      calls[[length(calls)+1L]] <<- args
      file.copy(file.path(f$source,'genome.fa.fai'),args[[3L]])
      list(status=0L,output=character())
    },.package='Neoantimon2')
  result <- neo_download_resources(resources='genome',dest=dest)
  expect_true(result$ok);expect_length(calls,1L)
  expect_identical(calls[[1L]][1:2],c('faidx','-o'))
  expect_identical(Neoantimon2:::setup_sha(result$resources[['reference.genome_fai']]$path),
    Neoantimon2:::setup_sha(file.path(f$source,'genome.fa.fai')))
})

test_that('advanced resources use installed runtime assets and read-only checking', {
  expect_error(Neoantimon2:::setup_python('/missing/python'),'Python 3')
  captured <- list()
  catalog <- Neoantimon2:::setup_catalog();catalog$resources$imgthla$tools$python <- list()
  testthat::local_mocked_bindings(setup_catalog=function()catalog,
    setup_python=function(p)'fake-python',setup_command=function(executable,args,timeout=0) {
    captured[[length(captured)+1L]] <<- args;list(status=0L,output=character())
  },.package='Neoantimon2')
  Neoantimon2:::setup_python_resource('imgthla',tempdir(),NULL,check=TRUE)
  expect_true(file.exists(captured[[1L]][[1L]]))
  expect_true('--check-only' %in% captured[[1L]])
  expect_true('--no-overwrite' %in% captured[[1L]])
  expect_true('--catalog' %in% captured[[1L]])
  expect_false(grepl('/scripts/',captured[[1L]][[1L]],fixed=TRUE))
})

test_that('catalog interpretation rejects malformed paths, checksums and profiles', {
  c <- Neoantimon2:::setup_catalog()
  wrong <- c;wrong$resources$mane$files[[1L]]$path <- '../escape'
  expect_error(Neoantimon2:::setup_validate_catalog(wrong),'catalog')
  wrong <- c;wrong$resources$mane$files[[1L]]$sha256 <- 'bad'
  expect_error(Neoantimon2:::setup_validate_catalog(wrong),'checksum')
  wrong <- c;wrong$profiles$full <- list('missing')
  expect_error(Neoantimon2:::setup_validate_catalog(wrong),'mapping')
  wrong <- c;wrong$resources$mane$files[[1L]]$url <- 'https://user@example.invalid/file'
  expect_error(Neoantimon2:::setup_validate_catalog(wrong),'credentials')
})

test_that('missing index tools fail before any download and do not create completed artifacts', {
  f <- setup_test_catalog();setup_test_bind(f);dest <- tempfile()
  testthat::local_mocked_bindings(setup_executable=function(...)NULL,
    setup_download_file=function(...)stop('download should not run'),.package='Neoantimon2')
  expect_error(neo_download_resources(resources='genome',dest=dest),'samtools')
  expect_false(dir.exists(file.path(dest,'reference')))
  expect_false(dir.exists(file.path(dest,'.neoantimon2-setup-lock')))
  report <- neo_check_resources(resources='imgthla',dest=dest,python='/missing/python')
  expect_identical(report$files$status,'DEPENDENCY_MISSING')
})

test_that('a missing checkout is classified separately from an invalid checkout', {
  testthat::local_mocked_bindings(setup_python_resource=function(...)stop('Missing checkout'),.package='Neoantimon2')
  r <- neo_check_resources(resources='imgthla',dest=tempfile())
  expect_identical(r$files$status,'MISSING')
})
