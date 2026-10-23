setup_test_catalog <- function() {
  catalog <- Neoantimon2:::setup_catalog()
  mirror <- tempfile('setup-mirror-');dir.create(mirror)
  files <- c(genome_fasta='genome.fa',genome_fai='genome.fa.fai',annotation_gtf='annotation.gtf',
    transcript_fasta='transcripts.fa',protein_coding_transcript_fasta='coding.fa',protein_fasta='proteins.fa',
    selenocysteine_metadata='sec.tsv',mane_summary='mane.tsv',mane_gtf='mane.gtf')
  source <- normalizePath(test_path('..','fixtures','reference'))
  fai <- read.delim(file.path(source,'genome.fa.fai'),header=FALSE)
  catalog$reference$genome_contigs <- stats::setNames(as.list(fai[[2L]]),fai[[1L]])
  for(id in c('genome','gencode','mane')) for(i in seq_along(catalog$resources[[id]]$files)) {
    f <- catalog$resources[[id]]$files[[i]]
    if(f$action!='download') next
    role <- if(id=='genome') 'genome_fasta' else f$role
    path <- file.path(mirror,paste0(id,'-',i,'.gz'))
    con <- gzfile(path,'wb');input <- file.path(source,files[[role]])
    writeBin(readBin(input,'raw',n=file.info(input)$size),con);close(con)
    absolute <- normalizePath(path,winslash='/')
    f$url <- paste0(if(grepl('^[A-Za-z]:',absolute)) 'file:///' else 'file://',
      utils::URLencode(absolute,reserved=FALSE))
    catalog$resources[[id]]$files[[i]] <- f
  }
  jsonlite::write_json(catalog,file.path(mirror,'catalog.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
  attr(catalog,'sha256') <- Neoantimon2:::setup_sha(file.path(mirror,'catalog.json'))
  list(catalog=catalog,mirror=mirror,source=source,files=files)
}
setup_test_bind <- function(fixture,env=parent.frame()) {
  testthat::local_mocked_bindings(setup_catalog=function()fixture$catalog,.package='Neoantimon2',.env=env)
}
setup_test_seed_fai <- function(fixture,dest) {
  path <- file.path(dest,'reference/genome/GRCh38.primary_assembly.genome.fa.fai')
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  file.copy(file.path(fixture$source,'genome.fa.fai'),path)
}
setup_test_config <- function() {
  files <- c(genome_fasta='genome.fa',genome_fai='genome.fa.fai',annotation_gtf='annotation.gtf',
    transcript_fasta='transcripts.fa',protein_coding_transcript_fasta='coding.fa',protein_fasta='proteins.fa',
    selenocysteine_metadata='sec.tsv',mane_summary='mane.tsv',mane_gtf='mane.gtf')
  root <- system.file('extdata','quickstart',package='Neoantimon2')
  resources <- lapply(names(files),function(role)list(path=file.path(root,files[[role]]),required=TRUE,type='file',
    resource_id=if(role=='annotation_gtf') Neoantimon2:::reference_scope_policy() else NULL,
    version=if(startsWith(role,'genome')) 'GRCh38.p14' else if(startsWith(role,'mane')) '1.5' else '50',source='synthetic',
    checksum=list(algorithm='sha256',scope='file',value=Neoantimon2:::setup_sha(file.path(root,files[[role]])))))
  names(resources) <- paste0('reference.',names(files))
  neo_config(env=character(),paths=list(run_root=tempfile('setup-config-runs-')),
    scientific=list(genome_assembly='GRCh38.p14'),resources=resources)
}
setup_test_executable <- function(version='NetMHCpan version 4.2',body=NULL) {
  testthat::skip_on_os('windows')
  root <- tempfile('tool with space ');dir.create(root)
  path <- file.path(root,'predictor')
  if(is.null(body)) body <- paste0("printf '%s\\n' ",shQuote(version))
  writeLines(c('#!/bin/sh',body),path);Sys.chmod(path,'0755')
  path
}
