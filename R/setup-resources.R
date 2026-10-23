setup_validate_file <- function(path,item,assembly='GRCh38.p14') {
  if(!file.exists(path) || dir.exists(path) || file.access(path,4L)!=0L || file.info(path)$size<=0)
    stop('Missing, empty or unreadable file')
  con <- file(path,'rb');magic <- readBin(con,'raw',2L);close(con)
  compressed <- identical(magic,as.raw(c(31,139)))
  if(grepl('[.]gz$',path)!=compressed) stop('Compression/extension mismatch')
  if(compressed) {
    con <- gzfile(path,'rb')
    tryCatch(withCallingHandlers({
      repeat { bytes <- readBin(con,'raw',8L*1024L*1024L);if(!length(bytes)) break }
    },warning=function(w) stop('Invalid compressed resource: ',conditionMessage(w))),finally=close(con))
  }
  con <- reference_connection(path)
  lines <- tryCatch(readLines(con,n=200L,warn=FALSE),finally=close(con))
  body <- lines[nzchar(lines) & !startsWith(lines,'#')]
  kind <- item$kind
  if(!length(body)) stop('No readable records')
  if(kind %in% c('fasta','transcript','protein')) {
    if(!startsWith(body[[1L]],'>') || !any(!startsWith(body,'>'))) stop('Invalid FASTA')
    alphabet <- if(kind=='protein') '^[ABCDEFGHIKLMNPQRSTVWXYZUO*abcdefghiklmnpqrstvwxyzuo]+$' else
      '^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+$'
    if(any(!grepl(alphabet,body[!startsWith(body,'>')]))) stop('Invalid FASTA sequence')
    if(kind %in% c('transcript','protein')) {
      fields <- strsplit(substring(body[startsWith(body,'>')],2L),'|',fixed=TRUE)
      if(any(lengths(fields)<if(kind=='protein') 8L else 7L)) stop('Malformed GENCODE FASTA headers')
    }
  } else if(kind %in% c('gtf','gencode-gtf')) {
    if(any(lengths(strsplit(body,'\t',fixed=TRUE))!=9L)) stop('Malformed GTF')
    descriptions <- grep('^##description:',lines,value=TRUE)
    if(kind=='gencode-gtf' && !length(descriptions)) stop('GENCODE identity header missing')
    if(length(descriptions) && !is.null(item$version)) {
      version <- paste0('version ',item$version,'\\b')
      if(!any(grepl(version,descriptions) & grepl(sub('[.]p[0-9]+$','',assembly),descriptions,fixed=TRUE)))
        stop('GENCODE release/assembly header mismatch')
    }
  } else if(kind=='mane-summary') {
    if(!any(grepl('NCBI_GeneID',lines,fixed=TRUE) & grepl('MANE_status',lines,fixed=TRUE) &
      grepl('GRCh38_chr',lines,fixed=TRUE))) stop('MANE/GRCh38 summary header missing')
    if(any(lengths(strsplit(body,'\t',fixed=TRUE))!=14L)) stop('Malformed MANE summary')
  } else if(kind=='sec') {
    fields <- strsplit(body,'\t',fixed=TRUE)
    if(any(lengths(fields)!=2L) || any(!vapply(fields,function(f)
      (isTRUE(item$configured) || startsWith(f[[1L]],'ENST')) && nzchar(f[[1L]]) &&
        grepl('^[0-9]+$',f[[2L]]),TRUE))) stop('Malformed Sec metadata')
  } else if(kind=='fai') {
    fields <- strsplit(body,'\t',fixed=TRUE)
    if(any(lengths(fields)!=5L) || any(!vapply(fields,function(f)
      all(grepl('^[0-9]+$',f[-1L])),TRUE))) stop('Malformed FASTA index')
  }
  invisible(TRUE)
}
setup_receipt <- function(dest,id) {
  path <- file.path(dest,'acquisition',paste0('r-setup-',id,'.json'))
  if(!file.exists(path)) return(NULL)
  x <- jsonlite::read_json(path)
  if(!identical(x$schema_version,'neo-r-setup-1') || !identical(x$resource_id,id))
    stop('Unrecognized R acquisition receipt retained: ',path)
  x
}
setup_expected_sha <- function(item,receipt=NULL) {
  if(!is.null(item$sha256)) return(list(value=item$sha256,kind='VERIFIED_DECLARED'))
  if(!is.null(receipt)) for(f in receipt$artifacts)
    if(identical(f$path,item$path)) return(list(value=f$sha256,kind='VERIFIED_RECEIPT'))
  NULL
}
setup_native_check <- function(root,spec,checksum=TRUE,receipt=NULL,assembly='GRCh38.p14') {
  lapply(spec$files,function(item) {
    path <- file.path(root,item$path)
    result <- tryCatch({
      setup_validate_file(path,item,assembly)
      if(!is.null(receipt) && !identical(receipt$release,spec$release))
        stop('Acquisition receipt release differs from selected profile')
      expected <- setup_expected_sha(item,receipt)
      sha <- if(checksum || !is.null(expected)) setup_sha(path) else NA_character_
      if(!is.null(expected) && !identical(sha,expected$value)) stop('SHA256 mismatch')
      list(status='OK',sha256=sha,integrity=if(!is.null(expected)) expected$kind else
        if(checksum) 'OBSERVED_ONLY' else 'FORMAT_CHECKED',message='')
    },error=function(e) list(status=if(!file.exists(path)) 'MISSING' else 'INVALID',
      sha256=NA_character_,integrity='NOT_VERIFIED',message=conditionMessage(e)))
    c(list(path=path,role=if(is.null(item$role)) '' else item$role),result)
  })
}
setup_genome_check <- function(path,fai,expected=NULL) {
  index <- reference_fai(fai,path)
  if(!is.null(expected)) {
    lengths <- unlist(expected,use.names=TRUE)
    found <- index$length[match(names(lengths),index$contig)]
    if(anyNA(found) || any(found!=lengths)) stop('Genome assembly/contig lengths mismatch')
  }
  invisible(index)
}
setup_gzip_matches <- function(source,plain) {
  input <- gzfile(source,'rb');output <- file(plain,'rb')
  on.exit({close(input);close(output)},add=TRUE)
  repeat {
    a <- readBin(input,'raw',8L*1024L*1024L);b <- readBin(output,'raw',8L*1024L*1024L)
    if(!identical(a,b)) return(FALSE)
    if(!length(a)) return(TRUE)
  }
}
setup_download_file <- function(url,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  status <- utils::download.file(url,path,quiet=TRUE,mode='wb',
    method=if(capabilities('libcurl')) 'libcurl' else 'auto')
  if(status!=0L) stop('Download failed: ',url)
  invisible(path)
}
setup_native_download <- function(id,root,catalog,overwrite,timeout) {
  spec <- catalog$resources[[id]]
  receipt <- tryCatch(setup_receipt(root,id),error=function(e)if(overwrite) NULL else stop(e))
  checks <- setup_native_check(root,spec,receipt=receipt,assembly=catalog$reference$assembly)
  if(id=='genome') {
    fa <- which(vapply(spec$files,function(f) identical(f$role,'genome_fasta'),TRUE))
    fi <- which(vapply(spec$files,function(f) identical(f$role,'genome_fai'),TRUE))
    gz <- which(vapply(spec$files,function(f) identical(f$action,'download'),TRUE))
    problem <- NULL
    if(all(vapply(checks[c(gz,fa)],function(c)c$status=='OK',TRUE)) &&
        !setup_gzip_matches(checks[[gz]]$path,checks[[fa]]$path))
      problem <- 'Compressed/uncompressed genome mismatch'
    if(all(vapply(checks[c(fa,fi)],function(c)c$status=='OK',TRUE))) {
      geometry <- tryCatch({reference_fai(checks[[fi]]$path,checks[[fa]]$path);NULL},error=function(e)conditionMessage(e))
      if(!is.null(geometry)) {
        checks[[fi]]$status <- 'INVALID';checks[[fi]]$message <- geometry
      } else problem <- tryCatch({setup_genome_check(checks[[fa]]$path,checks[[fi]]$path,
        catalog$reference$genome_contigs);problem},error=function(e)conditionMessage(e))
    }
    if(!is.null(problem)) for(i in seq_along(checks)) {
      checks[[i]]$status <- 'INVALID';checks[[i]]$message <- problem
    }
  }
  invalid <- which(vapply(checks,function(x) x$status!='OK' && file.exists(x$path),TRUE))
  if(length(invalid) && !overwrite) stop('Existing resource retained: ',checks[[invalid[[1L]]]]$path,
    '; ',checks[[invalid[[1L]]]]$message,'; inspect it or choose overwrite=TRUE')
  if(id=='genome') {
    at <- which(vapply(spec$files,function(f) identical(f$role,'genome_fai'),TRUE))
    if(checks[[at]]$status!='OK' && is.null(setup_executable(NULL,'samtools')))
      setup_dependency('samtools is required to prepare the missing or invalid genome FAI')
  }
  stage <- tempfile('.neo-setup-',tmpdir=root);dir.create(stage)
  on.exit(unlink(stage,recursive=TRUE),add=TRUE)
  old <- options(timeout=timeout);on.exit(options(old),add=TRUE)
  selected <- list();obtained <- list()
  for(i in seq_along(spec$files)) {
    item <- spec$files[[i]];target <- file.path(root,item$path)
    if(checks[[i]]$status=='OK') {
      selected[[item$path]] <- target
      obtained[[i]] <- 'validated_existing'
      next
    }
    candidate <- file.path(stage,item$path)
    dir.create(dirname(candidate),recursive=TRUE,showWarnings=FALSE)
    if(item$action=='download') setup_download_file(item$url,candidate)
    else if(item$action=='gunzip') {
      input <- gzfile(selected[[item$parent]],'rb');output <- file(candidate,'wb')
      tryCatch(withCallingHandlers({repeat {
        bytes <- readBin(input,'raw',8L*1024L*1024L)
        if(!length(bytes)) break
        writeBin(bytes,output)
      }},warning=function(w) stop(conditionMessage(w))),finally={close(input);close(output)})
    } else if(item$action=='faidx') {
      tool <- setup_executable(NULL,'samtools')
      if(is.null(tool)) setup_dependency('samtools is required to prepare the missing genome FAI')
      result <- setup_command(tool,c('faidx','-o',candidate,selected[[item$parent]]))
      if(result$status!=0L) stop('samtools faidx failed: ',paste(result$output,collapse='\n'))
    }
    setup_validate_file(candidate,item,catalog$reference$assembly)
    if(!is.null(item$sha256) && !identical(setup_sha(candidate),item$sha256)) stop('Downloaded SHA256 mismatch')
    selected[[item$path]] <- candidate
    obtained[[i]] <- if(item$action=='download') 'downloaded' else 'prepared'
  }
  if(id=='genome') {
    genome <- spec$files[[which(vapply(spec$files,function(f) identical(f$role,'genome_fasta'),TRUE))]]$path
    fai <- spec$files[[which(vapply(spec$files,function(f) identical(f$role,'genome_fai'),TRUE))]]$path
    setup_genome_check(selected[[genome]],selected[[fai]],catalog$reference$genome_contigs)
  }
  artifacts <- lapply(seq_along(spec$files),function(i) {
    item <- spec$files[[i]];from <- selected[[item$path]]
    target <- file.path(root,item$path)
    sha <- setup_sha(from)
    if(!identical(from,target)) setup_promote(from,setup_write_target(root,item$path),overwrite)
    list(path=item$path,sha256=sha,bytes=as.numeric(file.info(target)$size),
      retrieval_mode=obtained[[i]],selected_source_url=item$url,
      original_download_date=if(obtained[[i]]=='downloaded') format(Sys.Date()) else NULL)
  })
  if(is.null(receipt) || any(vapply(obtained,function(x) x!='validated_existing',TRUE))) {
    record <- list(schema_version='neo-r-setup-1',resource_id=id,release=spec$release,
      catalog_sha256=attr(catalog,'sha256'),observed_at=format(Sys.time(),tz='UTC',usetz=TRUE),
      semantic_scope=spec$semantic_scope,license_url=spec$license_url,
      checksum_policy='Observed/receipt hashes are not publisher authenticity certificates',artifacts=artifacts)
    setup_json(record,setup_write_target(root,file.path('acquisition',paste0('r-setup-',id,'.json'))),overwrite)
  }
  invisible(artifacts)
}

neo_download_resources <- function(profile=NULL,resources=NULL,dest,overwrite=FALSE,python=NULL,timeout=600) {
  setup_flag(overwrite,'overwrite')
  if(!is.numeric(timeout) || length(timeout)!=1L || !is.finite(timeout) || timeout<1) stop('Invalid timeout')
  catalog <- setup_catalog();ids <- setup_selection(catalog,profile,resources)
  root <- setup_dest(dest,TRUE)
  lock <- file.path(root,'.neoantimon2-setup-lock')
  if(!dir.create(lock,showWarnings=FALSE)) stop('Setup destination is locked; check for an active/interrupted acquisition: ',lock)
  on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  for(id in ids) {
    if(id %in% c('genome','gencode','mane')) setup_native_download(id,root,catalog,overwrite,timeout)
    else setup_python_resource(id,root,python,overwrite=overwrite)
  }
  report <- neo_check_resources(resources=ids,dest=root,python=python)
  if(!report$ok) stop('Acquisition incomplete: ',paste(report$files$message[report$files$status!='OK'],collapse='; '))
  report
}

setup_resource_specs <- function(root,catalog,ids,checksum=TRUE,valid_paths=NULL) {
  out <- list()
  for(id in ids) for(f in catalog$resources[[id]]$files) if(!is.null(f$role)) {
    path <- file.path(root,f$path)
    if(!is.null(valid_paths) && !path %in% valid_paths) next
    if(!file.exists(path) || dir.exists(path)) next
    s <- list(path=path,required=TRUE,type='file',resource_id=f$resource_id,
      version=f$version,source=catalog$resources[[id]]$source_name,
      checksum=if(checksum) list(algorithm='sha256',scope='file',value=setup_sha(path)) else NULL)
    out[[paste0('reference.',f$role)]] <- s
  }
  out
}

neo_check_resources <- function(profile=NULL,resources=NULL,dest=NULL,config=NULL,
                                checksum=TRUE,validate_reference=FALSE,python=NULL) {
  setup_flag(checksum,'checksum');setup_flag(validate_reference,'validate_reference')
  if(!is.null(config)) {
    if(!is.null(profile) || !is.null(resources) || !is.null(dest)) stop('Choose config or a profile/resource destination, not both')
    return(setup_check_config(config,checksum,validate_reference))
  }
  catalog <- setup_catalog();ids <- setup_selection(catalog,profile,resources)
  root <- setup_dest(dest)
  rows <- list()
  for(id in ids) {
    spec <- catalog$resources[[id]]
    if(id %in% c('genome','gencode','mane')) {
      receipt <- tryCatch(setup_receipt(root,id),error=function(e)e)
      checks <- setup_native_check(root,spec,checksum,
        if(inherits(receipt,'error')) NULL else receipt,catalog$reference$assembly)
      extra <- if(inherits(receipt,'error')) conditionMessage(receipt) else NULL
      if(id=='genome' && all(vapply(checks,function(c)c$status=='OK',TRUE))) {
        f <- Filter(function(f) identical(f$role,'genome_fasta'),spec$files)[[1L]]
        i <- Filter(function(f) identical(f$role,'genome_fai'),spec$files)[[1L]]
        extra <- tryCatch({setup_genome_check(file.path(root,f$path),file.path(root,i$path),
          catalog$reference$genome_contigs)
          source <- Filter(function(f) identical(f$action,'download'),spec$files)[[1L]]
          if(!setup_gzip_matches(file.path(root,source$path),file.path(root,f$path)))
            stop('Compressed/uncompressed genome mismatch')
          extra},error=function(e)conditionMessage(e))
      }
      if(!is.null(extra)) for(i in seq_along(checks)) {
        checks[[i]]$status <- 'INVALID';checks[[i]]$message <- extra
      }
      for(c in checks) rows[[length(rows)+1L]] <- c(list(resource_id=id),c)
    } else {
      result <- tryCatch({setup_python_resource(id,root,python,check=TRUE,checksum=checksum);NULL},error=function(e)e)
      absent <- if(!is.null(spec$checkout)) !dir.exists(file.path(root,spec$checkout$path)) else
        any(!vapply(spec$files,function(f)file.exists(file.path(root,f$path)),TRUE))
      status <- if(is.null(result)) 'OK' else if(inherits(result,'neo_setup_dependency_error')) 'DEPENDENCY_MISSING' else
        if(absent) 'MISSING' else 'INVALID'
      rows[[length(rows)+1L]] <- list(resource_id=id,path=root,role='',status=status,
        sha256=NA_character_,integrity=if(is.null(result)) 'PROVIDER_VALIDATED' else 'NOT_VERIFIED',
        message=if(is.null(result)) '' else conditionMessage(result))
    }
  }
  files <- do.call(rbind,lapply(rows,function(r)as.data.frame(r,stringsAsFactors=FALSE)))
  specs <- setup_resource_specs(root,catalog,ids,checksum,files$path[files$status=='OK'])
  # Config uses genome_assembly, while the catalog describes assembly.
  cfg <- neo_config(env=character(),paths=list(data_root=root),
    scientific=list(species=catalog$reference$species,genome_assembly=catalog$reference$assembly),resources=specs)
  reference <- if(validate_reference) setup_reference_validation(cfg) else list(status='NOT_RUN',message='Full cross-resource reference validation not requested')
  list(ok=all(files$status=='OK') && reference$status!='FAILED',selected_resources=ids,destination=root,
    files=files,resources=specs,scientific=cfg$scientific,reference_validation=reference,
    catalog_sha256=attr(catalog,'sha256'))
}

setup_reference_validation <- function(config) {
  tryCatch({
    project <- structure(list(resources=resolve_resources(config),config=config,
      analysis_id='setup-reference-validation'),class=c('neo_project','list'))
    bundle <- reference_bundle(project)
    list(status='PASSED',message='Existing full reference validator passed',profile=bundle$profile)
  },error=function(e)list(status='FAILED',message=conditionMessage(e)))
}
setup_check_config <- function(config,checksum,validate_reference) {
  if(!inherits(config,'neo_config')) stop('config must be a neo_config object')
  names_all <- names(config$resources)
  legacy <- any(paste0('reference.',reference_roles('legacy')) %in% names_all)
  roles <- reference_roles(if(legacy) 'legacy' else 'gencode')
  kinds <- c(genome_fasta='fasta',genome_fai='fai',annotation_gtf='gtf',
    transcript_fasta='transcript',protein_coding_transcript_fasta='transcript',
    protein_fasta='protein',selenocysteine_metadata='sec',mane_summary='mane-summary',mane_gtf='gtf',
    refFlat='text',refMrna='fasta')
  rows <- lapply(roles,function(role) {
    logical <- paste0('reference.',role);spec <- config$resources[[logical]]
    path <- if(is.null(spec) || is.null(spec$path)) NA_character_ else
      tryCatch(resolve_path(spec$path,config$paths$data_root),error=function(e)NA_character_)
    result <- tryCatch({
      if(is.null(spec) || is.na(path)) stop('Required role/path missing: ',role)
      if(!identical(spec$type,'file')) stop('Reference role must be a file')
      if(!legacy && (is.null(spec$version) || is.null(config$scientific$genome_assembly)))
        stop('Declare reference version and genome assembly')
      if(role=='annotation_gtf' && !identical(spec$resource_id,reference_scope_policy()))
        stop('Declare the supported primary-reference scope policy')
      setup_validate_file(path,list(kind=kinds[[role]],version=spec$version,configured=TRUE),config$scientific$genome_assembly)
      sha <- if(checksum || !is.null(spec$checksum)) setup_sha(path) else NA_character_
      if(!is.null(spec$checksum) && !identical(spec$checksum$value,sha)) stop('Declared SHA256 mismatch')
      list(status='OK',integrity=if(!is.null(spec$checksum)) 'VERIFIED_DECLARED' else
        if(checksum) 'OBSERVED_ONLY' else 'FORMAT_CHECKED',sha256=sha,message='')
    },error=function(e)list(status=if(is.na(path) || !file.exists(path)) 'MISSING' else 'INVALID',
      integrity='NOT_VERIFIED',sha256=NA_character_,message=conditionMessage(e)))
    c(list(resource_id=logical,role=role,path=path),result)
  })
  files <- do.call(rbind,lapply(rows,function(r)as.data.frame(r,stringsAsFactors=FALSE)))
  conflict <- character()
  if(legacy && any(paste0('reference.',reference_roles('gencode')) %in% names_all))
    conflict <- roles
  if(!legacy) {
    version <- function(role) config$resources[[paste0('reference.',role)]]$version
    for(role in c('genome_fasta','genome_fai'))
      if(!identical(version(role),config$scientific$genome_assembly)) conflict <- c(conflict,role)
    for(role in c('transcript_fasta','protein_coding_transcript_fasta','protein_fasta','selenocysteine_metadata'))
      if(!identical(version(role),version('annotation_gtf'))) conflict <- c(conflict,role)
    if(!identical(version('mane_gtf'),version('mane_summary'))) conflict <- c(conflict,'mane_gtf')
    if(!identical(config$scientific$species,'NCBITaxon:9606')) conflict <- roles
  }
  at <- files$role %in% conflict & files$status!='MISSING'
  files$status[at] <- 'INVALID';files$message[at] <- 'Declared reference family/species/release/assembly conflict'
  if(!legacy && all(files$status[files$role %in% c('genome_fasta','genome_fai')]=='OK')) {
    message <- tryCatch({setup_genome_check(files$path[files$role=='genome_fasta'],
      files$path[files$role=='genome_fai']);NULL},error=function(e)conditionMessage(e))
    if(!is.null(message)) {
      at <- files$role=='genome_fai';files$status[at]<-'INVALID';files$message[at]<-message
    }
  }
  reference <- if(validate_reference) setup_reference_validation(config) else
    list(status='NOT_RUN',message='Full cross-resource reference validation not requested')
  list(ok=all(files$status=='OK') && reference$status!='FAILED',selected_resources='configured',
    destination=config$paths$data_root,files=files,resources=config$resources,
    scientific=config$scientific,reference_validation=reference,scope='configured reference roles only')
}
