# Setup is opt-in and independent of neo_project()/neo_run().
setup_flag <- function(x, label) {
  if (!is.logical(x) || length(x)!=1L || is.na(x)) stop(label," must be TRUE or FALSE")
  x
}
setup_catalog <- function() {
  path <- system.file('setup','resource-catalog.json',package='Neoantimon2')
  if (!nzchar(path)) stop('Installed resource catalog is missing; reinstall this development build')
  x <- setup_validate_catalog(jsonlite::read_json(path))
  attr(x,'sha256') <- digest::digest(file=path,algo='sha256',serialize=FALSE)
  x
}
setup_validate_catalog <- function(x) {
  if (!identical(x$schema_version,'neo-resource-catalog-1') ||
      !is.list(x$resources) || !length(x$resources) || is.null(names(x$resources)) ||
      anyDuplicated(names(x$resources)) || any(!nzchar(names(x$resources)))) stop('Invalid resource catalog')
  paths <- character()
  for (spec in x$resources) for (f in spec$files) {
    scalar_text(f$path,'catalog path')
    parts <- strsplit(chartr('\\','/',f$path),'/',fixed=TRUE)[[1L]]
    if (grepl('^(/|[A-Za-z]:)',f$path) || grepl('/$',f$path) || any(parts %in% c('','..','.')) ||
        !f$action %in% c('download','gunzip','faidx','sra') ||
        !f$kind %in% c('fasta','transcript','protein','gtf','gencode-gtf','sec','mane-summary','fai','vcf','tbi','bam','bai','bed','sra','text'))
      stop('Invalid catalog path/action/kind')
    if (!is.null(f$url) && (!grepl('^(https?|file)://',f$url) ||
        grepl('[@\r\n?]',f$url))) stop('Catalog URL must be a public URL without credentials')
    if (!is.null(f$sha256) && !grepl('^[a-f0-9]{64}$',f$sha256)) stop('Invalid catalog checksum')
    paths <- c(paths,f$path)
  }
  if (anyDuplicated(paths) || any(!unlist(x$profiles,use.names=FALSE) %in% names(x$resources)))
    stop('Invalid catalog selection/path mapping')
  if(!setequal(unlist(x$profiles$full,use.names=FALSE),names(x$resources)) ||
      !identical(unlist(x$profiles$minimum,use.names=FALSE),c('genome','gencode','mane')))
    stop('Invalid catalog profile mapping')
  x
}
setup_selection <- function(catalog,profile,resources) {
  if (is.null(profile) && is.null(resources)) stop('Choose profile or resources explicitly')
  out <- character()
  if (!is.null(profile)) {
    scalar_text(profile,'profile')
    if (!profile %in% names(catalog$profiles)) stop('Unknown resource profile: ',profile)
    out <- unlist(catalog$profiles[[profile]],use.names=FALSE)
  }
  if (!is.null(resources)) {
    if (!is.character(resources) || !length(resources) || anyNA(resources) ||
        any(!resources %in% names(catalog$resources))) stop('Unknown resource selection')
    out <- c(out,resources)
  }
  unique(out)
}
setup_dest <- function(path,create=FALSE) {
  scalar_text(path,'destination')
  path <- resolve_path(path,getwd())
  if (file.exists(path) && !dir.exists(path)) stop('Destination is not a directory: ',path)
  if (create && !dir.exists(path) && !dir.create(path,recursive=TRUE,showWarnings=FALSE))
    stop('Cannot create destination: ',path)
  if (create && file.access(path,2L)!=0L) stop('Destination is not writable: ',path)
  path
}
setup_executable <- function(path,default=NULL) {
  if (is.null(path)) path <- if (is.null(default)) '' else unname(Sys.which(default))
  if (!is.character(path) || length(path)!=1L || is.na(path) || !nzchar(path)) return(NULL)
  if (!grepl('[/\\\\]',path)) path <- unname(Sys.which(path))
  if (!nzchar(path) || !file.exists(path) || dir.exists(path) || file.access(path,1L)!=0L) return(NULL)
  # Do not resolve a venv Python leaf symlink to its base interpreter.
  file.path(normalizePath(dirname(path),winslash='/',mustWork=TRUE),basename(path))
}
setup_command <- function(executable,args=character(),timeout=0) {
  out <- tryCatch(suppressWarnings(system2(executable,vapply(args,shQuote,''),
    stdout=TRUE,stderr=TRUE,timeout=timeout)),error=function(e)
      structure(conditionMessage(e),status=1L))
  status <- attr(out,'status')
  list(status=if(is.null(status)) 0L else as.integer(status),output=as.character(out))
}
setup_python <- function(python) {
  path <- setup_executable(python,'python3')
  if (is.null(path)) setup_dependency('Python 3 is required only for this selected operation; supply python= or install it')
  path
}
setup_dependency <- function(message) stop(structure(list(message=message,call=NULL),
  class=c('neo_setup_dependency_error','error','condition')))
setup_python_resource <- function(id,dest,python,check=FALSE,overwrite=FALSE,checksum=TRUE) {
  script <- system.file('setup','download_public_data.py',package='Neoantimon2')
  catalog <- system.file('setup','resource-catalog.json',package='Neoantimon2')
  if (!nzchar(script) || !nzchar(catalog)) stop('Installed acquisition assets missing; reinstall the package')
  py <- setup_python(python)
  tools <- unlist(setup_catalog()$resources[[id]]$tools$python,use.names=FALSE)
  if(check) tools <- setdiff(tools,'curl')
  missing <- tools[!nzchar(Sys.which(tools))]
  if(length(missing)) setup_dependency(paste('Required command(s):',paste(missing,collapse=', ')))
  args <- c(script,'--catalog',catalog,'--dest',dest,'--resource',id,
    if(check) '--check-only',if(!overwrite) '--no-overwrite',if(checksum) '--record-sha256')
  result <- setup_command(py,args)
  if(result$status!=0L) stop('Resource ',id,': ',paste(result$output,collapse='\n'),call.=FALSE)
  invisible(result)
}
setup_sha <- function(path) digest::digest(file=path,algo='sha256',serialize=FALSE)
setup_json <- function(x,path,overwrite=FALSE) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  stage <- tempfile('.neo-metadata-',tmpdir=dirname(path))
  on.exit(unlink(stage),add=TRUE)
  jsonlite::write_json(x,stage,auto_unbox=TRUE,pretty=TRUE,null='null')
  setup_promote(stage,path,overwrite)
}
setup_promote <- function(stage,target,overwrite=FALSE) {
  dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
  if (file.exists(target)) {
    if (!dir.exists(target) && identical(setup_sha(stage),setup_sha(target))) return(invisible(target))
    if (!overwrite) stop('Existing resource retained: ',target,'; inspect it or choose overwrite=TRUE')
    link <- Sys.readlink(target)
    if (dir.exists(target) || (!is.na(link) && nzchar(link))) stop('Refusing to replace a directory/symlink: ',target)
    # Retain the previous bytes until the replacement has been promoted.
    backup <- tempfile('.neo-previous-',tmpdir=dirname(target))
    if(!file.rename(target,backup)) stop('Cannot preserve previous resource: ',target)
    if(!file.rename(stage,target)) {
      file.rename(backup,target)
      stop('Could not promote verified resource: ',target)
    }
    unlink(backup)
  } else if(!file.rename(stage,target)) stop('Could not promote verified resource: ',target)
  invisible(target)
}
setup_write_target <- function(root,relative) {
  path <- file.path(root,relative)
  parent <- resolve_path(dirname(path),root)
  if(!within_path(parent,root)) stop('Target parent escapes destination through a symlink: ',relative)
  path
}
