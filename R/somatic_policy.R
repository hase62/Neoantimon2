somatic_evidence_version <- function() "1.0.0"
somatic_eligibility_version <- function() "somatic-technical-v1"
somatic_ccf_model_version <- function() "conditional-multiplicity-enumeration-v1"

somatic_evidence_policy <- function(x) {
  if (is.null(x)) return(NULL)
  strict_keys(x,c("sources","eligibility"),"somatic evidence")
  sources <- if (is.null(x$sources)) list() else x$sources
  if (!is.list(sources) || !is.null(names(sources)))
    stop("somatic evidence sources must be a record array")
  kinds <- c("dna_reads","rna_reads","expression","purity","facets_cn","ccf","loh")
  for (i in seq_along(sources)) {
    s <- sources[[i]]
    strict_keys(s,c("kind","sample_id","path","source","source_version",
      "reference_build","checksum"),"somatic evidence source")
    if (is.null(s$kind) || !s$kind %in% kinds) stop("Unknown somatic evidence source kind")
    scalar_text(s$sample_id,"evidence sample_id")
    scalar_text(s$path,"evidence path")
    if (is.null(s$source)) s$source <- s$kind
    scalar_text(s$source,"evidence source")
    for (key in c("source_version","reference_build","checksum"))
      if (!is.null(s[[key]])) scalar_text(s[[key]],key)
    if (!is.null(s$checksum) && !grepl("^[a-fA-F0-9]{64}$",s$checksum))
      stop("Evidence checksum must be SHA256 hex")
    if (s$kind=="facets_cn" && is.null(s$reference_build))
      stop("FACETS source requires reference_build")
    sources[[i]] <- s
  }
  keys <- vapply(sources,function(s) paste(s$kind,s$sample_id,sep="\037"),"")
  if (anyDuplicated(keys)) stop("Duplicate evidence source kind/sample")
  if (length(sources)) sources <- sources[order(keys,method="radix")]
  policy <- if (is.null(x$eligibility)) list() else x$eligibility
  strict_keys(policy,c("tumor_min_alt_reads","normal_max_alt_reads",
    "rna_min_alt_reads"),"somatic eligibility")
  defaults <- list(tumor_min_alt_reads=5L,normal_max_alt_reads=1L,rna_min_alt_reads=1L)
  for (name in names(defaults)) {
    value <- if (is.null(policy[[name]])) defaults[[name]] else policy[[name]]
    if (!is.numeric(value) || length(value)!=1L || !is.finite(value) ||
        value<0 || value!=floor(value)) stop("Invalid somatic eligibility threshold: ",name)
    defaults[[name]] <- as.integer(value)
  }
  list(sources=sources,eligibility=defaults)
}
