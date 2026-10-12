# Serialization is independent of R's serialize(), JSON formatting and locale.
# Numeric payloads use IEEE-754 binary64, big endian, 16 lowercase hex digits.
canonical_bytes <- function(x) {
  net <- function(tag, payload) paste0(tag, nchar(payload, type = "bytes"), ":", payload)
  if (is.null(x) || (is.atomic(x) && length(x) == 1L && is.na(x) && !is.nan(x)))
    return("N0:")
  if (is.list(x)) {
    if (is.null(names(x))) return(net("A", paste0(vapply(x, canonical_bytes, ""), collapse = "")))
    nm <- stringi::stri_trans_nfc(enc2utf8(names(x)))
    if (anyNA(nm) || anyDuplicated(nm) || any(!nzchar(nm))) stop("Invalid object keys")
    # Hex encodes UTF-8 bytes so radix sorting never depends on locale.
    hex <- vapply(nm, function(s) paste0(sprintf("%02x", as.integer(charToRaw(s))), collapse = ""), "")
    i <- order(hex, method = "radix")
    return(net("O", paste0(vapply(i, function(j)
      paste0(canonical_bytes(nm[j]), canonical_bytes(x[[j]])), ""), collapse = "")))
  }
  if (length(x) != 1L) stop("Canonical scalars require length one; arrays require lists")
  if (is.character(x)) {
    s <- stringi::stri_trans_nfc(enc2utf8(x))
    if (is.na(s)) stop("Invalid UTF-8")
    return(net("S", s))
  }
  if (is.logical(x)) return(net("B", if (x) "1" else "0"))
  if (is.integer(x)) return(net("I", as.character(x)))
  if (is.double(x)) {
    if (!is.finite(x)) stop("Non-finite number")
    if (x == 0) x <- 0 # normalize negative zero
    return(net("D", paste0(sprintf("%02x", as.integer(writeBin(x, raw(), size = 8L,
                                                               endian = "big"))), collapse = "")))
  }
  stop("Unsupported canonical type")
}

biological_id <- function(kind, values, contract = schema_contract()) {
  def <- contract$identities[[kind]]
  if (is.null(def)) stop("Unknown identity kind")
  validate_record(values, def, contract, paste0("identity.", kind))
  ordered <- lapply(def$fields, function(f) {
    v <- values[[f$name]]
    if (!is.null(v) && f$type == "integer") v <- as.integer(v)
    v
  })
  payload <- canonical_bytes(c(list("Neoantimon2", contract$id_schema_version, kind), ordered))
  paste0(kind, ":", digest::digest(charToRaw(payload), algo = "sha256", serialize = FALSE))
}

# Explicit scale conversion only; no inference from the magnitude of a value.
rank_percent <- function(value, input_unit) {
  if (missing(input_unit) || !input_unit %in% c("percent_0_100", "fraction_0_1"))
    stop("Declare input rank unit")
  lim <- if (input_unit == "fraction_0_1") 1 else 100
  if (!is.numeric(value) || any(!is.finite(value[!is.na(value)])) ||
      any(value < 0 | value > lim, na.rm = TRUE)) stop("Rank out of range")
  if (input_unit == "fraction_0_1") value * 100 else value
}
