phase_fixture <- function(name) normalizePath(testthat::test_path("..", "fixtures", "phase", name))

phase_write_vcf <- function(records, source = "Synthetic") {
  path <- tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.2", paste0("##source=", source),
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##FORMAT=<ID=PS,Number=1,Type=String,Description=\"Phase set\">",
    "##FORMAT=<ID=PID,Number=1,Type=String,Description=\"Physical block\">",
    "##FORMAT=<ID=PGT,Number=1,Type=String,Description=\"Physical genotype\">",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
      "INFO", "FORMAT", "sample-1"), collapse = "\t"), records), path)
  path
}

phase_record <- function(pos, gt, ps = NULL, chrom = "chr1", alt = "C", id = ".",
                         pid = NULL, pgt = NULL) {
  ref <- if (chrom == "chrM") c("A", "C", "G", "T")[(pos - 1L) %% 4L + 1L] else
    c("A", "C", "G", "T")[(pos - 1L) %% 4L + 1L]
  fmt <- c(GT = gt)
  if (!is.null(ps)) fmt <- c(fmt, PS = ps)
  if (!is.null(pid)) fmt <- c(fmt, PID = pid)
  if (!is.null(pgt)) fmt <- c(fmt, PGT = pgt)
  paste(c(chrom, pos, id, ref, alt, ".", "PASS", ".",
    paste(names(fmt), collapse = ":"), paste(fmt, collapse = ":")), collapse = "\t")
}

phase_rows <- function(paths, roles = rep("somatic_vcf", length(paths)),
                       ids = paste0("input-", seq_along(paths))) {
  p <- phase_event_test_project(paths, roles, ids)
  Neoantimon2:::phase_resolve_events(Neoantimon2:::event_import_project(p))
}

phase_at <- function(rows, pos, alt = NULL, role = NULL) {
  found <- Filter(function(x) identical(x$pos, as.integer(pos)) &&
    (is.null(alt) || x$alt == alt) && (is.null(role) || x$event_role == role), rows)
  expect_length(found, 1L)
  found[[1L]]
}

phase_relation <- function(rows, a, b) Neoantimon2:::phase_relative_events(rows,
  a$event_id, b$event_id)$relation

test_that("GT phase uses local PS, while slash, no PS, contig and source boundaries hold", {
  file <- phase_write_vcf(c(phase_record(1, "0|1", "44", alt = "C"),
    phase_record(2, "0|1", "44", alt = "T"),
    phase_record(3, "1|0", "44", alt = "A"),
    phase_record(4, "0|1", "45", alt = "C"),
    phase_record(5, "0/1", "44", alt = "C"),
    phase_record(6, "0|1", alt = "T"),
    phase_record(7, "0|1", alt = "A"),
    phase_record(1, "0|1", "44", chrom = "chrM", alt = "G")))
  rows <- phase_rows(file)
  at <- function(p, chrom = "chr1") Filter(function(x) x$pos == p && x$contig == chrom, rows)[[1L]]
  expect_identical(phase_relation(rows, at(1), at(2)), "SAME_HAPLOTYPE")
  expect_identical(phase_relation(rows, at(1), at(3)), "OPPOSITE_HAPLOTYPE")
  for (p in 4:7) expect_identical(phase_relation(rows, at(1), at(p)), "UNKNOWN")
  expect_identical(phase_relation(rows, at(6), at(7)), "UNKNOWN")
  expect_identical(phase_relation(rows, at(1), at(1, "chrM")), "UNKNOWN")
  other <- phase_write_vcf(phase_record(8, "0|1", "44", alt = "C"))
  isolated <- phase_rows(c(file, other), ids = c("a", "b"))
  a <- Filter(function(x) x$pos == 1L && x$contig == "chr1", isolated)[[1L]]
  b <- phase_at(isolated, 8)
  expect_identical(phase_relation(isolated, a, b), "UNKNOWN")
  expect_length(unlist(lapply(isolated, `[[`, "phase_links_json"), recursive = FALSE), 0L)
})

test_that("split ALT indexes and homozygous or absent alleles remain distinct", {
  one <- phase_write_vcf(c(phase_record(1, "1|2", "11", alt = "C,G"),
    phase_record(2, "2|1", "11", alt = "T,G"),
    phase_record(3, "1|1", "11", alt = "A"),
    phase_record(4, "0/1", "11", alt = "C,G")))
  rows <- phase_rows(one)
  c1 <- phase_at(rows, 1, "C"); g1 <- phase_at(rows, 1, "G")
  t2 <- phase_at(rows, 2, "T"); g2 <- phase_at(rows, 2, "G")
  expect_identical(phase_relation(rows, c1, g1), "OPPOSITE_HAPLOTYPE")
  expect_identical(phase_relation(rows, t2, g2), "OPPOSITE_HAPLOTYPE")
  expect_identical(phase_relation(rows, c1, g2), "SAME_HAPLOTYPE")
  expect_identical(phase_at(rows, 3)$event_details_json$phase_assertions[[1]]$genotype_state, "BOTH")
  expect_identical(phase_relation(rows, c1, phase_at(rows, 3)), "UNKNOWN")
  expect_identical(phase_relation(rows, c1, phase_at(rows, 4, "C")), "UNKNOWN")
  expect_identical(phase_at(rows, 4, "G")$event_details_json$phase_assertions[[1]]$genotype_state, "ABSENT")
  for (gt in c(".|.", ".|1", "1|.", "1", "0|1|1", "0|9")) {
    src <- list(format = "VCF", input_id = "input-1", sample_id = "sample-1",
      line = 1L, caller = NULL, original = list(chrom = "chr1", pos = 1L,
        alts = list("C")), selected_alt_index = 1L,
      raw_format = list(keys = "GT:PS", sample = paste0(gt, ":44")))
    row <- phase_at(rows, 1, "C")
    expect_null(Neoantimon2:::phase_source_assertions(row, src)[[1]]$local_side)
  }
})

test_that("PID/PGT can phase slash GT; contradictory evidence is retained", {
  f <- phase_write_vcf(c(phase_record(1, "0/1", pid = "physical-x", pgt = "0|1", alt = "C"),
    phase_record(2, "0/1", pid = "physical-x", pgt = "1|0", alt = "T"),
    phase_record(3, "0|1", ps = "44", pid = "physical-x", pgt = "1|0", alt = "A"),
    phase_record(4, "0/1", pid = "physical-x", pgt = "1|1", alt = "C")))
  rows <- phase_rows(f)
  expect_identical(phase_relation(rows, phase_at(rows, 1), phase_at(rows, 2)),
    "OPPOSITE_HAPLOTYPE")
  expect_identical(phase_at(rows, 3)$phase_status, "NOT_EVALUABLE")
  expect_identical(phase_at(rows, 4)$phase_status, "NOT_EVALUABLE")
  expect_identical(phase_relation(rows, phase_at(rows, 1), phase_at(rows, 3)), "UNKNOWN")
  expect_identical(phase_at(rows, 3)$event_details_json$source_records[[1]]$raw_format$sample,
    "0|1:44:physical-x:1|0")
  bridge <- phase_rows(phase_write_vcf(c(
    phase_record(1, "0|1", ps = "44", pid = "physical-x", pgt = "0|1", alt = "C"),
    phase_record(2, "0|1", ps = "44", alt = "T"),
    phase_record(3, "0/1", pid = "physical-x", pgt = "1|0", alt = "A"))))
  expect_identical(phase_relation(bridge, phase_at(bridge, 2), phase_at(bridge, 3)),
    "OPPOSITE_HAPLOTYPE")
})

test_that("exact shared markers orient sources and conflicting anchors block propagation", {
  a <- phase_fixture("source_a.vcf")
  for (variant in c("source_b_same.vcf", "source_b_flipped.vcf")) {
    rows <- phase_rows(c(a, phase_fixture(variant)), ids = c("a", "b"))
    relation <- phase_relation(rows, phase_at(rows, 1), phase_at(rows, 3))
    expect_identical(relation, if (variant == "source_b_same.vcf")
      "SAME_HAPLOTYPE" else "OPPOSITE_HAPLOTYPE")
    edges <- unlist(lapply(rows, `[[`, "phase_links_json"), recursive = FALSE)
    expect_length(edges, 1L)
    expect_identical(edges[[1]]$orientation, if (variant == "source_b_same.vcf")
      "SAME" else "FLIPPED")
    expect_identical(edges[[1]]$evidence, "SHARED_MARKER")
  }
  a2 <- phase_write_vcf(c(phase_record(1, "0|1", "44", alt = "C"),
    phase_record(2, "0|1", "44", alt = "T"),
    phase_record(4, "0|1", "44", alt = "C")))
  agree <- phase_write_vcf(c(phase_record(2, "0|1", "91", alt = "T"),
    phase_record(4, "0|1", "91", alt = "C"),
    phase_record(3, "0|1", "91", alt = "A")))
  rows <- phase_rows(c(a2, agree), ids = c("a", "b"))
  edge <- unlist(lapply(rows, `[[`, "phase_links_json"), recursive = FALSE)[[1L]]
  expect_identical(edge$orientation, "SAME")
  expect_length(edge$anchors, 2L)
  conflict <- phase_rows(c(a2, phase_fixture("source_b_conflict.vcf")), ids = c("a", "b"))
  expect_identical(phase_relation(conflict, phase_at(conflict, 1), phase_at(conflict, 3)), "UNKNOWN")
  expect_identical(phase_relation(conflict, phase_at(conflict, 1), phase_at(conflict, 2)),
    "SAME_HAPLOTYPE")
  edge <- unlist(lapply(conflict, `[[`, "phase_links_json"), recursive = FALSE)[[1L]]
  expect_identical(edge$orientation, "CONFLICT")
  expect_identical(edge$evidence, "SHARED_MARKER")
  opposite <- phase_rows(c(phase_write_vcf(c(
    phase_record(1, "0|1", "44", alt = "C"),
    phase_record(2, "1|0", "44", alt = "T"),
    phase_record(4, "1|0", "44", alt = "C"))),
    phase_fixture("source_b_conflict.vcf")), ids = c("a", "b"))
  expect_identical(phase_relation(opposite, phase_at(opposite, 1), phase_at(opposite, 2)),
    "OPPOSITE_HAPLOTYPE")
  expect_identical(phase_relation(opposite, phase_at(opposite, 1), phase_at(opposite, 3)),
    "UNKNOWN")
  removed <- phase_rows(c(phase_write_vcf(phase_record(1, "0|1", "44", alt = "C")),
    phase_fixture("source_b_same.vcf")), ids = c("a", "b"))
  expect_identical(phase_relation(removed, phase_at(removed, 1), phase_at(removed, 3)), "UNKNOWN")
})

test_that("compatible same-record source bridge survives cross-source conflict", {
  a <- phase_write_vcf(c(
    phase_record(1, "0|1", ps = "44", pid = "physical-x", pgt = "0|1", alt = "C"),
    phase_record(2, "0|1", ps = "44", alt = "T"),
    phase_record(3, "0/1", pid = "physical-x", pgt = "1|0", alt = "A"),
    phase_record(4, "0|1", ps = "44", alt = "C")))
  b <- phase_write_vcf(c(
    phase_record(1, "0|1", ps = "91", alt = "C"),
    phase_record(4, "1|0", ps = "91", alt = "C"),
    phase_record(5, "0|1", ps = "91", alt = "C")))
  rows <- phase_rows(c(a, b), ids = c("a", "b"))
  edges <- unlist(lapply(rows, `[[`, "phase_links_json"), recursive = FALSE)
  expect_true(any(vapply(edges, function(x) x$evidence == "SOURCE_ASSERTION" &&
    x$orientation == "SAME", TRUE)))
  expect_true(any(vapply(edges, function(x) x$evidence == "SHARED_MARKER" &&
    x$orientation == "CONFLICT", TRUE)))
  expect_identical(Neoantimon2:::phase_relative_events(rows,
    phase_at(rows, 2)$event_id, phase_at(rows, 3)$event_id),
    list(relation = "OPPOSITE_HAPLOTYPE", evidence = "SOURCE_ASSERTION"))
  expect_identical(phase_relation(rows, phase_at(rows, 2), phase_at(rows, 5)), "UNKNOWN")
})

test_that("somatic and germline PS collision is isolated; shared marker can orient them", {
  isolated <- phase_rows(c(phase_event_fixture("phase_somatic.vcf"),
    phase_event_fixture("phase_germline.vcf")),
    roles = c("somatic_vcf", "germline_vcf"), ids = c("somatic", "germline"))
  expect_identical(phase_relation(isolated, phase_at(isolated, 1, role = "somatic"),
    phase_at(isolated, 3, role = "germline")), "UNKNOWN")
  expect_length(unlist(lapply(isolated, `[[`, "phase_links_json"), recursive = FALSE), 0L)
  linked <- phase_rows(c(phase_fixture("source_a.vcf"), phase_fixture("source_b_same.vcf")),
    roles = c("somatic_vcf", "germline_vcf"), ids = c("somatic", "germline"))
  expect_identical(phase_relation(linked, phase_at(linked, 1, role = "somatic"),
    phase_at(linked, 3, role = "germline")), "SAME_HAPLOTYPE")
  # A coincident PS value in B changes no cross-source relation.
  b44 <- phase_write_vcf(c(phase_record(2, "0|1", "44", alt = "T"),
    phase_record(3, "0|1", "44", alt = "A")))
  still_linked <- phase_rows(c(phase_fixture("source_a.vcf"), b44),
    roles = c("somatic_vcf", "germline_vcf"), ids = c("somatic", "germline"))
  expect_identical(phase_relation(still_linked, phase_at(still_linked, 1, role = "somatic"),
    phase_at(still_linked, 3, role = "germline")), "SAME_HAPLOTYPE")
})

test_that("phase semantics are deterministic and IDs are unchanged", {
  records <- c(phase_record(1, "0|1", "44", alt = "C"),
    phase_record(2, "1|0", "44", alt = "T"),
    phase_record(3, "0|1", "44", alt = "A"))
  f1 <- phase_write_vcf(records); f2 <- phase_write_vcf(rev(records))
  p1 <- phase_event_test_project(f1, "somatic_vcf", "fixed")
  p2 <- phase_event_test_project(f2, "somatic_vcf", "fixed")
  raw1 <- Neoantimon2:::event_import_project(p1)
  raw2 <- Neoantimon2:::event_import_project(p2)
  a <- Neoantimon2:::phase_resolve_events(raw1)
  b <- Neoantimon2:::phase_resolve_events(raw2)
  expect_identical(vapply(a, `[[`, "", "event_id"), vapply(raw1, `[[`, "", "event_id"))
  expect_identical(vapply(a, `[[`, "", "event_key"), vapply(raw1, `[[`, "", "event_key"))
  signature <- function(rows) lapply(rows, function(x) lapply(x$event_details_json$phase_assertions,
    function(z) z[c("block_id", "local_side", "genotype_state")]))
  expect_identical(signature(a), signature(b))
  swapped <- phase_rows(c(phase_fixture("source_a.vcf"), phase_fixture("source_b_same.vcf")),
    ids = c("a", "b"))
  reversed <- phase_rows(c(phase_fixture("source_b_same.vcf"), phase_fixture("source_a.vcf")),
    ids = c("b", "a"))
  expect_identical(lapply(swapped, `[[`, "phase_links_json"),
    lapply(reversed, `[[`, "phase_links_json"))
  expect_identical(phase_relation(swapped, phase_at(swapped, 1), phase_at(swapped, 3)),
    phase_relation(reversed, phase_at(reversed, 1), phase_at(reversed, 3)))
})

test_that("one large block stores linear phase metadata and zero rows succeed", {
  records <- vapply(1:60, function(p) phase_record(p, "0|1", "44",
    alt = c("C", "T", "A", "C")[(p - 1L) %% 4L + 1L]), "")
  rows <- phase_rows(phase_write_vcf(records))
  expect_length(rows, 60L)
  expect_length(unlist(lapply(rows, `[[`, "phase_links_json"), recursive = FALSE), 0L)
  expect_identical(Neoantimon2:::phase_resolve_events(list()), list())
  bytes <- nchar(as.character(jsonlite::toJSON(lapply(rows, function(x)
    x$event_details_json$phase_assertions), auto_unbox = TRUE)), type = "bytes")
  expect_lt(bytes, 60000)
})

test_that("duplicate caller assertions remain independent and disagreements are unknown", {
  a <- phase_write_vcf(c(phase_record(1, "0|1", "44", alt = "C"),
    phase_record(2, "0|1", "44", alt = "T")), "CallerA")
  b <- phase_write_vcf(c(phase_record(1, "0|1", "88", alt = "C"),
    phase_record(2, "0|1", "88", alt = "T")), "CallerB")
  compatible <- phase_rows(c(a, b), ids = c("a", "b"))
  expect_length(phase_at(compatible, 1)$event_details_json$phase_assertions, 2L)
  expect_identical(phase_relation(compatible, phase_at(compatible, 1),
    phase_at(compatible, 2)), "SAME_HAPLOTYPE")
  reversed <- phase_write_vcf(c(phase_record(1, "0|1", "88", alt = "C"),
    phase_record(2, "1|0", "88", alt = "T")), "CallerB")
  incompatible <- phase_rows(c(a, reversed), ids = c("a", "b"))
  expect_length(phase_at(incompatible, 1)$event_details_json$phase_assertions, 2L)
  expect_identical(phase_relation(incompatible, phase_at(incompatible, 1),
    phase_at(incompatible, 2)), "UNKNOWN")
  contradictory_anchor <- phase_write_vcf(c(
    phase_record(1, "0|1", "44", alt = "C"),
    phase_record(2, "0|1", "44", alt = "T"),
    phase_record(2, "1|1", "44", alt = "T")))
  other <- phase_write_vcf(c(phase_record(2, "0|1", "88", alt = "T"),
    phase_record(3, "0|1", "88", alt = "A")))
  ambiguous <- phase_rows(c(contradictory_anchor, other), ids = c("a", "b"))
  expect_identical(phase_relation(ambiguous, phase_at(ambiguous, 1),
    phase_at(ambiguous, 3)), "UNKNOWN")
  expect_length(unlist(lapply(ambiguous, `[[`, "phase_links_json"), recursive = FALSE), 0L)
})

test_that("PS changes, local flips and filename changes obey relative invariants", {
  base <- c(phase_record(1, "0|1", "44", alt = "C"),
    phase_record(2, "1|0", "44", alt = "T"))
  original <- phase_rows(phase_write_vcf(base), ids = "same-id")
  changed <- phase_rows(phase_write_vcf(c(base[[1]],
    phase_record(2, "1|0", "45", alt = "T"))), ids = "same-id")
  flipped <- phase_rows(phase_write_vcf(c(phase_record(1, "1|0", "44", alt = "C"),
    phase_record(2, "0|1", "44", alt = "T"))), ids = "same-id")
  expect_identical(phase_relation(original, phase_at(original, 1), phase_at(original, 2)),
    "OPPOSITE_HAPLOTYPE")
  expect_identical(phase_relation(changed, phase_at(changed, 1), phase_at(changed, 2)),
    "UNKNOWN")
  expect_identical(phase_relation(flipped, phase_at(flipped, 1), phase_at(flipped, 2)),
    "OPPOSITE_HAPLOTYPE")
  copy <- tempfile(fileext = ".vcf")
  file.copy(phase_fixture("source_a.vcf"), copy)
  original <- phase_rows(phase_fixture("source_a.vcf"), ids = "stable-input")
  renamed <- phase_rows(copy, ids = "stable-input")
  expect_identical(lapply(original, function(x) x$event_details_json$phase_assertions[[1]]$block_id),
    lapply(renamed, function(x) x$event_details_json$phase_assertions[[1]]$block_id))
})

test_that("transitive orientation is labeled derived and impossible cycles are conflicts", {
  a <- phase_write_vcf(c(phase_record(1, "0|1", "11", alt = "C"),
    phase_record(2, "0|1", "11", alt = "T"),
    phase_record(4, "0|1", "11", alt = "C")))
  b <- phase_write_vcf(c(phase_record(2, "0|1", "22", alt = "T"),
    phase_record(3, "0|1", "22", alt = "A"),
    phase_record(6, "0|1", "22", alt = "T")))
  c <- phase_write_vcf(c(phase_record(3, "1|0", "33", alt = "A"),
    phase_record(5, "0|1", "33", alt = "C")))
  rows <- phase_rows(c(a, b, c), ids = c("a", "b", "c"))
  derived <- Neoantimon2:::phase_relative_events(rows, phase_at(rows, 1)$event_id,
    phase_at(rows, 5)$event_id)
  expect_identical(derived, list(relation = "OPPOSITE_HAPLOTYPE", evidence = "TRANSITIVE"))
  cycle <- phase_write_vcf(c(phase_record(3, "1|0", "33", alt = "A"),
    phase_record(4, "0|1", "33", alt = "C"),
    phase_record(5, "0|1", "33", alt = "C"),
    phase_record(6, "1|0", "33", alt = "T")))
  bad <- phase_rows(c(a, b, cycle), ids = c("a", "b", "c"))
  expect_identical(phase_relation(bad, phase_at(bad, 1), phase_at(bad, 5)), "UNKNOWN")
  expect_identical(phase_relation(bad, phase_at(bad, 1), phase_at(bad, 2)),
    "SAME_HAPLOTYPE")
  expect_identical(phase_relation(bad, phase_at(bad, 3), phase_at(bad, 5)),
    "OPPOSITE_HAPLOTYPE")
  edges <- unlist(lapply(bad, `[[`, "phase_links_json"), recursive = FALSE)
  expect_true(all(vapply(edges, function(x) x$orientation == "CONFLICT", TRUE)))
})

test_that("site-only germline and zero-event runtime complete without phase links", {
  site <- phase_write_vcf(phase_record(1, "0|1", "44", alt = "C"))
  # A site-only source has no FORMAT/sample columns.
  writeLines(c("##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t1\t.\tA\tC\t.\tPASS\t."), site)
  p <- phase_event_test_project(site, "germline_vcf", "site")
  p <- neo_run(p)
  row <- Neoantimon2:::read_runtime_snapshot(p)$events[[1L]]
  expect_identical(row$phase_status, "NOT_PROVIDED")
  expect_length(row$phase_links_json, 0L)
  expect_true("phase_resolution" %in% vapply(Neoantimon2:::read_manifest(p$path)$checkpoints,
    `[[`, "", "step_id"))
  writeLines("tampered", file.path(p$path, "intermediate", "phase_resolution.json"))
  expect_error(Neoantimon2:::read_runtime_snapshot(p), "checkpoint|mismatch|checksum")
  empty <- phase_write_vcf(character())
  p0 <- neo_run(phase_event_test_project(empty, "germline_vcf", "empty"))
  expect_length(Neoantimon2:::read_runtime_snapshot(p0)$events, 0L)
  expect_identical(Filter(function(x) x$step_id == "phase_resolution",
    Neoantimon2:::read_manifest(p0$path)$checkpoints)[[1L]]$status, "COMPLETED")
})
