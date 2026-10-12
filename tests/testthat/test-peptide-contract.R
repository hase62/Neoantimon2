test_that("legacy AA20 and explicit translated AA20+U contracts coexist", {
  old <- Neoantimon2:::schema_contract()
  new <- Neoantimon2:::peptide_contract()
  field <- function(c, table, name) Filter(function(f) f$name == name,
    c$tables[[table]]$fields)[[1L]]
  for (table in c("peptide_origins", "neoantigen_results")) {
    for (name in c("peptide_sequence", "ligand_sequence")) {
      expect_error(Neoantimon2:::validate_value("ACUDE", field(old, table, name), old, name), "pattern")
      expect_silent(Neoantimon2:::validate_value("ACUDE", field(new, table, name), new, name))
      expect_silent(Neoantimon2:::validate_value("ACDEFGHIKLMNPQRSTVWY", field(new, table, name), new, name))
      for (bad in c("X", "B", "Z", "J", "O", "*", "u", "A-C", ""))
        expect_error(Neoantimon2:::validate_value(bad, field(new, table, name), new, name), "pattern")
    }
  }
  for (name in c("wt_peptide_sequence", "upstream_flank", "downstream_flank"))
    expect_silent(Neoantimon2:::validate_value("U", field(new, "peptide_origins", name), new, name))
  expect_silent(Neoantimon2:::validate_value("ACUDE", field(new, "neoantigen_results", "binding_core"), new, "core"))
  core <- Filter(function(f) f$name == "binding_core", new$records$epitope_assignment$fields)[[1L]]
  expect_silent(Neoantimon2:::validate_value("ACUDE", core, new, "core"))
  expect_identical(new$identities, old$identities)
  expect_identical(new$id_schema_version, old$id_schema_version)
  for (sequence in c("ACDE", "ACUE")) {
    values <- list(ligand_sequence = sequence)
    expect_identical(Neoantimon2:::biological_id("peptide_key", values, old),
      Neoantimon2:::biological_id("peptide_key", values, new))
  }
  id <- function(s) Neoantimon2:::biological_id("peptide_key", list(ligand_sequence = s))
  expect_identical(id("ACUE"), id("ACUE"))
  expect_false(id("ACUE") == id("ACCE"))
  expect_false(id("ACUE") == id("ACXE"))
  fixture <- jsonlite::read_json(testthat::test_path("..","fixtures","truth","simple_somatic.json"))
  origin <- fixture$peptide_origins[[1L]]; candidate <- fixture$neoantigen_results[[1L]]
  expect_identical(Neoantimon2:::biological_id("origin_id",origin[c("sample_peptide_id","origin_identity")],new),origin$origin_id)
  expect_identical(Neoantimon2:::biological_id("candidate_key",candidate[c("species","hla_class","ligand_sequence","mhc_key")],new),candidate$candidate_key)
  expect_identical(Neoantimon2:::runtime_contract(fixture)$schema_version,"1.0.0")
  expect_null(Neoantimon2:::runtime_contract(fixture)$peptide_version)
})

test_that("Sec provenance survives the downstream translation envelope", {
  ct <- list(reconstruction_id = "reconstruction:fixture", translation = list(
    translation_id = "translation:fixture", translation_system = "nuclear_sec",
    translation_exception_ids = "annotated-sec:8", sec_sites_used = "annotated-sec:8"))
  p <- Neoantimon2:::peptide_translation_provenance(ct)
  expect_identical(p$translation_id, ct$translation$translation_id)
  expect_identical(p$translation_exception_ids, list("annotated-sec:8"))
  expect_identical(p$sec_sites_used, list("annotated-sec:8"))
})
