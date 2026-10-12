fixture <- function(name) jsonlite::read_json(testthat::test_path("..", "fixtures", "truth", paste0(name, ".json")))
contract <- Neoantimon2:::schema_contract()
valid <- function(x) Neoantimon2:::validate_snapshot(x, contract)
id <- function(kind, values) Neoantimon2:::biological_id(kind, values, contract)
