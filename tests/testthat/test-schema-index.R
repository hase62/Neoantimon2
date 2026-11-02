test_that("schema indexes preserve exact lookup multiplicity and record order", {
  rows <- list(list(key="b", value=1),list(key="a", value=2),
    list(key="b", value=3),list(key="é", value=4),list(key="", value=5))
  index <- Neoantimon2:::schema_row_index(rows,"key")
  for (value in list("a","b","é","missing","",NULL,NA_character_,1L)) {
    expected <- Filter(function(r)identical(r$key,value),rows)
    expect_identical(Neoantimon2:::schema_index_rows(index,value),expected)
  }
  expect_identical(Neoantimon2:::schema_index_rows(
    Neoantimon2:::schema_row_index(list(),"key"),"missing"),list())
  expect_identical(rows[[1L]]$value,1)
})

test_that("indexed validation agrees with independent scan lookups on truth fixtures", {
  cases <- names(fixture("expected")$cases)
  outcome <- function(x)tryCatch(valid(x),error=function(e)conditionMessage(e))
  for (name in cases) {
    x <- fixture(name)
    indexed <- outcome(x)
    scanned <- testthat::with_mocked_bindings(outcome(x),
      schema_index_rows=function(index,value)Filter(
        function(row)identical(row[[index$field]],value),index$rows),
      .package="Neoantimon2")
    expect_identical(indexed,scanned,info=name)
    expect_true(indexed,info=name)
  }
})

test_that("indexed validation retains foreign-key and all-Origin rejection", {
  missing <- fixture("simple_somatic")
  missing$peptide_origins[[1L]]$source_event_ids <- list("nonexistent")
  duplicate <- fixture("simple_somatic")
  duplicate$events[[2L]] <- duplicate$events[[1L]]
  lost <- fixture("multiple_origins")
  lost$neoantigen_results[[1L]]$origin_ids <- lost$neoantigen_results[[1L]]$origin_ids[1L]
  outcome <- function(x)tryCatch(valid(x),error=function(e)conditionMessage(e))
  for (x in list(missing,duplicate,lost)) {
    indexed <- outcome(x)
    scanned <- testthat::with_mocked_bindings(outcome(x),
      schema_index_rows=function(index,value)Filter(
        function(row)identical(row[[index$field]],value),index$rows),
      .package="Neoantimon2")
    expect_type(indexed,"character")
    expect_identical(indexed,scanned)
  }
  expect_error(valid(missing),"foreign key")
  expect_error(valid(duplicate),"duplicate primary key")
  expect_error(valid(lost),"All origins must be preserved")
})
