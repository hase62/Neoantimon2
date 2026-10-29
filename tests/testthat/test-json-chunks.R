test_that("chunked JSON preserves nested jsonlite value semantics", {
  values <- list(NULL, list(), structure(list(), names=character()),
    list(a=1L, 2, absent=NULL, missing=NA_real_, zero=0, fraction=0.123456789012345),
    list(text="quote \" slash \\ newline\n日本語", logical=c(TRUE,FALSE,NA),
      numeric=c(0,NA,Inf,-Inf,NaN), one=list("x"), empty=character()),
    setNames(list(1,2,3),c(NA,"","x")),
    list(frame=data.frame(x=c(1,NA),y=c("a","b")), matrix=matrix(1:4,2),
      date=as.Date("2026-10-30"), raw=as.raw(0:3)),
    list(events=rep(list(list(id="event", evidence=list(vaf=.2, counts=NULL))), 10),
      origins=list(list(source_event_ids=list("event")))))
  for (value in values) {
    expected <- jsonlite::fromJSON(jsonlite::toJSON(value,auto_unbox=TRUE,
      null="null",digits=NA,pretty=TRUE),simplifyVector=FALSE)
    f <- tempfile();con <- file(f,"wb")
    Neoantimon2:::json_write_chunks(value,con,leaf_bytes=1)
    close(con)
    expect_identical(jsonlite::read_json(f),expected)
    unlink(f)
  }
})

test_that("large atomic JSON retains validation and old destination on failure", {
  f <- tempfile();on.exit(unlink(f))
  Neoantimon2:::atomic_json(list(old=TRUE),f)
  # Repeated strings keep the fixture cheap while object.size selects streaming.
  value <- rep(list(list(payload=paste(rep("x",8192),collapse=""), vaf=.2)),9000)
  expect_gt(as.numeric(object.size(value)),64*1024^2)
  expect_error(Neoantimon2:::atomic_json(value,f,before_rename=function(...)stop("injected")),"injected")
  expect_identical(jsonlite::read_json(f),list(old=TRUE))
  Neoantimon2:::atomic_json(value,f)
  expect_identical(jsonlite::read_json(f),value)
})
