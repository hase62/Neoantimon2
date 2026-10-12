test_that("installation retains the original required licensing notice", {
  notice <- system.file("NOTICE", package = "Neoantimon2")
  expect_true(nzchar(notice) && file.exists(notice))
  expect_identical(digest::digest(file = notice, algo = "sha256", serialize = FALSE),
    "cf8adbaa2510fdfe9e157fcae9d2095676970481e5a3b025de5181f0aaec07f4")
})
