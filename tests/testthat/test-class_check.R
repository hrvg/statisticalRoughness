context("factors")

test_that("get_factors returns integers", {
  expect_is(get_factors(12), "integer")
})

test_that("get_all_R_L returns a list with 2 elements for all options", {
  expect_is(get_all_R_L(1024, 5, only = NULL, logfilter = TRUE, len = 32), "list")
  expect_is(get_all_R_L(1024, 5, only = 610, logfilter = TRUE, len = 32), "list")
  expect_is(get_all_R_L(1024, 5, only = 12, logfilter = TRUE, len = 32), "list")
  expect_is(get_all_R_L(1024, 5, only = NULL, logfilter = FALSE, len = 32), "list")
  expect_is(get_all_R_L(1024, 5, only = 610, logfilter = FALSE, len = 32), "list")
  expect_is(get_all_R_L(1024, 5, only = 12, logfilter = FALSE, len = 32), "list")
  expect_equal(length(get_all_R_L(1024, 5, only = NULL, logfilter = TRUE, len = 32)), 2)
  expect_equal(length(get_all_R_L(1024, 5, only = 610, logfilter = TRUE, len = 32)), 2)
  expect_equal(length(get_all_R_L(1024, 5, only = 12, logfilter = TRUE, len = 32)), 2)
  expect_equal(length(get_all_R_L(1024, 5, only = NULL, logfilter = FALSE, len = 32)), 2)
  expect_equal(length(get_all_R_L(1024, 5, only = 610, logfilter = FALSE, len = 32)), 2)
  expect_equal(length(get_all_R_L(1024, 5, only = 12, logfilter = FALSE, len = 32)), 2)
})

test_that("get_all_R_L returns the right number of requested factors", {
  expect_equal(length(get_all_R_L(1024, 5, only = NULL, logfilter = TRUE, len = 3)[[1]]), 3)
  expect_equal(length(get_all_R_L(1024, 5, only = 610, logfilter = TRUE, len = 3)[[1]]), 3)
  expect_equal(length(get_all_R_L(1024, 5, only = 12, logfilter = TRUE, len = 3)[[1]]), 3)
})

test_that("get_all_R_L throws an error when not enough factors returned", {
  expect_error(get_all_R_L(24, 3), 'Number of factors return is lower than requested.')
})