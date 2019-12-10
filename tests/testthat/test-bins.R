context("bins")

tmin <- rnorm(10, 0, 5)

test_that("bins works with tmin/tmax vectors", {
  expect_equal(nrow(bins_ss(0, 5, tmin, tmin + 5)), 10)
})

test_that("bins works with tmin/tmax vectors and parrallel = TRUE", {
  expect_equal(nrow(bins_ss(0, 5, tmin, tmin + 5, parallel = TRUE)), 10)
})
