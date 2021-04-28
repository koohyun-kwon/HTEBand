test_that("check positive part", {

  n <- 10
  mat <- cbind(rnorm(n), rnorm(n))
  k <- 2
  res.mat <- mat_sq(mat, arr.ret = FALSE)
  res.arr <- mat_sq(mat, arr.ret = TRUE)

  expect_equal(sum(res.mat[, 1] >= 0), n)
  expect_equal(sum(res.mat[, 4] >= 0), n)
  expect_equal(sum(res.arr[, 1, 1] >= 0), n)
  expect_equal(sum(res.arr[, 2, 2] >= 0), n)
})
