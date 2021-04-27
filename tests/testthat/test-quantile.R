test_that("positive average variance", {
  res <- avar(rep(1/50, 50), rep(1/50, 50), rep(1, 50), rep(1, 50), 1)
  expect_equal(res > 0, TRUE)

  w.1 <- w.0 <- matrix(rep(1/50, 100), ncol = 2)
  omega.1 <- omega.0 <- array(rep(c(1,0,0,1), each = 50), dim = c(50, 2, 2))
  T.grad <- c(1, 1)
  res <- avar(w.1, w.0, omega.1, omega.0, T.grad)
  expect_equal(res > 0, TRUE)
})
