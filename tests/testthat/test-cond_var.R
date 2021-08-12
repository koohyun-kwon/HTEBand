test_that("residual calculation works - no TE", {

  n <- 500
  x <- seq(from = -1, to = 1, length.out = n)
  y <- x^2 + stats::rnorm(n, 0, 0.1)
  expect_equal(sum(resid_calc(y, x, deg = 1)$resid.1^2 > 0), n)

})

test_that("residual calculation works - TE", {

  n <- 500
  x.1 <- seq(from = -1, to = 1, length.out = n/2)
  x.0 <- seq(from = -1, to = 1, length.out = n/2)
  y.1 <- x.1^2 + stats::rnorm(n/2, 0, 0.1)
  y.0 <- x.0 + stats::rnorm(n/2, 0, 0.1)
  d <- rep(c(1, 0), each = n/2)
  y <- c(y.1, y.0)
  x <- c(x.1, x.0)
  resid.res <- resid_calc(y, x, d, deg = 1)
  expect_equal(sum(resid.res$resid.1^2 > 0), n/2)
  expect_equal(sum(resid.res$resid.0^2 > 0), n/2)

})


test_that("residual calculation works - large n", {

  n <- 20000
  x <- seq(from = -1, to = 1, length.out = n)
  y <- x^2 + stats::rnorm(n, 0, 0.5)
  system.time({
    expect_equal(sum(resid_calc(y, x, deg = 1)$resid.1^2 > 0), n)
  })
})
