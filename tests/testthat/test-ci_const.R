test_that("valid confidence interval - Lipschitz regression", {

  n <- 250
  x.1 <- x.0 <- seq(-1, 1, length.out = n)
  sd.true <- 1/2 + x.1^2
  eps <- stats::rnorm(n, 0, sd.true)
  y.1 <- x.1 + eps
  y.0 <- x.1^2 + eps/2

  y <- c(y.1, y.0)
  x <- c(x.1, x.0)
  d <- rep(c(1, 0), each = n)

  res <- ci_reg_Lip(y, x, 0, 2, 0.99, TE = TRUE, d = d, kern = "tri",
                    bw.eq = FALSE, se.method = "resid")

  res.eq <- ci_reg_Lip(y, x, 0, 2, 0.99, TE = TRUE, d = d, kern = "tri",
                       bw.eq = TRUE, se.method = "resid")

  res
  res.eq

  expect_equal(!is.na(res), rep(TRUE, 4))
  expect_equal(!is.na(res.eq), rep(TRUE, 4))
})


test_that("valid confidence interval - Hölder regression", {

  n <- 250
  x.1 <- x.0 <- seq(-1, 1, length.out = n)
  sd.true <- 1/2 + x.1^2
  eps <- stats::rnorm(n, 0, sd.true)
  y.1 <- x.1 + eps
  y.0 <- x.1^2 + eps/2

  y <- c(y.1, y.0)
  x <- c(x.1, x.0)
  d <- rep(c(1, 0), each = n)

  res.eq <- ci_reg_Hol(y, x, 0, 2, 0.99, TE = TRUE, d = d,
                       bw.eq = TRUE)
  res.eq

  expect_equal(!is.na(res.eq), rep(TRUE, 4))
})
