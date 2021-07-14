test_that("valid weights and positive bandwidths: Holder", {
  n <- 500
  J <- 5
  x <- stats::runif(n, min = -1, max = 1)
  y <- x + rnorm(n, 0, 1/4)
  eval <- seq(from = -0.9, to = 0.9, length.out = J)
  res <- w_get_Hol(y, x, eval, 1, 0.95)
  w.pos <- is.na(res$w.mat)
  bw.pos <- res$bw.vec > 0
  expect_equal(sum(w.pos), 0)
  expect_equal(sum(bw.pos), J)
})

test_that("valid weights: Lipschitz", {

  n <- 250
  x.1 <- x.0 <- seq(-1, 1, length.out = n)
  sd.true <- 1/2 + x.1^2
  eps <- stats::rnorm(n, 0, sd.true)
  y.1 <- x.1 + eps
  y.0 <- x.1^2 + eps/2
  m = 5
  eval <- seq(from = -0.9, to = 0.9, length.out = m)

  y <- c(y.1, y.0)
  x <- c(x.1, x.0)
  d <- rep(c(1, 0), each = n)

  res <- w_get_Lip(y, x, eval, 2, 0.95, TE = TRUE, d = d, kern = "tri", bw.eq = FALSE,
                   deg = 0)
  res.eq <- w_get_Lip(y, x, eval, 2, 0.95, TE = TRUE, d = d, kern = "tri", bw.eq = TRUE,
                      deg = 0)

  expect_equal(as.numeric(res$w.mat.1 >= 0), rep(1, n * m))
  expect_equal(as.numeric(res$w.mat.0 >= 0), rep(1, n * m))
  expect_equal(as.numeric(res.eq$w.mat.1 >= 0), rep(1, n * m))
  expect_equal(as.numeric(res.eq$w.mat.0 >= 0), rep(1, n * m))

})
