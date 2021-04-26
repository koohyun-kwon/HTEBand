test_that("valid weights and positive bandwidths", {
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
