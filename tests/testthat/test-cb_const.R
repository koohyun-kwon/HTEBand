test_that("valid confidence band", {

  n <- 500
  M <- 500
  x <- stats::runif(n, min = -1, max = 1)
  n.T <- 10
  eval <- seq(from = -0.9, to = 0.9, length.out = n.T)
  T.grad.mat <- rep(1, n.T)
  level <- 0.95
  y <- x + rnorm(n, 0, 1/4)

  res <- NA

  system.time({
    opt.res <- cb_const("reg.Hol", 1, y, x, 0, eval, T.grad.mat, level,
                        1, "triangle", FALSE, M, seed = NULL, useloop = TRUE,
                        root.robust = TRUE)
    res <- opt.res$cb.data
    res.inc <- opt.res$increasing
  })

  expect_equal(is.double(res[, 2]), TRUE)
  expect_equal(res.inc, TRUE)
})
