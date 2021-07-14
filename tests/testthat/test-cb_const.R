test_that("valid confidence band (Holder)", {

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
    opt.res <- cb_const("reg.Hol", 1, y, x, NULL, eval, T.grad.mat, level,
                        1, "tri", M, seed = NULL, useloop = TRUE,
                        root.robust = FALSE)
  })

  expect_equal(is.na(opt.res[, 2]), rep(FALSE, n.T))

  # system.time({
  #   opt.res <- cb_const("reg.Hol", 1, y, x, NULL, eval, T.grad.mat, level,
  #                       1, "tri", M, seed = NULL, useloop = TRUE,
  #                       root.robust = TRUE)
  #   res <- opt.res$cb.data
  #   res.inc <- opt.res$increasing
  # })
  #
  # expect_equal(is.double(res[, 2]), TRUE)
})


test_that("valid confidence band (Lipschitz)", {

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
    opt.res <- cb_const("reg.Lip", 1, y, x, NULL, eval, T.grad.mat, level,
                        1, "tri", M, seed = NULL, useloop = TRUE,
                        root.robust = FALSE)
  })

  expect_equal(is.na(opt.res[, 2]), rep(FALSE, n.T))

  # system.time({
  #   opt.res <- cb_const("reg.Lip", 1, y, x, NULL, eval, T.grad.mat, level,
  #                       1, "tri", M, seed = NULL, useloop = TRUE,
  #                       root.robust = TRUE)
  #   res <- opt.res$cb.data
  #   res.inc <- opt.res$increasing
  # })
  #
  # expect_equal(is.double(res[, 2]), TRUE)
  # expect_equal(res.inc, TRUE)
})


test_that("valid confidence band (Lipschitz-TE)", {

  n <- 250
  x.1 <- x.0 <- seq(-1, 1, length.out = n)
  sd.true <- 1/2 + x.1^2
  eps <- stats::rnorm(n, 0, sd.true)
  y.1 <- x.1 + eps
  y.0 <- x.1^2 + eps/2

  y <- c(y.1, y.0)
  x <- c(x.1, x.0)
  d <- rep(c(1, 0), each = n)

  M <- 500
  n.T <- 10
  eval <- seq(from = -0.9, to = 0.9, length.out = n.T)

  T.grad.mat <- rep(1, n.T)
  level <- 0.95

  res <- NA

  system.time({
    opt.res <- cb_const("TE.Lip", 2, y, x, d, eval, T.grad.mat, level,
                        1, "tri", M, seed = NULL, useloop = TRUE,
                        root.robust = FALSE)
  })

  expect_equal(is.na(opt.res[, 2]), rep(FALSE, n.T))

  # system.time({
  #   opt.res <- cb_const("TE.Lip", 2, y, x, d, eval, T.grad.mat, level,
  #                       1, "tri", M, seed = NULL, useloop = TRUE,
  #                       root.robust = TRUE)
  #   res <- opt.res$cb.data
  #   res.inc <- opt.res$increasing
  # })
  #
  # expect_equal(is.double(res[, 2]), TRUE)
  # expect_equal(res.inc, TRUE)
})
