# test_that("valid optimal quantile value (Holder)", {
#
#   n <- 500
#   M <- 500
#   x <- stats::runif(n, min = -1, max = 1)
#   n.T <- 10
#   eval <- seq(from = -0.9, to = 0.9, length.out = n.T)
#   T.grad.mat <- rep(1, n.T)
#   level <- 0.95
#   y <- x + rnorm(n, 0, 1/4)
#
#   res <- NA
#
#   system.time({
#     opt.res <- opt_w("reg.Hol", 1, y, x, NULL, eval, T.grad.mat, level,
#                       1, "tri", M, seed = NULL, useloop = TRUE,
#                       root.robust = TRUE)
#     res <- opt.res$c.root
#     res.inc <- opt.res$increasing
#   })
#
#   expect_equal(!is.na(res), TRUE)
#   expect_equal(res.inc, TRUE)
#
#   # Optional test
#
#   # kern.reg <- "triangular"
#   # se.initial <- "EHW"
#   # se.method <- "nn"
#   # J <- 3
#   # C <- 1
#   #
#   # eq <- function(c){
#   #
#   #   level.int <- stats::pnorm(2 * c)
#   #   w.1 <- w_get_Hol(y, x, eval, C, level.int, kern.reg, se.initial, se.method, J)$w.mat
#   #   w.1 <- array(w.1, dim = c(length(y), 1, n.T))
#   #   w.0 <- array(0, dim = c(1, 1, n.T))
#   #
#   #   q.sim <- sup_quant_sim(y, 0, x, 0, w.1, w.0, rep(1, n.T), level,
#   #                          1, "tri", M, seed = NULL, useloop = TRUE)
#   #   return(c - q.sim)
#   # }
#   #
#   # c.min <- stats::qnorm(level)
#   # c.max <- stats::qnorm(1 - (1 - level)/(2 * n.T)) + 1
#   #
#   # res <- numeric(10)
#   # grid <- seq(from = c.min, to = c.max, length.out = 10)
#   # for(i in 1:10){
#   #   print(i)
#   #   res[i] <- eq(grid[i])
#   # }
#   #
#   # plot(grid, res, type = "l")
# })
#
#
# test_that("valid optimal quantile value (Lipschitz)", {
#
#   n <- 500
#   M <- 500
#   x <- stats::runif(n, min = -1, max = 1)
#   n.T <- 10
#   eval <- seq(from = -0.9, to = 0.9, length.out = n.T)
#   T.grad.mat <- rep(1, n.T)
#   level <- 0.95
#   y <- x + rnorm(n, 0, 1/4)
#
#   res <- NA
#
#   system.time({
#     opt.res <- opt_w("reg.Lip", 1, y, x, NULL, eval, T.grad.mat, level,
#                      1, "tri", M, seed = NULL, useloop = TRUE,
#                      root.robust = TRUE)
#     res <- opt.res$c.root
#     res.inc <- opt.res$increasing
#   })
#
#   expect_equal(!is.na(res), TRUE)
#   expect_equal(res.inc, TRUE)
# })
#
#
# test_that("valid optimal quantile value (Lipschitz-TE)", {
#
#   n <- 250
#   x.1 <- x.0 <- seq(-1, 1, length.out = n)
#   sd.true <- 1/2 + x.1^2
#   eps <- stats::rnorm(n, 0, sd.true)
#   y.1 <- x.1 + eps
#   y.0 <- x.1^2 + eps/2
#
#   y <- c(y.1, y.0)
#   x <- c(x.1, x.0)
#   d <- rep(c(1, 0), each = n)
#
#   M <- 500
#   n.T <- 10
#   eval <- seq(from = -0.9, to = 0.9, length.out = n.T)
#
#   T.grad.mat <- rep(1, n.T)
#   level <- 0.95
#
#   res <- NA
#
#   system.time({
#     opt.res <- opt_w("TE.Lip", 1, y, x, d, eval, T.grad.mat, level,
#                      1, "tri", M, seed = NULL, useloop = TRUE,
#                      root.robust = TRUE)
#     res <- opt.res$c.root
#     res.inc <- opt.res$increasing
#   })
#
#   expect_equal(!is.na(res), TRUE)
#   expect_equal(res.inc, TRUE)
# })
#
