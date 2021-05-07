# test_that("optional test", {
#   x <- rep(seq(-1, 1, length.out = 500), each = 2)
#   d <- rep(c(0, 1), 500)
#   y <- d * x^2 + (1 - d) * x + rnorm(500, 0, 1/4)
#   system.time({
#
#     res <- CATEBand(y, x, d, 2, 0.95, n.eval = 50)
#   })
#
#   plot(res$eval, res$cb.lower, type = "l")
#   lines(res$eval, res$cb.upper)
#   lines(res$eval, res$eval^2 - res$eval)
#
# })
