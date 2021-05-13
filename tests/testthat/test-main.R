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
#
#
# test_that(" running time", {
#
#   x <- seq(-1, 1, length.out = 500)
#   true.fun <- function(x) x^2
#   y <- true.fun(x) + rnorm(500, 0, abs(x)/4 + 1/16)
#
#   system.time({
#     res <- NpregBand(y, x, 2, 0.95, "L", n.eval = 50)
#   })
#
#   method = "reg.Lip"
#   eval = res$eval
#   C = 2
#   kern.reg = "triangle"
#   deg = 1
#   loo = FALSE
#   se.method = "resid"
#   resid.1 = eps_hat(y, x, deg, kern.reg, loo)
#   resid.0 <- 0
#   d = rep(1, length(y))
#   k = 1
#   n.T = length(eval)
#   T.grad.mat <- rep(1, n.T)
#
#   opt_w_2(2.5, method, C, y, x, d, eval, T.grad.mat, 0.95,
#                       1, "triangle", loo = FALSE, M = 500, seed = NULL, useloop = TRUE,
#                       root.robust = FALSE, ng = 10)
#
#
#   system.time({
#     res.H <- NpregBand(y, x, 2, 0.95, "H", n.eval = 50)
#   })
#
#   res.all <- cbind(rbind(res, res.H), fspace = rep(c("L", "H"), each = nrow(res)),
#                    truef = true.fun(res$eval))
#
#
#   library(tidyverse)
#
#   ggplot(data = res.all, aes(x = eval, group = fspace)) + geom_line(aes(y = cb.lower, color = fspace)) +
#     geom_line(aes(y = cb.upper, color = fspace)) + geom_line(aes(y = truef))
#
#   ggplot(data = res.all, aes(x = eval, group = fspace)) + geom_line(aes(y = h.t, color = fspace))
# })
