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
# test_that("wiggly band? + running time", {
#
#   x <- seq(-1, 1, length.out = 500)
#   true.fun <- function(x) x^2
#   y <- true.fun(x) + rnorm(500, 0, 1/4)
#
#   system.time({
#     res <- NpregBand(y, x, 0.9, 0.95, "L", n.eval = 50, q.int = 0.1)
#   })
#
#
#   system.time({
#     res.H <- NpregBand(y, x, 10, 0.95, "H", n.eval = 50, q.int = 0.1)
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
#   ggplot(data = res.all, aes(x = eval, group = fspace)) + geom_line(aes(y = h.t, color = fspace)) +
#     ylim(0, NA)
#
#   h.test <- 0.17
#
#   sd.res.L <- numeric(length(res.all$eval)/2)
#   sd.res.H <- numeric(length(res.all$eval)/2)
#
#   for(i in 1:((length(res.all$eval)) / 2)){
#
#     sd.res.L[i] <- sqrt(var_Lip(y, x, res.all$eval[i], "triangle", h.test, 1, FALSE))
#
#     d <- RDHonest::LPPData(as.data.frame(cbind(y, x)), point = res.all$eval[i])
#     d <- RDHonest::NPRPrelimVar.fit(d, se.initial = "EHW")
#     r <- RDHonest::NPRHonest.fit(d, 10, "triangular", c(p = h.test, m = h.test),
#                        alpha = 0.05, se.method = "supplied.var", sclass = "H",
#                        order = 1, T0 = 0, T0bias = TRUE)
#     sd.res.H[i] <- r$sd
#
#   }
#
#   res.all.2 <- cbind(res.all, sd = c(sd.res.L, sd.res.H))
#
#   ggplot(data = res.all.2, aes(x = eval, group = fspace)) + geom_line(aes(y = sd, color = fspace)) +
#     ylim(0, NA) + geom_hline(yintercept = 0.03134209)
#
#
#
# })
