test_that("positive worst case bias", {
  x <- seq(-1, 1, length.out = 100)
  expect_equal(bias_Lip(x, 0, 1, "triangle", 0.1) > 0, TRUE)
})


test_that("positive variance", {
  n <- 250
  x <- seq(-1, 1, length.out = n)
  sd.true <- 1/2 + x^2
  eps <- stats::rnorm(n, 0, sd.true)
  y <- x + eps
  expect_equal(var_Lip(y, x, 0, "triangle", 0.1, 0, FALSE) > 0, TRUE)

  # Additional test for difference between estimated sd and true sd
  # diff <- 0
  # true.val <- var_Lip_resid(x, 0, "triangle", 0.1, sd.true)
  # for(i in 1:100){
  #
  #   print(i)
  #   eps <- stats::rnorm(n, 0, sd.true)
  #   y <- x + eps
  #
  #   diff <- diff + (var_Lip(y, x, 0, "triangle", 0.1, 0, FALSE) - true.val) / true.val
  # }
  #
  # diff / 100
})

test_that("valid gradient values", {

  n <- 250
  x <- seq(-1, 1, length.out = n)
  sd.true <- 1/2 + x^2
  eps <- stats::rnorm(n, 0, sd.true)
  y <- x + eps

  # increasing bandwidth should increase bias
  expect_equal(bias_Lip_gr(x, 0, 1, "triangle", 0.1) >= 0, TRUE)
  # increasing bandwidth should decrease sd
  expect_equal(sd_Lip_gr(y, x, 0, "triangle", 0.5, 0, FALSE) <= 0, TRUE)

  # Optional tests

  # h.min <- min(abs(x))
  # h.grid <- seq(h.min, 0.5, length.out = 50)
  # sd.res <- numeric(length(h.grid))
  #
  # for(i in 1:length(h.grid)) sd.res[i] <-
  #   2 * sd_Lip_gr(y, x, 0, "triangle", h.grid[i], 0, FALSE) + bias_Lip_gr(x, 0, 1, "triangle", h.grid[i])
  # res <- data.frame(h.grid = h.grid, sd.res = sd.res)
  # library(tidyverse)
  # ggplot(res, aes(x = h.grid, y = sd.res)) + geom_line() + ylim(-1, NA)
})

# test_that("valid optimal bandwidth", {
#
#   n <- 250
#   x <- seq(-1, 1, length.out = n)
#   sd.true <- 1/2 + x^2
#   eps <- stats::rnorm(n, 0, sd.true)
#   y <- x + eps
#   M <- 1
#   t <- 0
#   kern <- "triangle"
#
#   res <- bw_Lip(y, x, t, TE = FALSE, d = NULL, M, kern, 0.05, bw.eq = TRUE,
#                 1, FALSE)
#
#   res
#
#   expect_equal(res$h.opt > 0, TRUE)
#   expect_equal(res$hl.opt > 0, TRUE)
#
#   obj.1 <- function(h){
#
#     bias <- M * bias_Lip(x, t, M, kern, h)
#     sd <- sqrt(var_Lip(y, x, t, kern, h, deg, loo))
#     c <- stats::qnorm(1 - alpha) / 2
#     return(bias + c * sd)
#   }
# })


test_that("valid optimal bandwidth(TE)", {

  n <- 250
  x.1 <- x.0 <- seq(-1, 1, length.out = n)
  sd.true <- 1/2 + x.1^2
  eps <- stats::rnorm(n, 0, sd.true)
  y.1 <- x.1 + eps
  y.0 <- x.1^2 + eps/2

  y <- c(y.1, y.0)
  x <- c(x.1, x.0)
  d <- rep(c(1, 0), each = n)

  res <- bw_Lip(y, x, 0, TE = TRUE, d = d, 1, "triangle", 0.05, bw.eq = FALSE,
                1, FALSE)
  res.eq <- bw_Lip(y, x, 0, TE = TRUE, d = d, 1, "triangle", 0.05, bw.eq = TRUE,
                   1, FALSE)

  res
  res.eq

  expect_equal(res$h.opt > 0, rep(TRUE, 2))
  expect_equal(res$hl.opt > 0, TRUE)
  expect_equal(res.eq$h.opt > 0, rep(TRUE, 2))
  expect_equal(res.eq$hl.opt > 0, TRUE)
})
