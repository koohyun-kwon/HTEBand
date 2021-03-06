test_that("positive worst case bias", {
  x <- seq(-1, 1, length.out = 100)
  expect_equal(bias_Lip(x, 0, 1, "tri", 0.1) > 0, TRUE)
})

test_that("positive worst case bias: Hölder", {
  x <- seq(-1, 1, length.out = 100)
  t <- 0
  expect_equal(bias_Hol(x, t, 1, "tri", 0.1) > 0, TRUE)
})


test_that("positive variance", {
  n <- 250
  x <- seq(-1, 1, length.out = n)
  sd.true <- 1/2 + x^2
  eps <- stats::rnorm(n, 0, sd.true)
  y <- x + eps
  expect_equal(var_Lip(y, x, 0, "tri", 0.1, 0) > 0, TRUE)

  # Additional test for difference between estimated sd and true sd
  # diff <- 0
  # true.val <- var_Lip_resid(x, 0, "tri", 0.1, sd.true)
  # for(i in 1:100){
  #
  #   print(i)
  #   eps <- stats::rnorm(n, 0, sd.true)
  #   y <- x + eps
  #
  #   diff <- diff + (var_Lip(y, x, 0, "tri", 0.1, 0) - true.val) / true.val
  # }
  #
  # diff / 100
})


test_that("positive variance (Hölder)", {
  n <- 250
  x <- seq(-1, 1, length.out = n)
  sd.true <- 1/2 + x^2
  eps <- stats::rnorm(n, 0, sd.true)
  y <- x + eps
  expect_equal(unname(var_Hol(y, x, 0, "tri", 0.1)) > 0, TRUE)
})

test_that("valid optimal bandwidth", {

  n <- 250
  x <- seq(-1, 1, length.out = n)
  sd.true <- 1/2 + x^2
  eps <- stats::rnorm(n, 0, sd.true)
  y <- x + eps

  res <- bw_opt(y, x, 0, TE = FALSE, d = NULL, 1, "tri", 0.05, bw.eq = TRUE,
                1)

  res

  expect_equal(res$h.opt > 0, TRUE)
  expect_equal(res$hl.opt > 0, TRUE)
})


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

  res <- bw_opt(y, x, 0, TE = TRUE, d = d, 1, "tri", 0.05, bw.eq = FALSE,
                1)
  res.eq <- bw_opt(y, x, 0, TE = TRUE, d = d, 1, "tri", 0.05, bw.eq = TRUE,
                   1)
  # res.orc <- bw_Lip_supp(c(sd.true, sd.true/2), x, 0, TE = TRUE, d = d, 1, "tri", 0.05, bw.eq = TRUE)

  res
  res.eq

  expect_equal(res$h.opt > 0, rep(TRUE, 2))
  expect_equal(res$hl.opt > 0, TRUE)
  expect_equal(res.eq$h.opt > 0, rep(TRUE, 2))
  expect_equal(res.eq$hl.opt > 0, TRUE)
})


test_that("valid optimal bandwidth(TE - Hölder)", {

  n <- 250
  x.1 <- x.0 <- seq(-1, 1, length.out = n)
  sd.true <- 1/2 + x.1^2
  eps <- stats::rnorm(n, 0, sd.true)
  y.1 <- x.1 + eps
  y.0 <- x.1^2 + eps/2

  y <- c(y.1, y.0)
  x <- c(x.1, x.0)
  d <- rep(c(1, 0), each = n)

  res <- bw_opt(y, x, 0, TE = TRUE, d = d, 1, "tri", 0.05, bw.eq = FALSE,
                p = 2)
  res.eq <- bw_opt(y, x, 0, TE = TRUE, d = d, 1, "tri", 0.05, bw.eq = TRUE,
                   p = 2)
  res
  res.eq

  expect_equal(unname(res$h.opt) > 0, rep(TRUE, 2))
  expect_equal(unname(res$hl.opt) > 0, TRUE)
  expect_equal(unname(res.eq$h.opt) > 0, rep(TRUE, 2))
  expect_equal(unname(res.eq$hl.opt) > 0, TRUE)
})

test_that("Check same variance estimator", {

  x <- seq(-1, 1, length.out = 500)
  y <- x^2 + rnorm(500, 0, 1/2)
  d <- RDHonest::LPPData(as.data.frame(cbind(y, x)), point = 0)
  d <- RDHonest::NPRPrelimVar.fit(d, se.initial = "EHW") # Add conditional variance values

  kern.rdh = "triangular"
  h = 1

  w <- lp_w(kern.rdh, order = 1, h, d$X)
  res1 <- sum(w^2 * d$sigma2)

  r1 <- RDHonest::NPRreg.fit(d, h, kern.rdh, order = 1, se.method = "supplied.var",
                             TRUE) # NPRreg.fit() accepts scalar h
  res2 <- r1$se["supplied.var"]^2
  expect_equal(unname(res1) - unname(res2) < 10^(-10), TRUE)
})
