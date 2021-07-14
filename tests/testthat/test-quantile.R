test_that("positive asymptotic variance", {
  res <- avar(rep(1/50, 50), rep(1/50, 50), rep(1, 50), rep(1, 50), 1)
  expect_equal(res > 0, TRUE)

  w.1 <- w.0 <- matrix(rep(1/50, 100), ncol = 2)
  omega.1 <- omega.0 <- array(rep(c(1,0,0,1), each = 50), dim = c(50, 2, 2))
  T.grad <- c(1, 1)
  res <- avar(w.1, w.0, omega.1, omega.0, T.grad)
  expect_equal(res > 0, TRUE)

})

test_that("Asymptotic variance when w.0 = 0", {

  res <- avar(rep(1/50, 50), 0, rep(1, 50), 0, 1)
  expect_equal(res > 0, TRUE)

  w.1 <- matrix(rep(1/50, 100), ncol = 2)
  w.0 <- matrix(0, nrow = 1, ncol = 2)
  omega.1 <- array(rep(c(1,0,0,1), each = 50), dim = c(50, 2, 2))
  omega.0 <- array(0, dim = c(1, 2, 2))
  T.grad <- c(1, 1)
  res <- avar(w.1, w.0, omega.1, omega.0, T.grad)
  expect_equal(res > 0, TRUE)

})

test_that("positive estimated variance (k = 1)",{

  x.1 <- x.0 <- seq(from = -1, to = 1, length.out = 500)
  y.1 <- x.1^2 + stats::rnorm(500, 0, 0.1)
  y.0 <- x.0^2 + stats::rnorm(500, 0, 0.1)
  w.1 <- w.0 <- rep(1/500, 500)
  z.1 <- rnorm(500 * 500)
  z.0 <- rnorm(500 * 500)
  res <- stud_err_sim(y.1, y.0, x.1, x.0, w.1, w.0, 1, 1, "triangle", z.1, z.0)$dnmnt
  omega.1 <- omega.0 <- rep(0.1^2, 500)
  res.true <- avar(w.1, w.0, omega.1, omega.0, 1)
  c(res, sqrt(res.true))

  expect_equal(!is.na(res), TRUE)

})

test_that("positive estimated variance (k = 1); w.0 = 0",{

  x.1 <- seq(from = -1, to = 1, length.out = 500)
  x.0 <- 0
  y.1 <- x.1^2 + stats::rnorm(500, 0, 0.1)
  y.0 <- 0
  w.1 <- rep(1/500, 500)
  w.0 <- 0
  z.1 <- rnorm(500 * 500)
  z.0 <- rep(0, 500)
  res <- stud_err_sim(y.1, y.0, x.1, x.0, w.1, w.0, 1, 1, "triangle", z.1, z.0)$dnmnt
  omega.1 <- rep(0.1^2, 500)
  omega.0 <- 0
  res.true <- avar(w.1, w.0, omega.1, omega.0, 1)
  c(res, sqrt(res.true))

  expect_equal(!is.na(res), TRUE)

})

test_that("positive estimated variance (k = 2)",{

  n <- 500
  sd <- 0.1
  x.1 <- x.0 <- seq(from = -1, to = 1, length.out = n)
  eps.1 <- stats::rnorm(n, 0, sd)
  y.11 <- x.1^2 + eps.1 + stats::rnorm(n, 0, sd)
  y.12 <- x.1 + eps.1 + stats::rnorm(n, 0, sd)
  eps.0 <- stats::rnorm(n, 0, sd)
  y.01 <- x.0^2 + eps.0 + stats::rnorm(n, 0, sd)
  y.02 <- x.0 + eps.0 + stats::rnorm(n, 0, sd)
  y.1 <- cbind(y.11, y.12)
  y.0 <- cbind(y.01, y.02)
  w.1 <- w.0 <- cbind(rep(1/n, n), rep(1/n, n))
  z.1 <- rnorm(n * 2 * 500)
  z.0 <- rnorm(n * 2 * 500)


  res <- stud_err_sim(y.1, y.0, x.1, x.0, w.1, w.0, c(1, 1), 1, "triangle", z.1, z.0)$dnmnt

  omega.1 <- omega.0 <- array(rep(c(2 * sd^2, sd^2, sd^2 , 2 * sd^2), each = n), dim = c(n, 2, 2))
  res.true <- avar(w.1, w.0, omega.1, omega.0, c(1, 1))
  c(res, sqrt(res.true))

  expect_equal(!is.na(res), TRUE)
})


test_that("positive estimated variance (k = 2); w.0 = 0",{

  n <- 500
  sd <- 0.1
  x.1 <- seq(from = -1, to = 1, length.out = n)
  x.0 <- 0
  y.11 <- x.1^2 + stats::rnorm(n, 0, sd)
  y.12 <- x.1 + stats::rnorm(n, 0, sd)
  y.01 <- 0
  y.02 <- 0
  y.1 <- cbind(y.11, y.12)
  y.0 <- cbind(y.01, y.02)
  w.1 <- cbind(rep(1/n, n), rep(1/n, n))
  w.0 <- cbind(0, 0)
  z.1 <- rnorm(n * 2 * 500)
  z.0 <- rep(0, 2 * 500)


  res <- stud_err_sim(y.1, y.0, x.1, x.0, w.1, w.0, c(1, 1), 1, "triangle", z.1, z.0)$dnmnt

  omega.1 <- array(rep(c(sd^2,0,0,sd^2), each = n), dim = c(n, 2, 2))
  omega.0 <- array(c(0, 0, 0, 0), dim = c(1, 2, 2))
  res.true <- avar(w.1, w.0, omega.1, omega.0, c(1, 1))
  c(res, sqrt(res.true))

  expect_equal(!is.na(res), TRUE)
})

test_that("Valid true and simulated quantile value", {

  n <- 500
  M <- 100
  x <- stats::runif(n, min = -1, max = 1)
  n.T <- 10
  eval <- seq(from = -0.9, to = 0.9, length.out = n.T)
  omega <- rep(1, n)
  T.grad.mat <- rep(1, n.T)
  level <- 0.95

  eps.1.mat <- rnorm(n * M)
  eps.0.mat <- rnorm(M)
  y <- x + rnorm(n, 0, 1)
  w <- array(w_get_Hol(y, x, eval, 1, 0.95)$w.mat, dim = c(n, 1, n.T))

  res <- sup_quant_orc(eps.1.mat, eps.0.mat, w, array(0, dim = c(1, 1, n.T)),
                       omega, 0, T.grad.mat, level, M, useloop = TRUE)
  res.est <- sup_quant_sim(y, 0, x, 0, w, array(rep(0, n.T), dim = c(1, 1, n.T)),
                           rep(1, n.T), level, 1, "triangle", 1000)
  c(res, res.est)

  expect_equal(as.numeric(res) > stats::qnorm(level), TRUE)
  expect_equal(as.numeric(res.est) > stats::qnorm(level), TRUE)

  # Optional test to investigate whether two values are close

  # test.val <- 0
  # for(i in 1:50){
  #
  #   print(i)
  #   eps.1.mat <- rnorm(n * M)
  #   eps.0.mat <- rnorm(M)
  #   y <- x + rnorm(n, 0, 1)
  #   w <- array(w_get_Hol(y, x, eval, 1, 0.95)$w.mat, dim = c(n, 1, n.T))
  #
  #   res <- sup_quant_orc(eps.1.mat, eps.0.mat, w, array(0, dim = c(1, 1, n.T)),
  #                        omega, 0, T.grad.mat, level, M, useloop = TRUE)
  #   res.est <- sup_quant_sim(y, 0, x, 0, w, array(rep(0, n.T), dim = c(1, 1, n.T)),
  #                            rep(1, n.T), level, 1, "triangle", 1000)
  #
  #   test.val = test.val + (res - res.est)/res
  # }
  #
  # test.val / 50

})


test_that("Valid true and simulated quantile value: heteroskedastic", {

  n <- 500
  M <- 100
  x <- stats::runif(n, min = -1, max = 1)
  n.T <- 10
  eval <- seq(from = -0.9, to = 0.9, length.out = n.T)
  omega <- 1/2 + x^2
  T.grad.mat <- rep(1, n.T)
  level <- 0.95

  eps.1.mat <- rnorm(n * M)
  eps.0.mat <- rnorm(M)
  y <- x + rnorm(n, 0, omega)
  w <- array(w_get_Hol(y, x, eval, 1, 0.95)$w.mat, dim = c(n, 1, n.T))

  res <- sup_quant_orc(eps.1.mat, eps.0.mat, w, array(0, dim = c(1, 1, n.T)),
                       omega, 0, T.grad.mat, level, M, useloop = TRUE)
  res.est <- sup_quant_sim(y, 0, x, 0, w, array(rep(0, n.T), dim = c(1, 1, n.T)),
                           rep(1, n.T), level, 1, "triangle", 1000)
  c(res, res.est)

  expect_equal(as.numeric(res) > stats::qnorm(level), TRUE)
  expect_equal(as.numeric(res.est) > stats::qnorm(level), TRUE)

  # Optional test to investigate whether two values are close

  # test.val <- 0
  # for(i in 1:50){
  #
  #   print(i)
  #   eps.1.mat <- rnorm(n * M)
  #   eps.0.mat <- rnorm(M)
  #   y <- x + rnorm(n, 0, 1)
  #   w <- array(w_get_Hol(y, x, eval, 1, 0.95)$w.mat, dim = c(n, 1, n.T))
  #
  #   res <- sup_quant_orc(eps.1.mat, eps.0.mat, w, array(0, dim = c(1, 1, n.T)),
  #                        omega, 0, T.grad.mat, level, M, useloop = TRUE)
  #   res.est <- sup_quant_sim(y, 0, x, 0, w, array(rep(0, n.T), dim = c(1, 1, n.T)),
  #                            rep(1, n.T), level, 1, "triangle", 1000)
  #
  #   test.val = test.val + (res - res.est)/res
  # }
  #
  # test.val / 50

})

