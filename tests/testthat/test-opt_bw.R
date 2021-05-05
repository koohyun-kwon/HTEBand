test_that("positive worst case bias", {
  x <- seq(-1, 1, length.out = 100)
  expect_equal(bias_Lip(x, 0, 1, "triangle", 0.1) > 0, TRUE)
})


test_that("positive standard deviation", {
  n <- 250
  x <- seq(-1, 1, length.out = n)
  sd.true <- 1/2 + x^2
  eps <- stats::rnorm(n, 0, sd.true)
  y <- x + eps
  expect_equal(var_Lip(y, x, 0, "triangle", 0.1, 0, FALSE) > 0, TRUE)

  # Additional test for difference between estimated sd and true sd
  # diff <- 0
  # true.val <- var_Lip_true(x, 0, "triangle", 0.1, sd.true)
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
