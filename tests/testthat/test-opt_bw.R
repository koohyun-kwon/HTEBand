test_that("positive worst case bias", {
  x <- seq(-1, 1, length.out = 100)
  expect_equal(bias_Lip(x, 0, 1, "triangle", 0.1) > 0, TRUE)
})


test_that("positive standard deviation", {
  x <- seq(-1, 1, length.out = 100)
  eps <- rnorm(100, 0, 1/2 + x^2)
  y <- x + eps
  expect_equal(sd_Lip(y, x, 0, "triangle", 0.1, 0, FALSE) > 0, TRUE)
})
