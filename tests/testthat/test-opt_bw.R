test_that("positive worst case bias", {
  x <- seq(-1, 1, length.out = 100)
  expect_equal(bias_Lip(x, 0, 1, "triangle", 0.1) > 0, TRUE)
})
