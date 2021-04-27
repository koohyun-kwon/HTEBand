test_that("multiplication works", {
  res <- avar(rep(1/50, 50), rep(1/50, 50), rep(1, 50), rep(1, 50), 1)
  expect_equal(res > 0, TRUE)
})
