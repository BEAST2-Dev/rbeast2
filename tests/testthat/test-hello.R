context("hello")

test_that("hello: use", {
  expect_output(hello(), "Hello, world!")
})
