context("test-things_break")

test_that("Error throwing", {
  expect_error(regularizeMxModel(mxModel(), "A", penalty_function = "GHG"))
})
