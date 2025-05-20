library(testthat)

# Basic numeric matrices
set.seed(1)
X <- matrix(rnorm(20), 5, 4)
Y <- matrix(rnorm(15), 5, 3)

# Expect an error when both force_row_likeness and force_col_likeness are set
# simultaneously

test_that("genplsc stops if both force_row_likeness and force_col_likeness are supplied", {
  expect_error(
    genplsc(X, Y, dual = TRUE, force_row_likeness = TRUE, force_col_likeness = TRUE),
    "force_row_likeness"
  )
})
