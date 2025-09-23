test_that("var_prod calculates correct values", {
  # Basic calculation
  expect_equal(var_prod(x = 2, y = 3, var_x = 0.5, var_y = 0.8),
               3^2 * 0.5 + 2^2 * 0.8 + 0.5 * 0.8)

  # When variances are zero
  expect_equal(var_prod(x = 5, y = 4, var_x = 0, var_y = 0), 0)

  # When one variance is zero
  expect_equal(var_prod(x = 5, y = 4, var_x = 2, var_y = 0), 4^2 * 2)
  expect_equal(var_prod(x = 5, y = 4, var_x = 0, var_y = 3), 5^2 * 3)

  # Larger numbers
  expect_equal(var_prod(x = 100, y = 200, var_x = 25, var_y = 36),
               200^2 * 25 + 100^2 * 36 + 25 * 36)
})

test_that("var_prod validates input types", {
  # Non-numeric inputs
  expect_error(var_prod(x = "2", y = 3, var_x = 0.5, var_y = 0.8))
  expect_error(var_prod(x = 2, y = TRUE, var_x = 0.5, var_y = 0.8))
  expect_error(var_prod(x = 2, y = 3, var_x = "0.5", var_y = 0.8))
  expect_error(var_prod(x = 2, y = 3, var_x = 0.5, var_y = NA))
})

test_that("var_prod validates input lengths", {
  # Vector inputs (should fail)
  expect_error(var_prod(x = c(1, 2), y = 3, var_x = 0.5, var_y = 0.8))
  expect_error(var_prod(x = 2, y = 3, var_x = c(0.5, 0.6), var_y = 0.8))
})

test_that("var_prod rejects negative variances", {
  expect_error(var_prod(x = 2, y = 3, var_x = -0.5, var_y = 0.8))
  expect_error(var_prod(x = 2, y = 3, var_x = 0.5, var_y = -0.8))
  expect_error(var_prod(x = 2, y = 3, var_x = -1, var_y = -2))
})

test_that("var_prod handles edge cases", {
  # Zero inputs
  expect_equal(var_prod(x = 0, y = 0, var_x = 0, var_y = 0), 0)
  expect_equal(var_prod(x = 0, y = 5, var_x = 2, var_y = 3), 5^2 * 2 + 0 + 2*3)
  expect_equal(var_prod(x = 5, y = 0, var_x = 2, var_y = 3), 0 + 5^2 * 3 + 2*3)

  # Very small variances
  expect_equal(var_prod(x = 1, y = 1, var_x = 1e-10, var_y = 1e-10),
               1e-10 + 1e-10 + 1e-20)
})

test_that("var_prod handles special values correctly", {
  # Infinite values
  expect_error(var_prod(x = Inf, y = 3, var_x = 0.5, var_y = 0.8))
  expect_error(var_prod(x = 2, y = 3, var_x = Inf, var_y = 0.8))

  # NA/NaN values
  expect_error(var_prod(x = NA, y = 3, var_x = 0.5, var_y = 0.8))
  expect_error(var_prod(x = 2, y = NaN, var_x = 0.5, var_y = 0.8))
})

test_that("var_prod handles negative x and y values", {
  # Negative means with positive variances
  expect_equal(var_prod(x = -2, y = 3, var_x = 0.5, var_y = 0.8),
               3^2 * 0.5 + (-2)^2 * 0.8 + 0.5 * 0.8)

  expect_equal(var_prod(x = 2, y = -3, var_x = 0.5, var_y = 0.8),
               (-3)^2 * 0.5 + 2^2 * 0.8 + 0.5 * 0.8)

  expect_equal(var_prod(x = -2, y = -3, var_x = 0.5, var_y = 0.8),
               (-3)^2 * 0.5 + (-2)^2 * 0.8 + 0.5 * 0.8)
})

test_that("var_prod returns non-negative results", {
  # Even with potential floating point issues, result should be non-negative
  result <- var_prod(x = 1e-10, y = 1e-10, var_x = 1e-20, var_y = 1e-20)
  expect_true(result >= 0)
})

test_that("var_prod matches theoretical expectations", {
  # For independent variables, var(X*Y) should be positive when variances are positive
  expect_true(var_prod(x = 10, y = 20, var_x = 1, var_y = 1) > 0)

  # When both variances are zero, result should be zero
  expect_equal(var_prod(x = 10, y = 20, var_x = 0, var_y = 0), 0)

  # The formula should be symmetric in a way
  result1 <- var_prod(x = 2, y = 3, var_x = 0.5, var_y = 0.8)
  result2 <- var_prod(x = 3, y = 2, var_x = 0.8, var_y = 0.5)
  expect_equal(result1, result2)
})

# Test for the warning on negative results (though theoretically shouldn't happen)
test_that("var_prod handles potential negative results gracefully", {
  # This test is for the edge case where numerical precision might cause issues
  # In practice, with the current formula, this shouldn't occur
  # But the function is prepared just in case
  expect_silent(var_prod(x = 1, y = 1, var_x = 0.1, var_y = 0.1))
})
