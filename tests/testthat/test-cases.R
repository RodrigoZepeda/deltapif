create_test_pif <- function(type = "PIF") {
  pif(0.3, 0.1, 0.12, 0.05, 0.02, label = "pif")
}

create_test_paf <- function() {
  create_test_pif(type = "PAF")
}

# Test averted_cases function
test_that("averted_cases returns a correct class", {
  pif_obj <- create_test_pif()
  result <- averted_cases(cases = 100, pif = pif_obj)

  expect_true(S7::S7_inherits(result, cases_class))
})

test_that("averted_cases calculates correct values with identity link", {
  pif_obj <- create_test_pif()
  result <- averted_cases(cases = 100, pif = pif_obj, link = "identity")

  expect_equal(result@cases, pif_obj@pif*100)
})

test_that("averted_cases calculates correct values with log link", {
  pif_obj <- create_test_pif()
  result <- averted_cases(cases = 100, pif = pif_obj, link = "log")

  # Should still calculate basic attributable cases correctly
  expect_equal(result@cases, pif_obj@pif*100)
})

test_that("averted_cases handles variance calculation", {
  pif_obj <- create_test_pif()
  result <- averted_cases(cases = 100, pif = pif_obj, variance = 25)

  # Test that variance is incorporated (exact formula would depend on implementation)
  expect_true(!is.null(result@variance))
  expect_true(result@variance > 0)
})

test_that("averted_cases validates link function", {
  pif_obj <- create_test_pif()

  expect_error(
    averted_cases(cases = 100, pif = pif_obj, link = "invalid_link"),
    "Invalid link invalid_link"
  )
})

test_that("averted_cases works with custom link functions", {
  pif_obj <- create_test_pif()

  # Test with custom inverse link and derivative
  custom_inv <- function(x) 1/x
  custom_deriv <- function(x) - 1/x^2

  result <- averted_cases(
    cases = 100,
    pif = pif_obj,
    link = log,
    link_inv = custom_inv,
    link_deriv = custom_deriv
  )

  expect_true(result@variance > 0)
  expect_equal(result@cases, pif_obj@pif*100)

})

test_that("averted_cases handles different confidence levels", {
  pif_obj <- create_test_pif()

  result_90 <- averted_cases(cases = 100, pif = pif_obj, conf_level = 0.90)
  result_99 <- averted_cases(cases = 100, pif = pif_obj, conf_level = 0.99)

  expect_equal(result_90@conf_level, 0.90)
  expect_equal(result_99@conf_level, 0.99)
})

# Test attributable_cases function
test_that("attributable_cases returns correct structure", {
  paf_obj <- paf(0.3, 0.1, 0.1, 0.02)
  result  <- attributable_cases(cases = 100, paf = paf_obj)
  result2 <- averted_cases(cases = 100, pif = paf_obj)

  expect_identical(result, result2)
})

test_that("attributable_cases validates PAF input", {
  pif_obj <- pif(0.3,0.2, 0.1, 0.1, 0.02) # This is PIF, not PAF

  expect_error(
    attributable_cases(cases = 100, paf = pif_obj),
    "Can only estimate.*with a population attributable fraction"
  )
})

test_that("attributable_cases passes parameters correctly to averted_cases", {
  paf_obj <- paf(0.3, 0.1, 0.1, 0.02)

  result <- attributable_cases(
    cases = 100,
    paf = paf_obj,
    variance = 25,
    conf_level = 0.90,
    link = "log"
  )

  expect_equal(result@cases, paf_obj@pif*100)
  expect_equal(result@conf_level, 0.90)
})

# Test edge cases
test_that("functions handle zero cases", {
  pif_obj <- create_test_pif()
  paf_obj <- paf(0.3, 0.1, 0.1, 0.02)

  result1 <- averted_cases(cases = 0, pif = pif_obj)
  result2 <- attributable_cases(cases = 0, paf = paf_obj)

  expect_equal(result1@cases, 0)
  expect_equal(result2@cases, 0)
})

test_that("functions handle negative cases", {
  pif_obj <- create_test_pif()
  paf_obj <- paf(0.3, 0.1, 0.1, 0.02)

  result1 <- averted_cases(cases = -50, pif = pif_obj)
  result2 <- attributable_cases(cases = -50, paf = paf_obj)

  expect_equal(result1@cases, -50*pif_obj@pif)
  expect_equal(result2@cases, -50*paf_obj@pif)
})

test_that("functions handle very large cases", {
  pif_obj <- create_test_pif()
  paf_obj <- paf(0.3, 0.1, 0.1, 0.02)

  large_cases <- 1e6
  result1 <- averted_cases(cases = large_cases, pif = pif_obj)
  result2 <- attributable_cases(cases = large_cases, paf = paf_obj)

  expect_equal(result1@cases, large_cases * pif_obj@pif)
  expect_equal(result2@cases, large_cases * paf_obj@pif)
})

test_that("functions handle boundary PIF values", {
  # Test with PIF = 0
  pif_zero <- paf(0.0, 0.1, 0.1, 0.02)

  result <- averted_cases(cases = 100, pif = pif_zero)
  expect_equal(result@cases, 0)

  # Test with PIF = 1
  pif_one <- paf(p = 1.0, beta = 100, 0.1, 0.02)
  result <- averted_cases(cases = 100, pif = pif_one)
  expect_equal(result@cases, 100)
})

# Test error conditions
test_that("functions validate input types", {
  pif_obj <- create_test_pif()
  paf_obj <- paf(p = 0.1, beta = 2, 0.1, 0.02)

  expect_error(averted_cases(cases = "100", pif = pif_obj))
  expect_error(attributable_cases(cases = "100", paf = paf_obj))

  expect_error(averted_cases(cases = 100, pif = "not_a_pif"))
  expect_error(attributable_cases(cases = 100, paf = "not_a_paf"))
})

test_that("functions validate confidence level bounds", {
  pif_obj <- create_test_pif()
  paf_obj <- create_test_paf()

  expect_error(averted_cases(cases = 100, pif = pif_obj, conf_level = 1.5))
  expect_error(averted_cases(cases = 100, pif = pif_obj, conf_level = -0.1))

  expect_error(attributable_cases(cases = 100, paf = paf_obj, conf_level = 1.5))
  expect_error(attributable_cases(cases = 100, paf = paf_obj, conf_level = -0.1))
})

