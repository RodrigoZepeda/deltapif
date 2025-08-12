test_that("check_links validates input combinations correctly", {
  # Test valid character link
  expect_silent(check_links("logit", NULL, NULL))
  expect_silent(check_links("log-complement", NULL, NULL))
  expect_silent(check_links("hawkins", NULL, NULL))

  # Test valid function links with proper derivatives/inverses
  expect_silent(check_links(logit, deriv_logit, inv_logit))
  expect_silent(check_links(log_complement, deriv_log_complement, inv_log_complement))
  expect_silent(check_links(hawkins, deriv_hawkins, inv_hawkins))

  # Test invalid combinations
  expect_error(
    check_links("logit", deriv_logit, NULL),
    "A logit link was specified but a `link_deriv` was given"
  )

  expect_error(
    check_links("logit", NULL, inv_logit),
    "A logit link was specified but a `link_inv` was given"
  )

  #Calculates its link deriv
  expect_silent(
    check_links(logit, link_deriv = NULL, inv_logit)
  )

  expect_error(
    check_links(logit, deriv_logit, NULL),
    "A functional link was specified but no `link_inv` was given"
  )

  # Test with custom functions
  custom_link <- function(x) x^2
  custom_inv <- function(x) sqrt(x)
  custom_deriv <- function(x) 2*x

  expect_silent(check_links(custom_link, custom_deriv, custom_inv))
  expect_error(
    check_links(custom_link, "log", custom_inv),
    "A functional link was specified but no `link_deriv` was given"
  )
})

test_that("check_rr_links validates input combinations correctly", {
  # Test valid character rr_link
  expect_silent(check_rr_links("log", NULL))
  expect_silent(check_rr_links("identity", NULL))

  # Test valid function rr_links
  expect_silent(check_rr_links(exp, function(x) exp(x)))
  expect_silent(check_rr_links(function(x) x^2, function(x) 2*x))

  # Test invalid combinations
  expect_error(
    check_rr_links("log", function(x) 1/x),
    "A log rr_link was specified but a `rr_link_deriv` was given"
  )

  # Test with custom functions
  custom_rr_link <- function(x) exp(x)
  custom_rr_deriv <- function(x) exp(x)

  expect_silent(check_rr_links(custom_rr_link, custom_rr_deriv))

  #No rr_link_deriv = calculates its link deriv
  expect_silent(check_rr_links(custom_rr_link, NULL))
})

test_that("check_links handles edge cases", {
  # NULL link
  expect_error(
    check_links(NULL, NULL, NULL),
    "No link was specified"
  )

  # NA inputs
  expect_error(
    check_links(NA, NULL, NULL),
    "No link was specified"
  )

  expect_error(
    check_links(matrix(1), NULL, NULL),
    "Cannot handle link of type"
  )

})

test_that("check_rr_links handles edge cases", {
  # NULL rr_link
  expect_error(
    check_rr_links(NULL, NULL),
    "No rr_link was specified"
  )

  # NA inputs
  expect_error(
    check_rr_links(NA, NULL),
    "No rr_link was specified"
  )

  # Other inputs
  expect_error(
    check_rr_links(08945, NULL),
    "Cannot handle rr_link of type"
  )

})

test_that("check_links validates function properties", {
  # Test with non-function inputs where functions are expected
  expect_error(
    check_links("logit", "not_a_function", NULL),
    "A logit link was specified but a `link_deriv` was given"
  )

  expect_error(
    check_links(logit, 22, inv_logit),
    "Cannot handle link_deriv of type "
  )

  expect_error(
    check_links(logit, deriv_logit, link_inv = 2),
    "Cannot handle link_inv of type "
  )
})

test_that("check_rr_links validates function properties", {
  # Test with non-function inputs where functions are expected
  expect_error(
    check_rr_links("log", "not_a_function"),
    "A log rr_link was specified but a `rr_link_deriv` was given"
  )

  expect_error(
    check_rr_links(exp, "not_a_function"),
    "rr_link_deriv must be a function when rr_link is a function"
  )

  expect_error(
    check_rr_links(2, NULL),
    "Cannot handle rr_link of type double"
  )

  expect_error(
    check_rr_links(exp, 2),
    "Cannot handle rr_link_deriv of type double"
  )

  #This one depends oncheckers but should pass:
  expect_silent(
    paf(p = 0.069, beta = 0.4637340, var_beta = 0.0273858,
        rr_link = exp, rr_link_deriv = exp,
        link = logit, link_inv = inv_logit, link_deriv = deriv_logit, var_p = 0)
  )
})
