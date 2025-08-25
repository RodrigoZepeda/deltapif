# Helper function to create valid pif_atomic_class instance
create_valid_pif_atomic <- function() {
  pif(
    p = c(0.3, 0.7),
    p_cft = c(0.2, 0.8),
    beta = c(0.5, 1.0),
    var_p = matrix(c(0.01, 0.005, 0.005, 0.01), nrow = 2),
    var_beta = matrix(c(0.02, 0.01, 0.01, 0.02), nrow = 2),
    rr_link = function(x) exp(x),
    rr_link_deriv = function(x) exp(x),
    link = identity,
    link_deriv = function(x) 1,
    link_inv = identity,
    conf_level = 0.95,
    type = "PIF"
  )
}

test_that("Confidence levels work",{

  pif_atomic <- create_valid_pif_atomic()

  expect_equal(
    pif_atomic@conf_level,
    0.95
  )

  expect_equal(
    set_conf_level(pif_atomic, 0.9)@conf_level,
    0.9
  )

  expect_error(set_conf_level(pif_atomic, 1.2), "Invalid confidence")
  expect_error(set_conf_level(pif_atomic, -1.2), "Invalid confidence")
  expect_silent(set_conf_level(pif_atomic, 0.99))

})
