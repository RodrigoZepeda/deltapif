#Create pif objects
not_a_pif <- S7::new_class(name = "not_a_pif")
a_pif <- pif_atomic_class(p = 0.3, p_cft = 0.2, beta = 1.2, var_p = 0,
                          var_beta = 0, rr_link = identity,
                          rr_link_deriv = deriv_identity, link = logit,
                          link_inv = inv_logit, link_deriv = deriv_logit,
                          type = "PIF", conf_level = 0.95,
                          upper_bound_p = T, upper_bound_beta = T)
another_pif <- pif_class(pif = 0.2, variance = 0, conf_level = 0.95,
                         link = identity, link_inv = identity,
                         link_deriv = deriv_identity)


some_pifs <- S7::new_class(name = "some_pifs",
                           properties = list(
                             pif_list = S7::class_list,
                             pif_weights = S7::class_numeric
                           ))



test_that("Validators show errors",{

  expect_silent(
    validate_global_ensemble(
      some_pifs(
        pif_list = list(a_pif, another_pif),
        pif_weights = c(1,1)
      )
    )
  )

  expect_error(
    validate_global_ensemble(
      some_pifs(
        pif_list = list(a_pif, not_a_pif),
        pif_weights = c(1,1)
      )
    ),
    "must be a 'pif_class'"
  )

  expect_error(
    validate_global_ensemble(
      some_pifs(
        pif_list = list(a_pif, another_pif),
        pif_weights = c(1, 2, 3)
      )
    ),
    "have length"
  )

})
