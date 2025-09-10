#Create pif objects
not_a_pif <- S7::new_class(name = "not_a_pif")
a_pif <- pif(p = 0.3, p_cft = 0.2, beta = 1.2, var_p = 0,var_beta = 0.2)
another_pif <-  pif(p = 0.2, p_cft = 0.1, beta = 1.3, var_p = 0,var_beta = 0.1)

some_pifs <- S7::new_class(name = "some_pifs",
                           properties = list(
                             pif_list = S7::class_list,
                             weights = S7::class_numeric
                           ))



test_that("Validators show errors",{

  piflist <- list(a_pif, another_pif)
  names(piflist) <- c(a_pif@label, another_pif@label)

  not_a_piflist <- list(not_a_pif, another_pif)
  names(not_a_piflist) <- c("not", another_pif@label)

  expect_silent(
    validate_global_ensemble(
      some_pifs(
        pif_list = piflist,
        weights = c(1,1)
      )
    )
  )

  expect_error(
    validate_global_ensemble(
      some_pifs(
        pif_list = not_a_piflist,
        weights = c(1,1)
      )
    ),
    "must be a 'pif_class'"
  )

  expect_error(
    validate_global_ensemble(
      some_pifs(
        pif_list = list(a_pif, another_pif),
        weights = c(1,1)
      )
    ),
    "named list"
  )

  expect_error(
    validate_global_ensemble(
      some_pifs(
        pif_list = piflist,
        weights = c(1, 2, 3)
      )
    ),
    "have length"
  )

})
