test_that("public PIF operations preserve the parent class", {
  first <- pif(
    p = c(0.3, 0.7),
    p_cft = c(0.2, 0.8),
    beta = c(0.5, 1),
    var_p = diag(0.01, 2),
    var_beta = diag(0.01, 2),
    link = "identity",
    label = "first"
  )
  second <- pif(
    p = c(0.3, 0.7),
    p_cft = c(0.2, 0.8),
    beta = c(0.5, 1),
    var_p = diag(0.01, 2),
    var_beta = diag(0.01, 2),
    link = "identity",
    label = "second"
  )
  total <- pif_total(
    first,
    second,
    weights = c(0.5, 0.5),
    var_weights = matrix(0, 2, 2),
    label = "total"
  )
  ensemble <- pif_ensemble(first, second)
  expect_equal(coef(total), 0.5 * coef(first) + 0.5 * coef(second))
  expect_equal(coef(ensemble), 1 - (1 - coef(first)) * (1 - coef(second)))
  expect_identical(total@label, "total")
})
