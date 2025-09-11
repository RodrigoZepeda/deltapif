test_that("print_pif_class works correctly", {
  # Create a mock pif_class object
  pif_obj <- pif_class(
    pif = 0.0421, # 4.21%
    variance = 0.0049, # sd = 0.07 or 7%
    conf_level = 0.90,
    type = "PIF",
    link = logit,
    link_inv = inv_logit,
    link_deriv = deriv_logit,
    label = "pif_obj"
  )

  # Capture the output
  output <- capture.output(print_pif_class(pif_obj, accuracy = 0.001), type = "message")

  # Verify the output contains expected elements
  expect_true(any(grepl("Potential Impact Fraction", output)))
  expect_true(any(grepl("PIF = 4.2", output)))
  expect_true(any(grepl("90% CI", output)))
  expect_true(any(grepl("0.25", output)))
  expect_true(any(grepl("43.30", output)))
  #expect_true(any(grepl("standard_deviation\\(pif.*= 7", output)))
  #expect_true(any(grepl("standard_deviation\\(link.* = 1.7", output)))

  # Test with PAF
  paf_obj <- pif_obj
  paf_obj@type <- "PAF"
  output <- capture.output(print_pif_class(paf_obj, accuracy = 0.001), type = "message")
  expect_true(any(grepl("Population Attributable Fraction", output)))
  #expect_true(any(grepl("standard_deviation\\(paf", output)))
})

test_that("coef.pif_class works correctly", {
  pif_obj <- pif_class(
    pif = 0.3,
    variance = 0.01,
    conf_level = 0.95,
    type = "PIF",
    link = identity,
    link_inv = identity,
    link_deriv = function(x) 1
  )

  expect_equal(coef(pif_obj), 0.3)
  expect_equal(coef(pif_obj, "ignored"), 0.3) # Additional params ignored
})

test_that("confint.pif_class works correctly", {
  pif_obj <- pif_class(
    pif = 0.3,
    variance = 0.01,
    conf_level = 0.95,
    type = "PIF",
    link = identity,
    link_inv = identity,
    link_deriv = function(x) 1
  )

  # Test default level
  expect_equal(length(confint(pif_obj)), 2)

  # Test custom level
  expect_equal(length(confint(pif_obj, level = 0.9)), 2)

  expect_equal(confint(pif_obj, level = pif_obj@conf_level), pif_obj@ci)
})

test_that("summary.pif_class works correctly", {
  pif_obj <- pif(
    p = 0.1,
    beta = 0.2,
    rr_link = "exp",
    var_p = 0.0042,
    var_beta = 0.012,
    conf_level = 0.9
  )

  result <- summary(pif_obj)

  expect_named(result, c("PIF", "standard_deviation", "ci_low", "ci_up", "confidence"))
  expect_equal(as.numeric(result["PIF"]), pif_obj@pif)
  expect_equal(as.numeric(result["standard_deviation"]), sqrt(pif_obj@variance)) # sqrt(0.01)
  expect_equal(as.numeric(result["ci_low"]), pif_obj@ci[1])
  expect_equal(as.numeric(result["ci_up"]), pif_obj@ci[2])
  expect_equal(as.numeric(result["confidence"]), pif_obj@conf_level)

  # Test with custom level
  result <- summary(pif_obj, level = 0.8)
  expect_equal(as.numeric(result["confidence"]), 0.8)

  # Test with PAF
  paf_obj <- pif_obj
  paf_obj@type <- "PAF"
  result <- summary(paf_obj)
  expect_named(result, c("PAF", "standard_deviation", "ci_low", "ci_up", "confidence"))
})

test_that("as.data.frame.pif_class works correctly", {
  pif_obj <- pif(
    p = 0.1,
    beta = 0.2,
    rr_link = "exp",
    var_p = 0.0042,
    var_beta = 0.012,
    conf_level = 0.9
  )

  df <- as.data.frame(pif_obj)

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 1)
  expect_equal(ncol(df), 7)
  expect_equal(df$value, pif_obj@pif)
  expect_equal(df$standard_deviation, sqrt(pif_obj@variance))

  # Test with custom level
  df <- as.data.frame(pif_obj, level = 0.7)
  expect_equal(df$confidence, 0.7)

  # Test with PAF
  paf_obj <- pif_obj
  paf_obj@type <- "PAF"
  df <- as.data.frame(paf_obj)
  expect_equal(names(df)[1], "value")
  expect_equal(df[1,"type"], "PAF")
})

test_that("names work", {

  #Play with assigning and changing name
  my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, var_p = 0.1,
                  var_beta = 0.2, label = "Old name")
  expect_equal(names(my_pif), "Old name")

  #A pif composed of others
  my_pif1 <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, var_p = 0.1,
                  var_beta = 0.2, label = "Test 1")
  my_pif2 <- pif(p = 0.4, p_cft = 0.1, beta = 1.3, var_p = 0.1,
                  var_beta = 0.2, label = "Test 2")
  pif_tot <- pif_total(my_pif1, my_pif2, weights = c(0.2, 0.8),
                  label = "Parent")

  #Check names
  expect_equal(as.vector(names(pif_tot)), c("Test 1", "Test 2"))
  expect_equal(names(pif_tot), children(pif_tot))


})
