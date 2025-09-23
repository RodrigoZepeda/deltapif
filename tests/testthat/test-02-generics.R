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


test_that("print.pif_global_ensemble_class prints components block", {
  p1 <- paf(0.1, 1.2, var_p = 0.01, var_beta = 0, label = "A",
            rr_link = identity, rr_link_deriv = function(x) 1)
  p2 <- paf(0.2, 1.3, var_p = 0.01, var_beta = 0, label = "B",
            rr_link = identity, rr_link_deriv = function(x) 1)
  ens <- pif_ensemble(p1, p2, weights = c(0.5, 0.5), label = "Ensemble")
  out <- capture.output(print(ens, accuracy = 0.001), type = "message")
  expect_true(any(grepl("Components:", out)))
  expect_true(any(grepl("A", out)))
  expect_true(any(grepl("B", out)))
})

test_that("print.cases_class works correctly", {
  p <- paf(0.2, 1.3, var_p = 0.001, var_beta = 0, label = "Label",
           rr_link = identity, rr_link_deriv = function(x) 1)
  cs <- averted_cases(1000, p, link = "identity", variance = 25, conf_level = 0.9)
  out <- capture.output(print(cs, accuracy = 1), type = "message")
  expect_true(any(grepl("Averted cases|Attributable cases", out)))
  expect_true(any(grepl("90% CI", out)))
})

test_that("print.covariance_structure_class summarizes shapes and values", {
  mat <- matrix(c(1, 0, 0, 2), ncol = 2, dimnames = list(c("r1","r2"), c("c1","c2")))
  covs <- as_covariance_structure(mat)
  out <- capture.output(print(covs))
  expect_true(any(grepl("r1", out)))
  expect_true(any(grepl("c2", out)))
  expect_true(any(grepl("1|2", out)))
})

test_that("coef/confint/summary/as.data.frame work for cases_class", {
  p <- paf(0.2, 1.3, var_p = 0.001, var_beta = 0, label = "P",
           rr_link = identity, rr_link_deriv = function(x) 1)
  cs <- averted_cases(1000, p, link = "identity", variance = 16, conf_level = 0.95)

  expect_equal(coef(cs), cs@cases)

  ci_default <- confint(cs)
  ci_custom  <- confint(cs, level = 0.90)
  expect_length(ci_default, 2)
  expect_length(ci_custom, 2)
  expect_true(ci_custom[1] <= ci_custom[2])

  s <- summary(cs)
  expect_named(s, c("cases", "standard_deviation", "ci_low", "ci_up", "confidence"))
  expect_equal(as.numeric(s["cases"]), cs@cases)
  expect_equal(as.numeric(s["standard_deviation"]), sqrt(cs@variance))

  df <- as.data.frame(cs)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 1)
  expect_equal(df$value, cs@cases)

  # multiple cases_class to as.data.frame
  cs2 <- averted_cases(500, p, link = "identity", variance = 9)
  df2 <- as.data.frame(cs, cs2)
  expect_equal(nrow(df2), 2)
})

test_that("names() for pif_global_ensemble_class returns children labels", {
  p1 <- paf(0.1, 1.2, var_p = 0.01, var_beta = 0, label = "X",
            rr_link = identity, rr_link_deriv = function(x) 1)
  p2 <- paf(0.2, 1.3, var_p = 0.01, var_beta = 0, label = "Y",
            rr_link = identity, rr_link_deriv = function(x) 1)
  tot <- pif_total(p1, p2, weights = c(0.4, 0.6), label = "T")
  expect_equal(names(tot), c("X","Y"))
})

test_that("children() returns NULL for atomic and labels for ensembles", {
  p <- paf(0.2, 1.1, var_p = 0.001, var_beta = 0, label = "Atom",
           rr_link = identity, rr_link_deriv = function(x) 1)
  p2 <- paf(0.2, 1.1, var_p = 0.001, var_beta = 0, label = "Atom2",
           rr_link = identity, rr_link_deriv = function(x) 1)
  expect_null(children(p))
  ens <- pif_ensemble(p, p2, weights = c(0.5, 0.5), label = "Ens")
  expect_equal(children(ens), c("Atom", "Atom2"))
})

test_that("row_names/col_names, length and as.matrix for covariance_structure", {
  m <- matrix(c(1,2,3,4), ncol = 2, dimnames = list(c("r1","r2"), c("c1","c2")))
  cs <- as_covariance_structure(m)
  expect_equal(row_names(cs), c("r1","r2"))
  expect_equal(col_names(cs), c("r1","r2")) # col_names returns names(x@cov_list)
  expect_equal(length(cs), 2L)

  mat <- as.matrix(cs)
  expect_true(is.matrix(mat))
  expect_equal(dim(mat), c(2,2))
  expect_equal(mat[1,1], 1)
})

test_that("subset() on covariance_structure supports numeric, character, rows/cols and negate", {
  m <- matrix(c(1,2,3,4), ncol = 2, dimnames = list(c("r1","r2"), c("c1","c2")))
  cs <- as_covariance_structure(m)

  # select by numeric
  cs12 <- subset(cs, 1:2)
  expect_equal(length(cs12), 2L)

  # select by character
  csc1 <- subset(cs, "r1")
  expect_equal(row_names(csc1), "r1")

  # rows/cols
  byrow <- subset(cs, rows = "r2")
  expect_equal(row_names(byrow), "r2")
  bycol <- subset(cs, cols = "c2")
  # after col subset, column names come from inner list names
  expect_true(all(names(bycol@cov_list[[1]]) == "c2"))

  # negate
  cs_neg <- subset(cs, select = "r1", negate = TRUE)
  expect_equal(row_names(cs_neg), "r2")

  # errors
  expect_error(subset(cs, select = 3), "outside the range")
  expect_error(subset(cs, select = 1.5), "non integer")
  expect_error(subset(cs, select = TRUE), "either numeric or character")
})

test_that("subset_row/subset_col internal helpers validate and subset", {
  m <- matrix(c(1,2,3,4), ncol = 2, dimnames = list(c("r1","r2"), c("c1","c2")))
  cs <- as_covariance_structure(m)

  expect_error(subset_row(list(), "r1"), "only supports `covariance_structures`")
  expect_error(subset_col(list(), "c1"), "only supports `covariance_structures`")
  expect_error(subset_row(cs, 1.2), "non integer")
  expect_error(subset_col(cs, 3), "outside the range")

  r <- subset_row(cs, 2)
  expect_equal(row_names(r), "r2")

  c <- subset_col(cs, "c1")
  expect_true(all(names(c@cov_list[[1]]) == "c1"))
})

test_that("as_covariance_structure handles scalars, vectors, data.frames and errors", {
  # scalar
  cs1 <- as_covariance_structure(2, col_names = "col", row_names = "row")
  expect_true(S7::S7_inherits(cs1, covariance_structure_class))
  expect_equal(row_names(cs1), "row")

  # vector
  v <- c(0.1, 0.2, 0.3)
  csv <- as_covariance_structure(v, row_names = c("r1","r2","r3"), col_names = "c1")
  expect_equal(row_names(csv), c("r1","r2","r3"))

  # data.frame
  df <- data.frame(a = c(1,2), b = c(3,4))
  csd <- as_covariance_structure(df, row_names = c("r1","r2"))
  expect_true(S7::S7_inherits(csd, covariance_structure_class))

  # errors on wrong col/row name lengths
  m <- matrix(1:4, ncol = 2)
  expect_error(as_covariance_structure(m, col_names = c("a","b","c")))
  expect_error(as_covariance_structure(m, row_names = c("r1")))
  expect_error(as_covariance_structure(list()))
})
