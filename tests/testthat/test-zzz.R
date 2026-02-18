# Tests for zzz.R — the [[ and [[<- S7 methods registered in .onLoad
# (methods_register() itself is excluded)

# [[ extraction ----------------------------------------------------------------
test_that("[[ extracts inner list by integer index", {
  cov <- covariance_structure_class(
    list(
      "pif1" = list("pif1" = 0.1, "pif2" = 0.2),
      "pif2" = list("pif1" = 0.2, "pif2" = 0.5)
    )
  )
  result <- cov[[1]]
  expect_type(result, "list")
  expect_named(result, c("pif1", "pif2"))
})

test_that("[[ extracts inner list by character name", {
  cov <- covariance_structure_class(
    list(
      "A" = list("A" = 1.0, "B" = 0.3),
      "B" = list("A" = 0.3, "B" = 2.0)
    )
  )
  result <- cov[["B"]]
  expect_type(result, "list")
  expect_equal(result[["A"]], 0.3)
  expect_equal(result[["B"]], 2.0)
})

test_that("[[ chaining works to reach a scalar value", {
  cov <- covariance_structure_class(
    list(
      "pif1" = list("pif1" = 0.21, "pif2" = 0.12, "pif3" = 0.31),
      "pif2" = list("pif1" = 0.12, "pif2" = 0.33, "pif3" = -0.01),
      "pif3" = list("pif1" = 0.31, "pif2" = -0.01, "pif3" = 0.80)
    )
  )
  expect_equal(cov[[1]][[3]], 0.31)
  expect_equal(cov[["pif3"]][["pif3"]], 0.80)
  expect_equal(cov[["pif2"]][["pif3"]], -0.01)
})

test_that("[[ chaining works when entries are matrices", {
  mat <- matrix(c(1, 2, 3, 4), ncol = 2)
  cov <- covariance_structure_class(
    list(
      "pif1" = list("pif1" = 0.5, "pif2" = mat),
      "pif2" = list("pif1" = mat, "pif2" = 0.9)
    )
  )
  expect_equal(cov[["pif1"]][["pif2"]], mat)
  expect_equal(cov[[2]][[1]], mat)
})

# [[<- assignment --------------------------------------------------------------

test_that("[[<- replaces an entry with a list of correct length", {
  cov <- covariance_structure_class(
    list(
      "pif1" = list("pif1" = 0.1, "pif2" = 0.2),
      "pif2" = list("pif1" = 0.2, "pif2" = 0.5)
    )
  )
  new_val <- list("pif1" = 99, "pif2" = 88)
  cov[[1]] <- new_val
  expect_equal(cov[[1]], new_val)
})

test_that("[[<- replaces entry by character name", {
  cov <- covariance_structure_class(
    list(
      "X" = list("X" = 1.0, "Y" = 0.0),
      "Y" = list("X" = 0.0, "Y" = 1.0)
    )
  )
  cov[["X"]] <- list("X" = 5.0, "Y" = 3.0)
  expect_equal(cov[["X"]][["X"]], 5.0)
  expect_equal(cov[["X"]][["Y"]], 3.0)
})

test_that("[[<- errors when value is not a list", {
  cov <- covariance_structure_class(
    list(
      "pif1" = list("pif1" = 0.1, "pif2" = 0.2),
      "pif2" = list("pif1" = 0.2, "pif2" = 0.5)
    )
  )
  expect_error(
    { cov[[1]] <- 42 },
    "should be a list"
  )
  expect_error(
    { cov[["pif2"]] <- "string" },
    "should be a list"
  )
})

test_that("[[<- errors when replacement list has wrong length", {
  cov <- covariance_structure_class(
    list(
      "pif1" = list("pif1" = 0.1, "pif2" = 0.2),
      "pif2" = list("pif1" = 0.2, "pif2" = 0.5)
    )
  )
  # inner lists should have length 2, but we supply length 1
  expect_error(
    { cov[[1]] <- list("pif1" = 0.9) },
    "length"
  )
  # and length 3
  expect_error(
    { cov[["pif2"]] <- list("a" = 1, "b" = 2, "c" = 3) },
    "length"
  )
})

test_that("[[<- preserves unchanged rows after assignment", {
  cov <- covariance_structure_class(
    list(
      "pif1" = list("pif1" = 0.1, "pif2" = 0.2),
      "pif2" = list("pif1" = 0.2, "pif2" = 0.5)
    )
  )
  original_row2 <- cov[[2]]
  cov[[1]] <- list("pif1" = 7, "pif2" = 8)
  expect_equal(cov[[2]], original_row2)
})

test_that("repeated [[<- assignments accumulate correctly", {
  cov <- covariance_structure_class(
    list(
      "pif1" = list("pif1" = 0.1, "pif2" = 0.2),
      "pif2" = list("pif1" = 0.2, "pif2" = 0.5)
    )
  )
  cov[["pif1"]] <- list("pif1" = 10, "pif2" = 20)
  cov[["pif2"]] <- list("pif1" = 30, "pif2" = 40)

  expect_equal(cov[["pif1"]][["pif1"]], 10)
  expect_equal(cov[["pif1"]][["pif2"]], 20)
  expect_equal(cov[["pif2"]][["pif1"]], 30)
  expect_equal(cov[["pif2"]][["pif2"]], 40)
})

