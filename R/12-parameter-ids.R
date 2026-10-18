#' Parameter identifiers
#'
#' Functions to handle the identifiers of the parameters `p` and `beta`.
#' A parameter is given an identifier by passing it as a named list with
#' a single element, for example `beta = list("smoking" = 0.46)` or
#' `p = list("survey_2020" = c(0.12, 0.30))`. Fractions whose parameters
#' share an identifier are assumed to use the same estimate (and hence
#' are correlated). Parameters without an identifier are independent of
#' the parameters of every other fraction.
#'
#' @param x Either a numeric vector or a named list with one element.
#' @param parameter Name of the parameter (`"p"` or `"beta"`) for the messages.
#'
#' @return A list with the `value` of the parameter and its `id`
#' (`NA_character_` if no identifier was given).
#'
#' @keywords internal
#' @noRd
parse_parameter_id <- function(x, parameter) {

  if (!is.list(x)) {
    return(list(value = x, id = NA_character_))
  }

  if (length(x) != 1 || is.null(names(x)) || is.na(names(x)) || names(x) == "") {
    cli::cli_abort(
      paste0(
        "`{parameter}` should be a numeric vector or a named list with a single ",
        "element whose name is the identifier of the estimate, for example ",
        "{.code {parameter} = list(\"my_id\" = 0.12)}."
      )
    )
  }

  if (!is.numeric(x[[1]])) {
    cli::cli_abort(
      "The element {.val {names(x)}} of `{parameter}` should be numeric."
    )
  }

  list(value = x[[1]], id = names(x))
}

#' Get the identifier of a parameter of an atomic fraction
#' @keywords internal
#' @noRd
parameter_id <- function(pif, parameter) {
  if (parameter == "beta") pif@beta_id else pif@p_id
}

#' Get the value of a parameter of an atomic fraction
#' @keywords internal
#' @noRd
parameter_value <- function(pif, parameter) {
  if (parameter == "beta") pif@beta else pif@p
}

#' Get the covariance of a parameter of an atomic fraction
#' @keywords internal
#' @noRd
parameter_variance <- function(pif, parameter) {
  if (parameter == "beta") pif@var_beta else pif@var_p
}

#' Atomic fractions of a fraction
#'
#' Returns the atomic fractions that make up a fraction (without duplicates).
#'
#' @keywords internal
#' @noRd
atomic_fractions <- function(pif) {
  flat <- as.list(pif)
  if (S7::S7_inherits(flat, pif_class)) {
    flat <- list(flat)
  }
  flat <- flat[sapply(flat, function(x) S7::S7_inherits(x, pif_atomic_class))]
  flat[!duplicated(sapply(flat, function(x) x@label))]
}

#' Collect the covariance matrices between parameter identifiers
#'
#' Collects the `var_p_ids` (respectively `var_beta_ids`) matrices given to
#' the totals and ensembles that make up a fraction.
#'
#' @keywords internal
#' @noRd
collect_id_covariances <- function(pif, parameter) {

  if (S7::S7_inherits(pif, pif_atomic_class)) {
    return(list())
  }

  own <- if (parameter == "beta") pif@var_beta_ids else pif@var_p_ids
  res <- if (is.null(own)) list() else list(own)

  for (child in pif@pif_list) {
    res <- c(res, collect_id_covariances(child, parameter))
  }

  res
}

#' Covariance block between two parameter identifiers
#'
#' Looks for the covariance between the parameters with identifiers
#' `id1` and `id2` in a list of matrices whose row and column names are the
#' identifiers (repeated once per element of vector-valued parameters).
#'
#' @return The covariance block as a matrix or `NULL` if not found.
#'
#' @keywords internal
#' @noRd
id_covariance_block <- function(id_matrices, id1, id2, n1, n2) {

  for (id_matrix in id_matrices) {
    rows <- which(rownames(id_matrix) == id1)
    cols <- which(colnames(id_matrix) == id2)

    if (length(rows) > 0 && length(cols) > 0) {
      if (length(rows) != n1 || length(cols) != n2) {
        cli::cli_abort(
          paste0(
            "The covariance between identifiers {.val {id1}} and {.val {id2}} ",
            "should have {n1} row{?s} and {n2} column{?s} but has ",
            "{length(rows)} and {length(cols)}. Repeat the identifier of a ",
            "vector-valued parameter once per element in the row and column ",
            "names of the matrix."
          )
        )
      }
      return(unname(id_matrix[rows, cols, drop = FALSE]))
    }
  }

  NULL
}

#' Check that fractions sharing an identifier use the same estimate
#'
#' @keywords internal
#' @noRd
check_shared_ids <- function(pif_list) {

  atomic <- unlist(lapply(pif_list, atomic_fractions), recursive = FALSE)
  atomic <- atomic[!duplicated(sapply(atomic, function(x) x@label))]

  for (parameter in c("p", "beta")) {
    ids <- sapply(atomic, parameter_id, parameter = parameter)
    for (id in unique(ids[!is.na(ids)])) {
      same_id <- atomic[which(ids == id)]
      first   <- same_id[[1]]
      for (other in same_id[-1]) {
        same_value <- isTRUE(all.equal(
          unname(parameter_value(first, parameter)),
          unname(parameter_value(other, parameter))
        ))
        same_var <- isTRUE(all.equal(
          unname(as.matrix(parameter_variance(first, parameter))),
          unname(as.matrix(parameter_variance(other, parameter)))
        ))
        if (!same_value || !same_var) {
          cli::cli_abort(
            paste0(
              "Fractions {.val {first@label}} and {.val {other@label}} use the ",
              "identifier {.val {id}} for `{parameter}` but have different ",
              "values or variances. The same identifier must refer to the same ",
              "estimate."
            )
          )
        }
      }
    }
  }

  invisible(TRUE)
}

#' Validate a covariance matrix between parameter identifiers
#'
#' Checks that `id_matrix` is a symmetric matrix whose row and column names
#' are identifiers of the parameters of the fractions in `pif_list` (repeated
#' once per element) and that its diagonal blocks equal the variances given
#' to each fraction.
#'
#' @keywords internal
#' @noRd
check_id_covariance <- function(id_matrix, pif_list, parameter) {

  if (is.null(id_matrix)) {
    return(NULL)
  }

  arg <- paste0("var_", parameter)

  if (!is.matrix(id_matrix) || !is.numeric(id_matrix)) {
    cli::cli_abort(
      "`{arg}` should be a numeric matrix whose row and column names are the identifiers of `{parameter}`."
    )
  }

  if (is.null(rownames(id_matrix)) || is.null(colnames(id_matrix)) ||
      !identical(rownames(id_matrix), colnames(id_matrix))) {
    cli::cli_abort(
      "`{arg}` should have identical row and column names given by the identifiers of `{parameter}`."
    )
  }

  if (!isSymmetric(unname(id_matrix))) {
    cli::cli_abort("`{arg}` should be a symmetric matrix.")
  }

  atomic <- unlist(lapply(pif_list, atomic_fractions), recursive = FALSE)
  ids    <- sapply(atomic, parameter_id, parameter = parameter)

  unknown <- setdiff(unique(rownames(id_matrix)), ids)
  if (length(unknown) > 0) {
    example <- unknown[1]
    cli::cli_abort(
      paste0(
        "The following identifiers in `{arg}` are not used by any of the ",
        "fractions: {.val {unknown}}. Give `{parameter}` as a named list, e.g. ",
        "{.code {parameter} = list(\"{example}\" = ...)}."
      )
    )
  }

  for (pif in atomic[!is.na(ids) & ids %in% rownames(id_matrix)]) {
    id    <- parameter_id(pif, parameter)
    n     <- length(parameter_value(pif, parameter))
    block <- id_covariance_block(list(id_matrix), id, id, n, n)
    if (!isTRUE(all.equal(block, unname(as.matrix(parameter_variance(pif, parameter)))))) {
      cli::cli_abort(
        paste0(
          "The variance of {.val {id}} in `{arg}` differs from the variance ",
          "`{arg}` given to fraction {.val {pif@label}}."
        )
      )
    }
  }

  id_matrix
}
