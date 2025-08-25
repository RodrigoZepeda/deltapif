#' Print a `pif_class`
#'
#' Prints a `pif_class` object.
#'
#' @param x A `pif_class`
#'
#' @param accuracy The accuracy parameter for [`scales::percent`].
#'
#' @return Called for its side-effects of printing to the console
#'
#' @keywords internal
print_pif_class <- function(x, accuracy){

  # Printed text looks like:
  #   ── Potential Impact Fraction ──
  #
  # PIF = 4.421% [95% CI: 0.179% to 54.364%]
  # standard_deviation(pif %) = 7.003
  # standard_deviation(link(pif)) = 1.657

  pif_val <- scales::percent(x@pif, accuracy = accuracy)
  cilow   <- scales::percent(x@ci[1], accuracy = accuracy)
  cihigh  <- scales::percent(x@ci[2], accuracy = accuracy)
  title   <- ifelse(x@type == "PIF", "Potential Impact Fraction", "Population Attributable Fraction")

  cli::cli_h3("{title}: {.val {x@label}}")
  cli::cli_text(
    "{x@type} = {pif_val} ",
    "[{.emph {scales::percent(x@conf_level)} CI}: {cilow} to {cihigh}]"
  )
  # cli::cli_text(
  #   "standard_deviation({tolower(x@type)} %) = {scales::comma(100*sqrt(x@variance), accuracy = accuracy)}"
  # )
  # cli::cli_text(
  #   "standard_deviation(link({tolower(x@type)})) = {scales::comma(sqrt(x@link_variance), accuracy = accuracy)}"
  # )

  return(invisible())
}


#' Print or show a potential impact fraction
#'
#' Function to print or show a potential impact fraction object
#'
#' @param x A `pif_class`
#'
#' @param ... Additional arguments to pass to `print` or `show`.
#'
#' @param accuracy The accuracy of the printed value
#'
#' @examples
#' my_pif <- pif(p = 0.2, beta = 1.3, var_beta = 0.1)
#' print(my_pif)
#'
#' # Change the ammount of digits to show just 1
#' print(my_pif, accuracy = 0.1)
#' @name print
#' @export
S7::method(print, pif_class) <- function(x, ..., accuracy = 0.001) {
  print_pif_class(x, accuracy)
  # cli::cli_h3("Parameters:")
  # cli::cli_ul()
  # cli::cli_li("Observed prevalence ({.code p_obs}) = {x@p}")
  # cli::cli_li("Counterfactual prevalence ({.code p_cft}) = {x@p_cft}")
  # cli::cli_li("Relative risks ({.code rr}) = {x@rr}")
  # cli::cli_end()
}

#' Print or show a covariance structure class
#'
#' Function to print or show a `covariance_structure_class`
#'
#' @param x A `covariance_structure_class`
#'
#' @param ... Additional arguments to pass to `print`
#'
#' @examples
#' pif_lead_women <- paf(0.27, 2.2, quiet = TRUE, var_p = 0.001, var_beta = 0.015,
#'                       label = "Women lead")
#' pif_rad_women  <- paf(0.12, 1.2, quiet = TRUE, var_p = 0.001, var_beta = 0.022,
#'                       label = "Women radiation")
#' pif_women      <- pif_ensemble(pif_lead_women, pif_rad_women, label = "Women",
#'                                weights = c(0.8, 0.72),
#'                                var_weights = matrix(c(0.3, 0.1, 0.1, 0.4), ncol = 2))
#'
#' pif_lead_men   <- paf(0.30, 2.2, quiet = TRUE, var_p = 0.001, var_beta = 0.015,
#'                        label = "Men lead")
#' pif_rad_men    <- paf(0.10, 1.2, quiet = TRUE, var_p = 0.001, var_beta = 0.022,
#'                         label = "Men radiation")
#' pif_men        <- pif_ensemble(pif_lead_men, pif_rad_men, label = "Men",
#'                         weights = c(0.65, 0.68),
#'                         var_weights = matrix(c(0.1, -0.2, -0.2, 0.5), ncol = 2))
#' pif_tot        <- pif_total(pif_men, pif_women,
#'                         weights = c(0.49, 0.51), label = "Population",
#'                         var_weights = matrix(c(0.22, 0.4, 0.4, 0.8), ncol = 2))
#'
#' print(covariance_structure(pif_lead_women))
#' print(covariance_structure2(pif_lead_women, pif_lead_men))
#' print(default_weight_covariance_structure2(pif_men, pif_women))
#' print(default_parameter_covariance_structure(pif_tot, parameter = "beta"))
#' @name print
#' @export
S7::method(print, covariance_structure_class) <- function(x, ..., quote = FALSE) {
  ndim <- length(x@cov_list)
  mdim <- length(x@cov_list[[1]])

  mat <- matrix(".", nrow = ndim, ncol = mdim)
  rownames(mat) <- names(x@cov_list)
  colnames(mat) <- names(x@cov_list[[1]])
  for (k in 1:ndim){
    for (j in 1:mdim){
      if (is.matrix(x@cov_list[[k]][[j]]) && (ncol(x@cov_list[[k]][[j]]) > 1 || nrow(x@cov_list[[k]][[j]]) > 1)){
        mat[k,j] <- paste0(nrow(x@cov_list[[k]][[j]]), "x", ncol(x@cov_list[[k]][[j]]))
      } else if (is.vector(x@cov_list[[k]][[j]]) && length(x@cov_list[[k]][[j]]) > 1){
        mat[k,j] <- paste0(1, "x", length(x@cov_list[[k]][[j]]))
      } else if (is.numeric(x@cov_list[[k]][[j]]) && as.double(x@cov_list[[k]][[j]]) != 0){
        mat[k,j] <- as.double(x@cov_list[[k]][[j]])
      }
    }
  }
  print(mat, ..., quote = quote)
}



#' Extract coefficients of a pif object
#'
#' Gets the potential impact fraction value
#'
#' @param object A `pif_class` object.
#' @param ... Additional parameters to pass to `coef` (ignored)
#'
#' @examples
#' my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, var_p = 0.1, var_beta = 0.2)
#' coef(my_pif)
#' @name coef
#' @export
S7::method(coef, pif_class) <- function(object, ...) {
  object@pif
}

#' Extract confidence intervals of a pif object
#'
#' Gets the confidence interval for the potential impact fraction
#'
#' @param object A `pif_class` object.
#' @param level Level of confidence desired.
#' @param ... Additional parameters to pass to `confint` (ignored)
#'
#' @examples
#' my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, var_p = 0.1, var_beta = 0.2)
#' #Default 95% CI
#' confint(my_pif)
#'
#' #Custom 90% ci:
#' confint(my_pif, level = 0.90)
#' @name confint
#' @export
S7::method(confint, pif_class) <- function(object, ..., level = object@conf_level) {
  #Set the level
  object@conf_level <- level
  return(object@ci)

}

#' Extract weights of a pif_global_ensemble
#'
#' Gets the weights of a `pif_global_ensemble` object
#'
#' @param object A `pif_global_ensemble` object.
#' @param ... Additional parameters to pass to `weights` (ignored)
#'
#' @examples
#' my_pif1 <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, var_p = 0.1, var_beta = 0.2)
#' my_pif2 <- pif(p = 0.3, p_cft = 0.1, beta = 1.5, var_p = 0.1, var_beta = 0.2)
#' my_pif  <- pif_total(my_pif1, my_pif2, weights = c(0.8, 0.2))
#' weights(my_pif)
#'
#' @name weights
#' @export
S7::method(weights, pif_global_ensemble_class) <- function(object, ...) {
  return(object@weights)
}

#' Summary of a pif object
#'
#' Gets the potential impact fraction summary
#'
#' @param object A `pif_class` object.
#' @param ... Additional parameters to pass to `summary` (ignored)
#'
#' @examples
#' my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, var_p = 0.1, var_beta = 0.2)
#' summary(my_pif)
#' @name summary
#' @export
S7::method(summary, pif_class) <- function(object, level = object@conf_level, ...) {
  conf_interval <- confint(object, level = level)

  #Build the return vector
  return_vec <- c("value"      = coef(object),
                  "standard_deviation" = standard_deviation(object),
                  "ci_low"     = conf_interval[1],
                  "ci_up"      = conf_interval[2],
                  "confidence" = level)

  #Assign the name
  names(return_vec)[1] <- fraction_type(object)

  return(return_vec)
}

#' Transform a pif object into a data.frame
#'
#' Gets the potential impact fraction value, the link_variance and the confidence
#' interval values
#'
#' @param x A `pif_class` object.
#' @param ... Additional parameters (ignored)
#'
#' @examples
#' my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, var_p = 0.1, var_beta = 0.2)
#' as.data.frame(my_pif)
#' @name as.data.frame
#' @export
S7::method(as.data.frame, pif_class) <- function(x, ..., level = 0.95) {
  as.data.frame(t(summary(x, level = level)))
}

#' Get the label of a PIF or PAF
#'
#' Gets the label of a potential impact fraction or a population
#' attributable fraction
#'
#' @param x A `pif_class` object.
#' @param ... Additional parameters (ignored)
#'
#' @examples
#' #A simple example
#' my_pif <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, var_p = 0.1,
#'                 var_beta = 0.2, label = "Test")
#' names(my_pif)
#'
#' #A pif composed of others
#' my_pif1 <- pif(p = 0.5, p_cft = 0.25, beta = 1.3, var_p = 0.1,
#'                 var_beta = 0.2, label = "Test 1")
#' my_pif2 <- pif(p = 0.4, p_cft = 0.1, beta = 1.3, var_p = 0.1,
#'                 var_beta = 0.2, label = "Test 2")
#' pif_tot <- pif_total(my_pif1, my_pif2, weights = c(0.2, 0.8),
#'                 label = "Parent")
#' names(pif_tot)
#'
#' @name names
NULL

#' @name names
#' @export
S7::method(names, pif_class) <- function(x, ...) {
  x@label
}


#' @name names
#' @export
S7::method(names, pif_global_ensemble_class) <- function(x, ...) {
  return(children(x))
}


#' Length of a `covariance_structure`
#'
#' Gets the length of a covariance structure
#'
#' @param x A `covariance_structure`
#'
#' @name length
#' @export
S7::method(length, covariance_structure_class) <- function(x) {
  length(x@cov_list)
}

#' Length of a `pif` ensemble
#'
#' Gets the length of a `pif_ensemble_class`
#'
#' @param x A `pif_ensemble_class`
#'
#' @name length
#' @export
S7::method(length, pif_global_ensemble_class) <- function(x) {
  length(x@pif_list)
}

#' Length of a `pif` ensemble
#'
#' Gets the length of a `pif_atomic_class`
#'
#' @param x A `pif_atomic_class`
#'
#' @name length
#' @export
S7::method(length, pif_atomic_class) <- function(x) {
  1
}

#' Convert a `covariance_structure` to `matrix`
#'
#' Transforms a `covariance_structure` into a `matrix`.
#'
#' @param x A `covariance_structure`
#' @param ... Additional parameters (ignored)
#'
#' @name as.matrix
#'
#' @examples
#' as.matrix(covariance_structure_class(list(b = list(a = 1:3))))
#'
#'
#' @export
S7::method(as.matrix, covariance_structure_class) <- function(x, ...) {
  flag <- TRUE
  nrow <- length(x)
  ncol <- length(x@cov_list[[1]])
  mat  <- matrix(0, nrow = nrow, ncol = ncol)
  for (k in 1:ncol){
    for (j in 1:nrow){
      if (!is.matrix(x@cov_list[[j]][[k]]) && !is.vector(x@cov_list[[j]][[k]])){
        mat[j,k] <- x@cov_list[[j]][[k]]
      } else if (is.matrix(x@cov_list[[j]][[k]]) && ncol(x@cov_list[[j]][[k]]) == 1 && nrow(x@cov_list[[j]][[k]] == 1)){
        mat[j,k] <- as.numeric(x@cov_list[[j]][[k]])
      } else if (is.vector(x@cov_list[[j]][[k]]) && length(x@cov_list[[j]][[k]]) == 1){
        mat[j,k] <- as.numeric(x@cov_list[[j]][[k]])
      } else {
        flag <- TRUE
      }
    }
  }

  if (!flag){
    return(mat)
  }

  #Go through each of the columns and get the maximum column
  max_col <- 0
  for (k in 1:ncol){
    for (j in 1:nrow){
      if (is.matrix(x@cov_list[[j]][[k]])){
        max_col <- max(ncol(x@cov_list[[j]][[k]]), max_col)
      } else if (is.vector(x@cov_list[[j]][[k]])){ #Assume vectors are column-vectors
        max_col <- max(1, max_col)
      }
    }
  }
  #Go through each of the rows and get the maximum number of rows
  max_row <- 0
  for (k in 1:ncol){
    for (j in 1:nrow){
      if (is.matrix(x@cov_list[[j]][[k]])){
        max_row <- max(nrow(x@cov_list[[j]][[k]]), max_row)
      } else if (is.vector(x@cov_list[[j]][[k]])){ #Assume vectors are column-vectors
        max_row <- max(length(x@cov_list[[j]][[k]]), max_row)
      }
    }
  }

  #Create a matrix of that size
  mat <- matrix(0, ncol = max_col*ncol, nrow = max_row*nrow)
  for (k in 1:ncol){
    for (j in 1:nrow){
      init_row <- ((j - 1)*max_row + 1)
      init_col <- ((k - 1)*max_col + 1)
      if (is.matrix(x@cov_list[[j]][[k]])){
        mat[
          init_row:(init_row + nrow(x@cov_list[[j]][[k]]) - 1),
          init_col:(init_col + ncol(x@cov_list[[j]][[k]]) - 1)
        ] <- x@cov_list[[j]][[k]]
      } else if (is.vector(x@cov_list[[j]][[k]])){
        mat[
          init_col:(init_col + length(x@cov_list[[j]][[k]]) - 1),
          1
        ] <- x@cov_list[[j]][[k]]
      }
    }
  }

  return(mat)

}

#' Subset a `covariance_structure`
#'
#' Obtains a smaller `covariance_structure` with the entries
#' given by the `select` option as a vector.
#'
#' @param x A `covariance_structure`
#' @param select A vector of covariates to keep
#' @param cols A vector of column names to keep
#' @param rows A vector of row names to keep
#' @param ... Additional parameters (ignored)
#'
#' @examples
#' pif_lead_women <- paf(0.27, 2.2, quiet = TRUE, var_p = 0.001, var_beta = 0.015,
#'                       label = "Women lead")
#' pif_rad_women  <- paf(0.12, 1.2, quiet = TRUE, var_p = 0.001, var_beta = 0.022,
#'                       label = "Women radiation")
#' covstr <- default_parameter_covariance_structure2(pif_lead_women, pif_rad_women, parameter = "beta")
#' subset(covstr, "Women lead")
#' subset(covstr, c("Women radiation", "Women lead"))
#' subset(covstr, 2)
#' subset(covstr, 1:2)
#'
#' @name subset
#' @export
S7::method(subset, covariance_structure_class) <- function(x, select = NULL, cols = NULL, rows = NULL, ...) {

  if (is.null(select) && is.null(cols) && is.null(rows)){
    return(x)
  }

  if (!is.null(select) && is.null(cols) && is.null(rows)){

    if (is.numeric(select) && max(select) > length(x) || min(select) < 0){
      cli::cli_abort(
        "Cannot `select` values outside the range `1:{length(x)}`"
      )
    }

    if (is.numeric(select) && any(round(select) != select)){
      cli::cli_abort(
        "Cannot select non integer numeric values"
      )
    }

    #Get the names
    names_cov_str <- c()
    if (is.character(select)){
      names_cov_str <- names(x@cov_list)[which(names(x@cov_list) %in% select)]
    } else if (is.numeric(select)){
      names_cov_str <- names(x@cov_list)[select]
    } else {
      cli::cli_abort(
        "Parameter `select` of `subset` should be either numeric or character."
      )
    }

    if (length(names_cov_str) < 1){
      cli::cli_abort(
        "Subset to `select` not found."
      )
    }

    #Create a covariance structure
    cov_str <- vector("list", length = length(names_cov_str))
    names(cov_str) <- names_cov_str
    for (k in 1:length(cov_str)){
      cov_str[[k]] <- vector("list", length = length(names_cov_str))
      names(cov_str[[k]]) <- names_cov_str
    }

    #Loop using numbers
    for (var in names_cov_str){
      for (var2 in names_cov_str){
        cov_str[[var]][[var2]] <- x@cov_list[[var]][[var2]]
      }
    }

    #Add to class
    return(covariance_structure_class(cov_str))

  } else if (!is.null(select) && (!is.null(cols) || !is.null(rows))){
    cli::cli_abort(
      "Don't specify `select` at the same time as `cols` or `rows` as I don't know how to proceed."
    )
  } else if (is.null(select)){

    cov_str <- x

    if (!is.null(rows)){
      cov_str <- subset_row(cov_str, select = rows)
    }

    if (!is.null(cols)){
      cov_str <- subset_col(cov_str, select = cols)
    }

    return(cov_str)

  }
}


#' Function to get the rows out of a `covariance_structure`
#' @keywords internal
subset_col <- function(x, select, ...) {

  if (!S7::S7_inherits(x, covariance_structure_class)){
    cli::cli_abort(
      "Currently `subset_col` only supports `covariance_structures`"
    )
  }

  if (!is.vector(select)){
    cli::cli_abort(
      "`select` must be a vector of length > 0"
    )
  }

  if (length(select) < 1){
    return(x)
  }

  if (is.numeric(select) && max(select) > length(x) || min(select) < 0){
    cli::cli_abort(
      "Cannot `select` values outside the range `1:{length(x)}`"
    )
  }

  if (is.numeric(select) && any(round(select) != select)){
    cli::cli_abort(
      "Cannot select non integer numeric values"
    )
  }

  #Get the names
  names_cov_str <- c()
  if (is.character(select)){
    names_cov_str <- names(x@cov_list[[1]])[which(names(x@cov_list[[1]]) %in% select)]
  } else if (is.numeric(select)){
    names_cov_str <- names(x@cov_list[[1]])[select]
  } else {
    cli::cli_abort(
      "Parameter `select` of `subset` should be either numeric or character."
    )
  }
  names_unselected <- names(x@cov_list)


  if (length(names_cov_str) < 1){
    cli::cli_abort(
      "Subset to `select` not found."
    )
  }

  #Create a covariance structure
  cov_str <- vector("list", length = length(names_unselected))
  names(cov_str) <- names_unselected
  for (k in 1:length(cov_str)){
    cov_str[[k]] <- vector("list", length = length(names_cov_str))
    names(cov_str[[k]]) <- names_cov_str
  }

  #Loop using numbers
  for (var in names_unselected){
    for (var2 in names_cov_str){
      cov_str[[var]][[var2]] <- x@cov_list[[var]][[var2]]
    }
  }

  #Add to class
  covariance_structure_class(cov_str)

}

#' Function to get the rows out of a `covariance_structure`
#' @keywords internal
subset_row <- function(x, select, ...) {

  if (!S7::S7_inherits(x, covariance_structure_class)){
    cli::cli_abort(
      "Currently `subset_col` only supports `covariance_structures`"
    )
  }

  if (!is.vector(select)){
    cli::cli_abort(
      "`select` must be a vector of length > 0"
    )
  }

  if (length(select) < 1){
    return(x)
  }

  if (is.numeric(select) && max(select) > length(x) || min(select) < 0){
    cli::cli_abort(
      "Cannot `select` values outside the range `1:{length(x)}`"
    )
  }

  if (is.numeric(select) && any(round(select) != select)){
    cli::cli_abort(
      "Cannot select non integer numeric values"
    )
  }

  #Get the names
  names_cov_str <- c()
  if (is.character(select)){
    names_cov_str <- names(x@cov_list)[which(names(x@cov_list) %in% select)]
  } else if (is.numeric(select)){
    names_cov_str <- names(x@cov_list)[select]
  } else {
    cli::cli_abort(
      "Parameter `select` of `subset` should be either numeric or character."
    )
  }
  names_unselected <- names(x@cov_list[[1]])


  if (length(names_cov_str) < 1){
    cli::cli_abort(
      "Subset to `select` not found."
    )
  }

  #Create a covariance structure
  cov_str <- vector("list", length = length(names_cov_str))
  names(cov_str) <- names_cov_str
  for (k in 1:length(cov_str)){
    cov_str[[k]] <- vector("list", length = length(names_unselected))
    names(cov_str[[k]]) <- names_unselected
  }

  #Loop using numbers
  for (var in names_cov_str){
    for (var2 in names_unselected){
      cov_str[[var]][[var2]] <- x@cov_list[[var]][[var2]]
    }
  }

  #Add to class
  covariance_structure_class(cov_str)

}


#' Convert a `covariance_structure` to `list`
#'
#' Transforms a `covariance_structure` into a `list`.
#'
#' @param x A `covariance_structure`
#' @param ... Additional parameters (ignored)
#'
#' @name as.list
#' @export
S7::method(as.list, covariance_structure_class) <- function(x, ...) {
  x@cov_list
}

#' Get the children labels from a `pif_class`
#'
#' Gets the labels of the elemebts in the `pif_list` of an ensemble class.
#'
#' @param x A `pif_global_ensemble_class`
#' @param ... Additional arguments (currently ignored)
#' @name children
NULL

#' @rdname children
#' @export
children <- S7::new_generic(
  "children", "x",
  function(x, ...) {
    S7::S7_dispatch()
  }
)
S7::method(children, pif_global_ensemble_class) <- function(x, ...) {
  sapply(x@pif_list, function(x) x@label)
}
S7::method(children, pif_atomic_class) <- function(x, ...) {
  NULL
}


#' Transform into a covariance structure
#'
#' Transforms either a matrix or a vector into a covariance structure.
#'
#' @param x A matrix or a vector
#' @param col_names Names to assign to the columns
#' @param row_names Names to assign to the rows
#' @param ... Additional arguments (currently ignored)
#' @name as_covstr
#'
#' @examples
#' mat <- matrix(c(1,3,2,4), ncol = 2,
#'           dimnames = list(list("pif1", "pif2"), list("pif1", "pif2")))
#' as_covariance_structure(mat)
#' @export

#' @rdname as_covstr
#' @export
as_covariance_structure <- function(x, col_names = NULL, row_names = NULL, ...) {

  #Already a covariance structure
  if (S7::S7_inherits(x, covariance_structure_class)){
    return(x)
  }

  if (!is.matrix(x) && !S7::S7_inherits(x, covariance_structure_class)){
    cli::cli_abort(
      "Method is only available for objects of class matrix"
    )
  }

  if (is.null(col_names)){
    if (is.null(colnames(x))){
      col_names <- paste0(1:ncol(x))
    } else {
      col_names <- colnames(x)
    }
  } else {
    if (length(col_names) != ncol(x)){
      cli::cli_abort(
        "{length(col_names)} col_names were given for a matrix of {ncol(x)} columns."
      )
    }
  }

  if (is.null(row_names)){
    if (is.null(rownames(x))){
      row_names <- paste0(1:nrow(x))
    } else {
      row_names <- rownames(x)
    }
  } else {
    if (length(row_names) != nrow(x)){
      cli::cli_abort(
        "{length(row_names)} row_names were given for a matrix of {nrow(x)} rows"
      )
    }
  }

  #Create the empty structure
  row_list <- vector("list", nrow(x))
  names(row_list) <- row_names
  for (i in 1:nrow(x)){
    col_list <- vector("list", ncol(x))
    names(col_list) <- col_names

    #Loop assigning the value
    for (j in 1:ncol(x)){
      col_list[[j]] <- x[i,j]
    }

    row_list[[i]] <- col_list
  }

  covariance_structure_class(cov_list = row_list)
}




