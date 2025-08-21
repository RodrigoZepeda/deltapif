#' Function to flatten a `pif_global_ensemble_class` and obtain the atomic pifs
#'
#' Returns a list of all the potential impact fractions
#' contained in a `pif` object.
#'
#' @param pif A potential impact fraction of either `pif_atomic_class`
#' or `pif_global_ensemble_class`
#'
#'
#' @return A list of `pif_atomic_class` elements
#'
#' @examples
#' #Potential impact fraction for women
#' pif_lead_women <- paf(0.27, 2.2, quiet = TRUE, var_p = 0.001)
#' pif_rad_women  <- paf(0.12, 1.2, quiet = TRUE, var_p = 0.001)
#' pif_women      <- pif_ensemble(pif_lead_women, pif_rad_women)
#' pif_lead_men   <- paf(0.30, 2.2, quiet = TRUE, var_p = 0.001)
#' pif_rad_men    <- paf(0.10, 1.2, quiet = TRUE, var_p = 0.001)
#' pif_men        <- pif_ensemble(pif_lead_men, pif_rad_men)
#'
#' flatten(pif_women)
#'
#' pif_tot <- pif_total(pif_men, pif_women, weights = c(0.49, 0.51))
#'
#' flatten(pif_tot)
#'
#' @export
flatten <- function(pif){

  if (S7::S7_inherits(pif, pif_atomic_class)){
    pifvec <- list(pif)
    names(pifvec) <- pif@label
    return(pifvec)
  }

  if (!S7::S7_inherits(pif, pif_global_ensemble_class)){
    cli::cli_abort(
      "Invalid `pif` is not a `pif_global_ensemble_class` not a `pif_atomic_class`."
    )
  }

  #Assuming its a pif_global_ensemble_class
  pif_flattened_list <- list(pif)
  names(pif_flattened_list) <- pif@label
  for (k in 1:length(pif@pif_list)){
    pif_flattened_list <- pif_flattened_list |>
      append(flatten(pif@pif_list[[k]]))
  }

  return(pif_flattened_list)

}

#' Function to flatten a pif's names
#'
#' Returns a character vector such that the names are
#' `parent_child_grandchild`.
#'
#' @param pif A potential impact fraction of either `pif_atomic_class`
#' or `pif_global_ensemble_class`
#'
#'
#' @return A character vector
#'
#' @examples
#' pif_lead_women <- paf(0.27, 2.2, quiet = TRUE, var_p = 0.001, label = "Women lead")
#' pif_rad_women  <- paf(0.12, 1.2, quiet = TRUE, var_p = 0.001, label = "Women radiation")
#' pif_women      <- pif_ensemble(pif_lead_women, pif_rad_women, label = "Women")
#' pif_lead_men   <- paf(0.30, 2.2, quiet = TRUE, var_p = 0.001, label = "Men lead")
#' pif_rad_men    <- paf(0.10, 1.2, quiet = TRUE, var_p = 0.001, label = "Men radiation")
#' pif_men        <- pif_ensemble(pif_lead_men, pif_rad_men, label = "Men")
#' pif_tot <- pif_total(pif_men, pif_women, weights = c(0.49, 0.51), label = "Population")
#'
#' flatten_names(pif_lead_women)
#' flatten_names(pif_tot)
#'
#' @export
flatten_names <- function(pif, sep = "_>_", add_parents = TRUE){

  if (S7::S7_inherits(pif, pif_atomic_class)){
    return(pif@label)
  } else {
    return(
      c(pif@label, unlist(sapply(pif@pif_list, flatten_names)))
    )
  }
}



#' Removes the name parent_name from the cov_str list
#'
#' The `cov_str` is a list of covariances
#'
#' @keywords internal
remove_name_from_covariance <- function(covariance_str, parent_name){

  if (!S7::S7_inherits(covariance_str, covariance_structure_class)){
    cli::cli_abort(
      "Object `covariance_str` should be a covariance structure of `covariance_structure_class`"
    )
  }

  cov_str <- covariance_str@cov_list
  k <- 1
  while (k <= length(cov_str)){
    if (names(cov_str)[k] == parent_name){
      cov_str[[k]] <- NULL
    } else {
      j <- 1
      while (j <= length(cov_str[[k]])){
        if (names(cov_str[[k]])[j] == parent_name){
          cov_str[[k]][[j]] <- NULL
        } else {
          j <- j + 1
        }
      }
      k <- k + 1
    }
  }

  covariance_str@cov_list <- cov_str
  return(covariance_str)
}


