#' Covariance structure class
#'
#' The covariance structure class represents a collection of
#' matrices and 0's such that `cov[[i]][[j]]` contains the
#' covariance of elements of the i-th potential impact fraction
#' and the j-th potential impact fraction.
#'
#' @param cov_list List such that `list[[i]][[j]]` represents the covariance
#' between elements of the i-th potential impact fraction
#' and the j-th potential impact fraction.
#'
#' @param dim_list List of dimensions required such that `dim_list[[i]][[j]]`
#' has a vector of 2 entries indicating the number of rows and columns
#' that can be inputed in cov_list (used for validation).
#'
#' @export
covariance_structure_class <- S7::new_class(
  name = "covariance_structure_class",
  package = "deltapif",
  properties = list(
    cov_list = S7::class_list
    #dim_list = S7::class_list
  ),
  validator = function(self){

    #Validate the size of cov_list
    n <- length(self@cov_list)

    #Check that the list has sublists
    if (any(sapply(self@cov_list, function(x) !is.list(x)))){
      cli::cli_abort(
        "Entry {which(sapply(self@cov_list,  function(x) !is.list(x)))} of `cov_list` is not a list."
      )
    }

    #Check that they have sublists of correct length
    if (any(sapply(self@cov_list, length) != length(self@cov_list[[1]]))){
      cli::cli_abort(
        paste0(
          "Element {which(sapply(self@cov_list, length) != length(self@cov_list[[1]]))} of",
          "`cov_list` has different length than {length(self@cov_list[[1]])}."
        )
      )
    }
  }
)


#' Default covariance structures
#'
#' These are multidimensional arrays of lists where
#' the entry `list[[i]][[j]]` exists represents the covariance
#' between elements of the i-th potential impact fraction
#' and the j-th potential impact fraction. This is particularly
#' useful when handling ensembles and totals.
#'
#' @param pif A potential impact fraction
#' @param pif1 A potential impact fraction to obtain a covariance structure with `pif2`.
#' @param pif2 A potential impact fraction to obtain a covariance structure with `pif1`,
#' @param sep Separation for the names in case `add_parents = TRUE`
#' @param add_parents add_parents Whether to add the parent impact fractions
#' to the name of the fraction.
#' @param parameter Either `beta` or `p`. Indicating which parameter
#' we are calculating covariance for.
#' @param parent_name For the ones ending in `2`. The name of the covariance
#' associated to `pif1` against `pif2`.
#'
#' @note The `covariance_structure`s ending in `2` are meant to
#' obtain the default covariance structure between two fractions `pif1` and `pif2`
#' while the ones that don't end in `2` are meant to obtain
#' the covariance structure of a fraction with itself.
#'
#' @return A nested list of lists with the entry `[[i]][[j]]` representing
#' the covariance between elements `i` and `j`.
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
#' #This is the default constructor of a covariance. Use it for custom covariances
#' covariance_structure(pif_lead_women)
#' covariance_structure2(pif_lead_women, pif_lead_men)
#' default_weight_covariance_structure2(pif_men, pif_men)
#' default_weight_covariance_structure(pif_tot)
#' default_weight_covariance_structure2(pif_men, pif_women)
#' default_parameter_covariance_structure(pif_tot, parameter = "beta")
#' default_parameter_covariance_structure2(pif_lead_women, pif_lead_men, parameter = "beta")
#' @name covariance_structures
NULL

#' @rdname covariance_structures
#' @export
covariance_structure <- function(pif, sep = "_>_", add_parents = FALSE, is_variance = FALSE) {

  if (S7::S7_inherits(pif, pif_atomic_class)){
    is_variance <- TRUE
  }

  names_of_list <- flatten_names(pif, sep = sep, add_parents = add_parents)

  #Check the names are unique
  if (!is_variance && length(names_of_list) != length(unique(names_of_list))){
    dups <- names_of_list[duplicated(names_of_list)]
    cli::cli_abort(
      paste0(
        "Some of the labels used for the fractions are not unique. ",
        "Labels {.val {dups}} were used multiple times."
      )
    )
  } else if (is_variance){
    names_of_list <- unique(names_of_list)
  }

  npifs         <- length(names_of_list)
  cov_str       <- vector("list", length = npifs)

  names(cov_str) <- names_of_list

  for (k in 1:npifs){
    cov_str[[k]] <- vector("list", length = npifs)
    names(cov_str[[k]]) <- names_of_list
    for (j in 1:npifs){
      cov_str[[k]][[j]] <- 0
    }
  }

  covariance_structure_class(cov_str)

}

#' @rdname covariance_structures
#' @export
covariance_structure2 <- function(pif1, pif2, sep = "_>_", add_parents = FALSE, parent_name = "Global"){

  #This is just a hack to apply the same function as before
  together_pif <- pif_ensemble(pif1, pif2, weights = c(0, 0), var_weights = 0, label = parent_name)

  is_variance <- identical(pif1, pif2)
  cov_str     <- covariance_structure(together_pif, sep = sep, add_parents = add_parents, is_variance = is_variance)

  #Eliminate the one that says global
  remove_name_from_covariance(cov_str, parent_name)

}

#' @rdname covariance_structures
#' @export
default_weight_covariance_structure <- function(pif, sep = "_>_", add_parents = FALSE,  is_variance = FALSE) {

  #Get the default covariance structure
  cov_str  <- covariance_structure(pif, sep = sep, add_parents = add_parents,  is_variance = is_variance)
  npifs    <- length(cov_str)
  cov_str  <- as.list(cov_str)

  #Get the flattened impact fractions
  flat_pif <- flatten(pif)

  if (npifs > 1){
    for (k in 1:npifs){
      for (j in k:npifs){

        #Case one of them is atomic then the weight covariance should be NULL
        if (S7::S7_inherits(flat_pif[[k]], pif_atomic_class) || S7::S7_inherits(flat_pif[[j]], pif_atomic_class)) {
          cov_str[[k]][[j]] <- 0

        #Case both of them are the same hence assign var_weights
        } else if (j == k){
          cov_str[[k]][[k]] <- flat_pif[[j]]@var_weights
        } else if (S7::S7_inherits(flat_pif[[k]], pif_global_ensemble_class) && S7::S7_inherits(flat_pif[[j]], pif_global_ensemble_class)) {

          #Compare the variances of the weights and the weights such that if they are the same they add to covariance:
          if (!is.null(flat_pif[[j]]@var_weights) && !is.null(flat_pif[[k]]@var_weights) &&
              !is.null(flat_pif[[j]]@weights) && !is.null(flat_pif[[k]]@weights) &&
              length(flat_pif[[j]]@var_weights) == length(flat_pif[[k]]@var_weights) &&
              length(flat_pif[[j]]@weights) == length(flat_pif[[k]]@weights) &&
              all(flat_pif[[j]]@var_weights == flat_pif[[k]]@var_weights) &&
              all(flat_pif[[j]]@weights == flat_pif[[k]]@weights)){

           cov_str[[k]][[j]] <- flat_pif[[j]]@var_weights

          } else {

            cov_str[[k]][[j]] <- 0

          }

        } else {
          cli::cli_abort(
            "Unknown class given to `default_weight_covariance_structure`."
          )
        }

        #Then fill the inverse value with the same for symmetry
        if (k != j){
          cov_str[[j]][[k]] <- cov_str[[k]][[j]]
        }
      }
    }
  }

  covariance_structure_class(cov_str)
}

#' @rdname covariance_structures
#' @export
default_weight_covariance_structure2 <- function(pif1, pif2, sep = "_>_", add_parents = FALSE, parent_name = "Global"){

  #This is just a hack to apply the same function as before
  together_pif <- pif_ensemble(pif1, pif2, weights = c(0, 0), var_weights = 0, label = parent_name)

  is_variance <- identical(pif1, pif2)
  cov_str <- default_weight_covariance_structure(together_pif, sep = sep, add_parents = add_parents,  is_variance = is_variance)

  #Eliminate the one that says global
  remove_name_from_covariance(cov_str, parent_name)

}

#' @rdname covariance_structures
#' @export
default_parameter_covariance_structure <- function(pif, sep = "_>_", add_parents = FALSE, parameter = "p", is_variance = FALSE) {

  if (parameter != "p" && parameter != "beta"){
    cli::cli_abort(
      "Invalid `parameter` = {parameter}. Set it to either `p` or `beta`."
    )
  }

  #Get the default covariance structure
  cov_str  <- covariance_structure(pif, sep = sep, add_parents = add_parents,  is_variance = is_variance)
  npifs    <- length(cov_str)

  #Get the flattened impact fractions
  flat_pif <- flatten(pif)


  if (length(flat_pif) == 1){
    flat_pif <- list(flat_pif)
  }


  for (k in 1:npifs){
    for (j in k:npifs){

      name_k <- flat_pif[[k]]@label
      name_j <- flat_pif[[j]]@label

      #By default we don't need those that have the beta's elsewhere
      if (j == k && S7::S7_inherits(flat_pif[[k]], pif_atomic_class) && S7::S7_inherits(flat_pif[[j]], pif_atomic_class)){

        cov_str[[name_k]][[name_k]] <- ifelse(parameter == "beta", flat_pif[[j]]@var_beta, flat_pif[[j]]@var_p)

      } else if (S7::S7_inherits(flat_pif[[k]], pif_atomic_class) && S7::S7_inherits(flat_pif[[j]], pif_atomic_class)) {

        #Compare the variances of the weights and the weights such that if they are the same they add to covariance:
        if (parameter == "beta"){

          if (!is.null(flat_pif[[j]]@var_beta) && !is.null(flat_pif[[k]]@var_beta) &&
              !is.null(flat_pif[[j]]@beta) && !is.null(flat_pif[[k]]@beta) &&
              length(flat_pif[[j]]@var_beta) == length(flat_pif[[k]]@var_beta) &&
              length(flat_pif[[j]]@beta) == length(flat_pif[[k]]@beta) &&
              all(flat_pif[[j]]@var_beta == flat_pif[[k]]@var_beta) &&
              all(flat_pif[[j]]@beta == flat_pif[[k]]@beta)){

            cov_str[[name_k]][[name_j]] <- flat_pif[[j]]@var_beta
            cov_str[[name_j]][[name_k]] <- flat_pif[[j]]@var_beta

          }

        } else if (parameter == "p") {

          if (!is.null(flat_pif[[j]]@var_p) && !is.null(flat_pif[[k]]@var_p) &&
              !is.null(flat_pif[[j]]@p) && !is.null(flat_pif[[k]]@p) &&
              length(flat_pif[[j]]@var_p) == length(flat_pif[[k]]@var_p) &&
              length(flat_pif[[j]]@p) == length(flat_pif[[k]]@p) &&
              all(flat_pif[[j]]@var_p == flat_pif[[k]]@var_p) &&
              all(flat_pif[[j]]@p == flat_pif[[k]]@p)){

            cov_str[[name_k]][[name_j]] <- flat_pif[[j]]@var_p
            cov_str[[name_j]][[name_k]] <- flat_pif[[j]]@var_p
          }
        }
      }
    }
  }


  #covariance_structure_class(cov_str)
  cov_str
}

#' @rdname covariance_structures
#' @export
default_parameter_covariance_structure2 <- function(pif1, pif2, sep = "_>_", add_parents = FALSE, parameter = "p", parent_name = "Global"){

  #This is just a hack to apply the same function as before
  together_pif <- pif_ensemble(pif1, pif2, weights = c(0, 0),
                               var_weights = 0, label = parent_name)

  is_variance <- identical(pif1, pif2)
  cov_str <- default_parameter_covariance_structure(together_pif, sep = sep, add_parents = add_parents, parameter = parameter,
                                                    is_variance = is_variance)

  #Eliminate the one that says global
  remove_name_from_covariance(cov_str, parent_name)

}

#' @rdname covariance_structures
#' @export
default_pif_covariance_structure <- function(pif, sep = "_>_", add_parents = FALSE, is_variance = FALSE) {


  #Get the default covariance structure
  cov_str  <- covariance_structure(pif, sep = sep, add_parents = add_parents,  is_variance = is_variance)
  npifs    <- length(cov_str)

  #Get the flattened impact fractions
  flat_pif <- flatten(pif)


  if (length(flat_pif) == 1){
    flat_pif <- list(flat_pif)
  }


  for (k in 1:npifs){
    name_k <- flat_pif[[k]]@label
    cov_str[[name_k]][[name_k]] <- flat_pif[[k]]@variance
  }

  cov_str
}

#' @rdname covariance_structures
#' @export
default_pif_covariance_structure2 <- function(pif1, pif2, sep = "_>_", add_parents = FALSE, parent_name = "Global"){

  #This is just a hack to apply the same function as before
  together_pif <- pif_ensemble(pif1, pif2, weights = c(0, 0),
                               var_weights = 0, label = parent_name)

  is_variance <- identical(pif1, pif2)
  cov_str <- default_pif_covariance_structure(together_pif, sep = sep, add_parents = add_parents,
                                                    is_variance = is_variance)

  #Eliminate the one that says global
  remove_name_from_covariance(cov_str, parent_name)

}

#' @rdname covariance_structures
#' @export
default_weight_pif_covariance_structure <- function(pif, sep = "_>_", add_parents = FALSE, is_variance = FALSE) {


  #Get the default covariance structure
  cov_str  <- covariance_structure(pif, sep = sep, add_parents = add_parents,  is_variance = is_variance)
  npifs    <- length(cov_str)

  #Get the flattened impact fractions
  flat_pif <- flatten(pif)


  if (length(flat_pif) == 1){
    flat_pif <- list(flat_pif)
  }


  for (k in 1:npifs){
    name_k <- flat_pif[[k]]@label
    if (S7::S7_inherits(flat_pif[[k]], pif_ensemble_class)){
      cov_str[[name_k]][[name_k]] <- flat_pif[[k]]@var_weights
    }
  }

  cov_str
}

#' @rdname covariance_structures
#' @export
default_weight_pif_covariance_structure2 <- function(pif1, pif2, sep = "_>_", add_parents = FALSE, parent_name = "Global"){

  #This is just a hack to apply the same function as before
  together_pif <- pif_ensemble(pif1, pif2, weights = c(0, 0),
                               var_weights = 0, label = parent_name)

  is_variance <- identical(pif1, pif2)
  cov_str <- default_weight_pif_covariance_structure(together_pif, sep = sep, add_parents = add_parents,
                                              is_variance = is_variance)

  #Eliminate the one that says global
  remove_name_from_covariance(cov_str, parent_name)

}

