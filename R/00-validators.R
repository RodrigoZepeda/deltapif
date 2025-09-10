#' Validators
#' Functions used in the validations of classes
#'
#' @inheritParams getters
#'
#' @name validators
#' @keywords internal
NULL

#' Validate any generic ensemble including pif total
#' @rdname validators
validate_global_ensemble <- function(self) {

  #Check they are pif class
  for (i in seq_along(self@pif_list)) {
    if (!S7::S7_inherits(self@pif_list[[i]], pif_class)) {
      cli::cli_abort(
        "Element {i} of `pif_list` must be a 'pif_class'."
      )
    }
  }

  #Check that they all gave names
  if (is.null(names(self@pif_list))){
    cli::cli_abort(
      "`pif_list` should be a named list of `pif_class` elements. No names were given."
    )
  }

  #Check the length of the weights
  if (length(self@weights) != length(self@pif_list)){
    cli::cli_abort(
      paste0(
        "weights provided have length {length(self@weights)} but ",
        "{length(self@pif_list)} fractions were provided."
      )
    )
  }
}


