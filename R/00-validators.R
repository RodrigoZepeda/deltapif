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

  for (i in seq_along(self@pif_list)) {
    if (!S7::S7_inherits(self@pif_list[[i]], pif_class)) {
      cli::cli_abort(
        "Element {i} of `pif_list` must be a 'pif_class'."
      )
    }
  }

  if (length(self@weights) != length(self@pif_list)){
    cli::cli_abort(
      paste0(
        "weights provided have length {length(self@weights)} but ",
        "{length(self@pif_list)} fractions were provided."
      )
    )
  }
}
