# enable usage of <S7_object>@name in package code
#' @rawNamespace if (getRversion() < "4.3.0") importFrom("S7", "@")
NULL

.onLoad <- function(...) {
  S7::methods_register()

  S7::method(`[[`, covariance_structure_class) <- function(x, i, ...) {
    x@cov_list[[i]]
  }

  S7::method(`[[<-`, covariance_structure_class) <- function(x, i, ..., value) {
    if (!is.list(value)){
      cli::cli_abort(
        "Value to be assigned should be a list"
      )
    }

    if (length(value) != length(x)){
      cli::cli_abort(
        "Value to be assigned should be a list of length {length(x)}. But a value of length={length(value)} was given."
      )
    }
    x@cov_list[[i]] <- value
    x
  }



}
