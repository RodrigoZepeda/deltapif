#' Setters
#'
#' Collection of functions for the potential impact fraction class
#' so that they work as setters of the properties. Each
#' function is constructed as `set_property`
#'
#' @param self A `pif_class` object created with S7.
#'
#' @name setters
#'
#' @keywords internal
NULL

set_conf_level <- function(self, value){
  if (!is.numeric(value) | length(value) != 1 | value < 0 | value > 1){
    cli::cli_abort(
      "Invalid value {value} given to confidence level."
    )
  }
  self@conf_level <- value
  self
}
