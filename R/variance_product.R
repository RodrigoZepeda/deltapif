#' Variance product formula
#'
#' Implements the variance product formula for the product of two random variables.
#' This is the formula for `var(X*Y)` given by:
#' \deqn{
#' \text{Var}[XY] = Y^2\cdot\text{Var}[X]  + X^2\cdot \text{Var}[Y]
#' + \text{Var}[X] \cdot \text{Var}[Y]
#' }
#'
#' @param x Numeric value representing the mean of the first random variable X.
#' @param y Numeric value representing the mean of the second random variable Y.
#' @param var_x Numeric value representing the variance of X. Must be non-negative.
#' @param var_y Numeric value representing the variance of Y. Must be non-negative.
#'
#' @return Numeric value representing the variance of the product X*Y.
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' var_prod(x = 2, y = 3, var_x = 0.5, var_y = 0.8)
#'
#' # When one variable has zero variance
#' var_prod(x = 5, y = 4, var_x = 0, var_y = 1.2)
#'
#' # When both variables have zero variance
#' var_prod(x = 10, y = 20, var_x = 0, var_y = 0)
#'
#' # With larger variances
#' var_prod(x = 100, y = 200, var_x = 25, var_y = 36)
#' }
#' @keywords internal
var_prod <- function(x, y, var_x, var_y) {

  # Input validation
  if (!is.numeric(x) || length(x) != 1) {
    cli::cli_abort("x must be a single numeric value")
  }

  if (!is.numeric(y) || length(y) != 1) {
    cli::cli_abort("y must be a single numeric value")
  }

  if (!is.numeric(var_x) || length(var_x) != 1 || var_x < 0) {
    cli::cli_abort("var_x must be a single non-negative numeric value")
  }

  if (!is.numeric(var_y) || length(var_y) != 1 || var_y < 0) {
    cli::cli_abort("var_y must be a single non-negative numeric value")
  }

  # Handle infinite values
  if (any(is.infinite(c(x, y, var_x, var_y)))) {
    cli::cli_abort("Inputs cannot be infinite")
  }

  # Handle NA/NaN values
  if (any(is.na(c(x, y, var_x, var_y))) || any(is.nan(c(x, y, var_x, var_y)))) {
    cli::cli_abort("Inputs cannot be NA or NaN")
  }

  # Calculate variance using the product formula
  result <- (y^2 * var_x) + (x^2 * var_y) + (var_x * var_y)

  # Ensure non-negative result (theoretically should be, but check for numerical precision)
  if (result < 0) {
    cli::cli_alert_warning("Variance result is negative due to numerical precision, returning absolute value")
    result <- abs(result)
  }

  return(result)
}
