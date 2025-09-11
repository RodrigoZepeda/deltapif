#' Variance product formula
#'
#' Implements the variance product formula. This is the formula for
#' `var(X*Y)` given by:
#' \deqn{
#' \text{Var}[XY] = Y^2\cdot\text{Var}[X]  + X^2\cdot \text{Var}[Y]
#' + \text{Var}[X] \cdot \text{Var}[Y]
#' }
#'
#' @keywords internal
var_prod <- function(x, y, var_x, var_y){
  x^2*var_y + var_y*x^2 + var_x*var_y
}
