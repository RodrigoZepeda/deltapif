#' Function that checks the link
#'
#' @inheritParams pifpaf
#'
#' @returns Invisible. Called for its side effects
#' @keywords internal
check_links <- function(link, link_deriv, link_inv){
  if (is.character(link) && !is.null(link_inv)){
    cli::cli_abort(
      paste0(
        "A {link} link was specified but a `link_inv` was given. ",
        "Set `link_inv = NULL`."
      )
    )
  }

  if (is.character(link) && !is.null(link_deriv)){
    cli::cli_abort(
      paste0(
        "A {link} link was specified but a `link_deriv` was given. ",
        "Set `link_deriv = NULL`."
      )
    )
  }

  if (is.function(link) & is.null(link_inv)){
    cli::cli_abort(
      paste0(
        "A functional link was specified but no `link_inv` was given. Set",
        "`link_inv =` to the inverse of your function."
      )
    )
  }

  if (is.function(link) & is.null(link_deriv)){
    cli::cli_abort(
      paste0(
        "A functional link was specified but no `link_deriv` was given. Set",
        "`link_deriv =` the derivative of your `link` function."
      )
    )
  }

  return(invisible())
}


#' Function that checks the relative risk link
#'
#' @inheritParams pifpaf
#'
#' @returns Invisible. Called for its side effects
#'
#' @keywords internal
check_rr_links <- function(rr_link, rr_link_deriv){


  if (is.character(rr_link) && !is.null(rr_link_deriv)){
    cli::cli_abort(
      paste0(
        "A {rr_link} rr_link was specified but a `rr_link_deriv` was given. ",
        "Set `rr_link_deriv = NULL`."
      )
    )
  }

  return(invisible())
}
