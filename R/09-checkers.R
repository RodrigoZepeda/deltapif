#' Function that checks the link
#'
#' Takes the link or the relative risk link of a `pif` or `paf` and checks
#' that it has been correctly specified.
#'
#' @inheritParams pifpaf
#'
#' @returns Invisible. Called for its side effects
#' @name linkcheck
#' @keywords internal
check_links <- function(link, link_deriv, link_inv){

  if (!is.function(link) && (is.null(link) || is.na(link))){
    cli::cli_abort(
      "No link was specified. See {.help pif} for details"
    )
  }

  if (!is.character(link) && !is.function(link)){
    cli::cli_abort(
      "Cannot handle link of type {typeof(link)}"
    )
  }

  if (!is.character(link_deriv) && !is.function(link_deriv) && !is.null(link_deriv) && !is.na(link_deriv)){
    cli::cli_abort(
      "Cannot handle link_deriv of type {typeof(link_deriv)}"
    )
  }

  if (!is.character(link_inv) && !is.function(link_inv) && !is.null(link_inv) && !is.na(link_inv)){
    cli::cli_abort(
      "Cannot handle link_inv of type {typeof(link_inv)}"
    )
  }

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

  if (is.function(link) & (!is.function(link_deriv) && !is.null(link_deriv) && !is.na(link_deriv))){
    cli::cli_abort(
      paste0(
        "A functional link was specified but no `link_deriv` was given. Set",
        "`link_deriv =` the derivative of your `link` function."
      )
    )
  }

  return(invisible())
}


#' @rdname linkcheck
#' @keywords internal
check_rr_links <- function(rr_link, rr_link_deriv){

  if (!is.function(rr_link) && (is.null(rr_link) || is.na(rr_link))){
    cli::cli_abort(
      "No rr_link was specified. See {.help pif} for details"
    )
  }

  if (!is.character(rr_link_deriv) && !is.function(rr_link_deriv) && !is.null(rr_link_deriv)  && !is.na(rr_link_deriv)){
    cli::cli_abort(
      "Cannot handle rr_link_deriv of type {typeof(rr_link_deriv)}"
    )
  }


  if (!is.character(rr_link) && !is.function(rr_link)){
    cli::cli_abort(
      "Cannot handle rr_link of type {typeof(rr_link)}"
    )
  }

  if (is.character(rr_link) && !is.null(rr_link_deriv)){
    cli::cli_abort(
      paste0(
        "A {rr_link} rr_link was specified but a `rr_link_deriv` was given. ",
        "Set `rr_link_deriv = NULL`."
      )
    )
  }

  if (is.function(rr_link) && (!is.function(rr_link_deriv) && !is.na(rr_link_deriv) && !is.null(rr_link_deriv))) {
    cli::cli_abort(
      "rr_link_deriv must be a function when rr_link is a function"
    )
  }

  return(invisible())
}
