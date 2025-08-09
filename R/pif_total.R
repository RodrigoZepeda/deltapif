#' Total PIF
#' @examples
#' beta <- 1.3
#' sigma_beta <- 0.1
#' pif1 <- pif(0.5, 0.2, beta, sigma_p = 0.5 * (1 - 0.5) / 100, sigma_beta = sigma_beta)
#' pif2 <- pif(0.3, 0.1, beta, sigma_p = 0.3 * (1 - 0.3) / 100, sigma_beta = sigma_beta)
#' pif3 <- pif(0.7, 0.3, beta, sigma_p = 0.7 * (1 - 0.7) / 100, sigma_beta = sigma_beta)
#' pif_total(pif1, pif2, pif3, q = c(0.3, 0.2, 0.5))
#' @export
pif_total <- function(pif1, ..., q, sigma_q){

  #Get covariance matrix
  sigma_pif <- cov(pif1, ...)

  pif_list_class <- append(list(pif1), list(...))
  npifs    <- length(pif_list_class)
  coefs    <- sapply(pif_list_class, coef)


  as.numeric(t(q) %*% coefs)

}

