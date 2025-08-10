#' Exposure and Relative Risk data from Lee et Al
#'
#' Relative risk and exposure data for dementia risk-factors
#'
#' @format
#' A data frame with 12 rows and 11 columns:
#'
#' \describe{
#'   \item{risk_factor}{The risk factor for dementia}
#'   \item{RR}{The relative risk of dementia associated to the risk factor}
#'   \item{lower_CI, upper_CI}{Lower and upper bounds for the 95% confidence interval}
#'   \item{logrr}{The logarithm of the relative risk `log(RR)`}
#'   \item{sdlog}{The variance of the logarithm of the relative risk `variance(log(RR))`}
#'   \item{total}{The proportion of individuals exposed in the overall population}
#'   \item{hispanic}{The proportion of hispanic individuals exposed}
#'   \item{asian}{The proportion of non-hispanic asian individuals exposed}
#'   \item{black}{The proportion of non-hispanic black individuals exposed}
#'   \item{white}{The proportion of non-hispanic white individuals exposed}
#' }
#'
#' @references Lee, Mark, et al. "Variation in population attributable fraction
#' of dementia associated with potentially modifiable risk factors by race
#' and ethnicity in the US." JAMA network open 5.7 (2022): e2219672-e2219672.
#'
#' @seealso [dementiacov] for the covariance between the risk factors
"dementiarisk"

#' Relative risk covariance from Lee et al
#'
#' Covariances between the relative risks of Lee et al
#'
#' @format
#' A covariance matrix with 12 rows and 12 columns. Each entry
#' represents the covariance between them.
#'
#' @references Lee, Mark, et al. "Variation in population attributable fraction
#' of dementia associated with potentially modifiable risk factors by race
#' and ethnicity in the US." JAMA network open 5.7 (2022): e2219672-e2219672.
#'
#' @seealso [dementiarisk] for the relative risks and prevalence estimates
"dementiacov"

#' Alcohol consumption in Australia from Pandeya et al
#'
#' Alcohol consumption among Australian adults in grams/day
#'
#' @format
#' A data frame with 12 rows and 11 columns:
#'
#' \describe{
#'   \item{sex}{Whether individuals were male or female}
#'   \item{alcohol_g}{Category for the measure in grams of alcohol consumption}
#'   \item{median_intake}{The median intake for the group}
#'   \item{age_18_24, age_25_34, age_35_44, age_45_54, age_55_64, age_65_74, age_75_plus, age_18_plus}{Proportion of adults in each of the `alcohol_g` categories (by age group).}
#' }
#'
#' @references Pandeya, Nirmala, et al. "Cancers in Australia in 2010
#' attributable to the consumption of alcohol." Australian and New Zealand
#' journal of public health 39.5 (2015): 408-413.
#'
#' @seealso [dementiarisk] for the relative risks and prevalence estimates
"alcohol"

#' Relative risks for cancer from Pandeya et al
#'
#' Relative risks for cancer given different levels of alcohol consumption
#'
#' @format
#' A data frame with 24 rows and 5 columns:
#'
#' \describe{
#'   \item{cancer_type}{The type of cancer associated to the relative risk}
#'   \item{dose}{Dose in grams/day for which the relative risk was estimated}
#'   \item{RR}{Relative risk for cancer given the dose of alcohol}
#'   \item{lower_CI, upper_CI}{Lower and upper bounds for the 95% confidence interval}
#' }
#'
#' @references Pandeya, Nirmala, et al. "Cancers in Australia in 2010
#' attributable to the consumption of alcohol." Australian and New Zealand
#' journal of public health 39.5 (2015): 408-413.
#'
#' @seealso [dementiarisk] for the relative risks and prevalence estimates
"cancer_rr"
