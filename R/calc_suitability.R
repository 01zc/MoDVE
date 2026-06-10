
#' Calculate suitability score
#'
#' Given some environmental value, compute a suitability score between 0 and 1
#' indicating how well a species with the specified niche can develop in this condition.
#'
#' @param val the environmental value for which suitability is to be found
#' @param niche_min minimum value of the niche for this variable
#' @param niche_max maximum value of the niche for this variable
#' @param niche_opt optimal value of the niche for this variable
#'
#' @details
#' Suitability (\deqn{s}) is computed using the cardinal temperature response
#' introduced in Yan and Hunt (1999):
#'
#' \eqn{s = \frac{x_{max} - x}{x_{max} - x_{opt}}(\frac{x - x_{min}}{x_{opt} - x_{min}})^\frac{x_{opt} - x_{min}}{x_{max} - x_{opt}}}
#'
#' @references Weikei Yan, L. A. Hunt, An Equation for Modelling the Temperature
#' Response of Plants using only the Cardinal Temperatures, Annals of Botany,
#' Volume 84, Issue 5, November 1999, Pages 607–614,
#' https://doi.org/10.1006/anbo.1999.0955
#'
#' @export
calc_suitability <- function(val, niche_min, niche_max, niche_opt) {

  # Compute suitability only for valid entries
  num <- (niche_max - val) / (niche_max - niche_opt)
  denom <- (val - niche_min) / (niche_opt - niche_min)
  expo  <- (niche_opt - niche_min) / (niche_max - niche_opt)

  suitability <- num * denom ^ expo

  return(suitability)
}
