#' Derive parameters of the growth-light response from light niche parameters
#'
#' The growth response is a parabolic function of light intensity:
#' \deqn{y = a * light^2 + b*light + c}.
#'
#' This function finds a, b, c such that \deqn{y = 0} when \deqn{x = MinLight}
#' or \deqn{x = MaxLight} and \deqn{y = 1} when \deqn{x = OptimumLight}
#'
#' @param MinLight numeric, the minimum light in which the species can survive
#' @param MaxLight numeric, the maximum light at which the species can suvive
#' @param OptimumLight numeric, light intensity at which growth is maximized
#'
#' @returns a vector of three numerics, parameters a, b, c of the parabolic function
#' @export
#'
get_light_resp_params <- function(MinLight, MaxLight, OptimumLight) {

  x1 <- MinLight
  y1 <- 0

  x2 <- MaxLight
  y2 <- 0

  x3 <- OptimumLight
  y3 <- 1

  # Derive coefficients from the three known points
  # see Petter et al. 2021 Appendix A2
  a <- (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) /
    ((x1 - x2) * (x1 - x3) * (x3 - x2)) # eq 24

  b <- (x1^2 * (y2 - y3) + x2^2 * (y3 - y1) + x3^2 * (y1 - y2)) /
    ((x1 - x2) * (x1 - x3) * (x2 - x3)) # eq 25

  c <- (x1^2 * (x2 * y3 - x3 * y2) +
          x1 * (x3^2 * y2 - x2^2 * y3) +
          x2 * x3 * y1 * (x2 - x3)
  ) / ((x1 - x2) * (x1 - x3) * (x2 - x3)) # eq 26

  return(c(a, b, c))
}

#' Get light parabolic response
#'
#' Transform a given light incidence into a suitability score between 0 and 1
#' using light niche parameters a, b and c.
#'
#' @param light a light incidence value
#' @param a parabolic parameter characterising the light niche
#' @param b parabolic parameter characterising the light niche
#' @param c parabolic parameter characterising the light niche
#'
#' @return a suitability score between 0 and 1
#'
#'@export
get_parabolic_resp <- function(light, a, b, c) {
  return(pmax(0, (a * light^2) + (b * light) + c))
}

