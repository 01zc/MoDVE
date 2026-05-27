calc_suitability <- function(val, niche_min, niche_max, niche_opt) {

  # Compute suitability only for valid entries
  num <- (niche_max - val) / (niche_max - niche_opt)
  denom <- (val - niche_min) / (niche_opt - niche_min)
  expo  <- (niche_opt - niche_min) / (niche_max - niche_opt)

  suitability <- num * denom ^ expo

  return(suitability)
}
