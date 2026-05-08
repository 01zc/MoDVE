calc_suitability <- function (MinEnvVar, MaxEnvVar, OptEnvVar, EnvVar) {

  # Pre-compute denominators
  MaxOptDiff <- MaxEnvVar - OptEnvVar
  OptMinDiff <- OptEnvVar - MinEnvVar

  # Compute suitability only for valid entries
  num <- (MaxEnvVar - EnvVar) / MaxOptDiff
  denom <- (EnvVar - MinEnvVar) / OptMinDiff
  expo  <- OptMinDiff / MaxOptDiff

  suitability <- num * denom^expo

  return(suitability)  # shape: e.g. [50, 50, 60, 100, 2]
}
