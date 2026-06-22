check_suitability <- function(SuitabilityMat, dims, nb_species) {
  dims_suit <- dim(SuitabilityMat)
  if (length(dims_suit) != 4) {
    stop("SuitabilityMat must have strictly 4 dimensions.")
  }
  if (any(dims_suit[1:3] != dims)) {
    stop("The x, y, and/or z dimensions of the Microhabitat and Suitability matrices differ.")
  }
  if (dims_suit[4] != nb_species) {
    stop("The fourth dimension of Suitability does not match the number of species.")
  }

  if (any(is.na(SuitabilityMat))) {
    stop("Suitability score matrix contains one or more NAs.")
  }
  if (any(is.nan(SuitabilityMat))) {
    stop("Suitability score matrix contains one or more NaN values.")
  }
  if (any(SuitabilityMat < 0)) {
    stop("Suitability scores must be strictly between 0 and 1.")
  }
  if (any(SuitabilityMat > 1)) {
    stop("Suitability scores must be strictly between 0 and 1.")
  }
}
