check_init_dist <- function(InitDist, SpeciesPool, dimensions) {

  if (nrow(InitDist) < 1) {
    stop("InitDist is empty.")
  }

  exptd_cols <- inds_input_names()

  missing_cols <- exptd_cols[!exptd_cols %in% names(InitDist)]
  if (length(missing_cols > 0)) {
    stop(err_msg_missing_params("InitDist", missing_cols))
  }

  colnames <- names(InitDist)

  # Check for missing values
  is_missing_val <- sapply(InitDist, function(x)  {
    any(is.na(x)) || any(is.nan(x)) || any(is.null(x))
  })
  if (any(is_missing_val)) {
    wrong_params <- colnames[is_missing_val]
    stop(paste(c(
      "The following columns of InitDist contain NAs, NaN, or NULL:",
      wrong_params), rep(" ", length(wrong_params) + 1)))
  }

  if (any(InitDist$X < 1) || any(InitDist$X > dimensions[1]) ||
          any(InitDist$Y < 1) || any(InitDist$Y > dimensions[2]) ||
          any(InitDist$Z < 1) || any(InitDist$Z > dimensions[3])) {
    stop("InitDist contains coordinates that fall outside of the microhabitat matrix")
  }

  if (any(!InitDist$SpeciesID %in% SpeciesPool$SpeciesID)) {
    stop("One or more species in InitDist are absent from SpeciesPool.")
  }

  if (any(InitDist < 0)) {
    stop("InitDist contains one or more individuals with negative mass.")
  }

  for (i in 1:nrow(InitDist)) {
    mass <- InitDist$Mass[i]
    sp <- InitDist$SpeciesID[i]
    max_mass <- SpeciesPool$MaximumMass[which(SpeciesPool$SpeciesID == sp)]
    if (mass > max_mass) {
      stop(paste(
        "InitDist contains an individual larger than is allowed by its species (row", i, ")"))
    }
  }

  if (any(duplicated(InitDist$IndividualID))) {
    stop("Multiple individuals in InitDist share the same ID.")
  }

  if (any(InitDist$Age < 0)) {
    stop("Some individuals in InitDist have a negative age.")
  }

  if (any(InitDist$SurfaceAreaOccupied < 0)) {
    stop("Some individuals in InitDist have a negative surface area requirement.")
  }

  if (any(!InitDist$Status %in% c(1, 2))) {
    stop("All individuals in InitDist must have Status 1 or 2.")
  }

}
