#' Draw the trait value and location of initial individuals
#'
#' Draw the initial mass, age, occupied surface area of initial individuals
#' and distirbute them along eligible voxels in the habitat matrix.
#'
#' @param distr_params a list containing at least the following parameters:
#' * `IndividualsPerSpecies` (scalar) integer, the number of individuals of per
#' species.
#' * `ScalingPerHa` logical, defines whether the NumberSpecies are total numbers
#' irrespective of the model area (`TRUE`), or if the NumberSpecies or given per
#' hectar and are scaled to the model area (`FALSE`)
#' * `PercentageMaturePerSpecies` a number between 0 and 100 indicating the
#' percentage of initial individuals that should start as mature individuals
#' able to reproduce
#' * `SurfaceBiomassScaling` a strictly positive parameter (\eqn{g_S}) scaling
#' the surface area requirement as a function of the mass of the individual:
#' \deqn{S = M^{2/3} / g_S}
#' @param species_df a data frame containing the species traits, as created by
#' `draw_species_traits`
#' @param microhab_mat a matrix containing the surface area and light conditions
#' in each voxel, as created by [create_microhabitat_mat()]
#' @param path_to_output a string, where (folder and name) should the output
#' be saved? Must be `.csv`. If `NULL` (the default), no output is saved and the
#' output is returned as a data frame instead.
#' @param largest_inds_first logical, should larger individuals be allocated
#' first (`TRUE`) to represent competitive advantage, or randomly (`FALSE`)
#' @param most_surf_area_first logical, should individuals first be allocated to
#' the largest voxel available (`TRUE`) or a random one (`FALSE`)?
#'
#' @return if `path_to_output = NULL`, returns a data frame with one row per
#' initial individual and the following columns:
#' * `X` the x-coordinate of the individual
#' * `Y` the y-coordinate of the individual
#' * `Z` the z-coordinate of the individual
#' * `Mass` the mass of the individual
#' * `Status`, the status of the individual either 1 (alive) or 2 (dead, due to
#' a lack of available surface area to allocate this individual)
#' * `IndividualID` a unique identifier for this individual
#' * `SurfaceAreaOccupied` the amount of surface area that this individual
#' requires and uses
#' * `Age` age of the individual in years
#' * `SpeciesID` which species this individual belongs to.
#'
#' @export
#'
draw_initial_individuals <- function(distr_params, species_df,
                                        microhab_mat, path_to_output = NULL,
                                        largest_inds_first = FALSE,
                                        most_surf_area_first = FALSE
                                        ) {

  check_distr_params(distr_params)
  check_species_df(species_df)
  check_microhabitat(microhab_mat)

  nb_inds_per_sp <- distr_params$IndividualsPerSpecies

  dimPlot <- dim(microhab_mat)[1:3]
  if (distr_params$ScalingPerHa) {
    # TODO: unclear from doc what this parameter actually does
    # doc mentions it's a scaling nb species yet this suggest nb individuals
    nb_inds_per_sp <- nb_inds_per_sp * dimPlot[1] * dimPlot[2] / 10000
  }

  NumberSpecies <- length(species_df$SpeciesID)
  TotalIndividuals <- NumberSpecies * distr_params$IndividualsPerSpecies
  nb_mature_inds <- round(
    nb_inds_per_sp * distr_params$PercentageMaturePerSpecies / 100
  )
  nb_immature_inds <- nb_inds_per_sp - nb_mature_inds

  # Initialize individual matrix
  col_names <- inds_input_names()
  init_ind_mat <- array(
    rep(0, TotalIndividuals * length(col_names)),
    dim = c(TotalIndividuals, length(col_names))
  )
  # Pre-assign matrix columns to variables
  col_x <- match("X", col_names)
  col_y <- match("Y", col_names)
  col_z <- match("Z", col_names)
  col_mass <- match("Mass", col_names)
  col_status <- match("Status", col_names)
  col_id <- match("IndividualID", col_names)
  col_sa <- match("SurfaceAreaOccupied", col_names)
  col_age <- match("Age", col_names)
  col_sp <- match("SpeciesID", col_names)

  # Individual ID
  init_ind_mat[, col_id] <- seq_len(TotalIndividuals)

  # Species depdt qualities
  for (sp in seq_len(NumberSpecies)) {

    first_row <- (sp - 1) * nb_inds_per_sp + 1
    last_row <- nb_inds_per_sp * sp

    # Species ID
    spID <- species_df$SpeciesID[sp]
    init_ind_mat[first_row:last_row, col_sp] <- spID

    # Mass
    massMaturity <- species_df$MassAtMaturity[sp]
    maxMass <- species_df$MaximumMass[sp]
    row_matures <- first_row + nb_mature_inds - 1
    init_ind_mat[int_seq(first_row, row_matures), col_mass] <- stats::runif(
      nb_mature_inds, min = massMaturity, max = maxMass
    )
    init_ind_mat[int_seq((row_matures + 1), last_row), col_mass] <- stats::runif(
      nb_immature_inds, min = 0, max = massMaturity
    )

    # Age
    growthRate <- species_df$GrowthRate[sp]
    init_ind_mat[int_seq(first_row, last_row), col_age] <- round(
      -log(1 - (init_ind_mat[int_seq(first_row, last_row), col_mass] / maxMass)) /
        growthRate
    )
  } # species loop

  # Surface area occupied
  init_ind_mat[, col_sa] <- (init_ind_mat[, col_mass]^(2 / 3)) /
    distr_params$SurfaceBiomassScaling

  # Allocate voxels (x, y, z, status)

  # Available surface decreases with each allocated individual
  AvailableSurfaceArea <- array(
    microhab_mat[,,,1], dim = dim(microhab_mat)[1:3] # keep same dimensions
    )

  # Schedule individual allocation priority
  if (largest_inds_first) { # largest individuals win competition
    ind_queue <- order(init_ind_mat[, col_mass], decreasing = TRUE)
  } else { # random
    ind_queue <- sample(seq_len(TotalIndividuals), TotalIndividuals, replace=FALSE)
  }

  for (ind in ind_queue) {

    sp <- init_ind_mat[ind, col_sp]

    # Find all suitable voxels for this individual
    AreaNeededInd <- init_ind_mat[ind, col_sa]
    hasEnoughSurface <- AvailableSurfaceArea > AreaNeededInd
    hasEnoughLight <- array(
      microhab_mat[, , , 3] >= species_df$MinLight[sp],
      dim = dim(microhab_mat)[1:3] # keep same dimensions
    )
    hasEnoughShade <- array(
      microhab_mat[, , , 3] <= species_df$MaxLight[sp],
      dim = dim(microhab_mat)[1:3] # keep same dimensions
    )
    SuitableVoxels <- which(hasEnoughLight & hasEnoughShade & hasEnoughSurface)

    if (length(SuitableVoxels) > 0) {

      # Select one voxel
      if (most_surf_area_first) { # voxel with the most available surface area
        whichVoxel <- which(AvailableSurfaceArea[SuitableVoxels] ==
                              max(AvailableSurfaceArea[SuitableVoxels]))[1]
      } else { # random voxel
        whichVoxel <- sample.int(length(SuitableVoxels), size=1)
      }
      ids <- arrayInd(SuitableVoxels[whichVoxel], dim(AvailableSurfaceArea))
      x <- ids[, 1]
      y <- ids[, 2]
      z <- ids[, 3]

      # Update available Surface Area
      AvailableSurfaceArea[x, y, z] <- AvailableSurfaceArea[x, y, z] - AreaNeededInd

      # Update Initial Epiphyte Matrix
      init_ind_mat[ind, col_x] <- x
      init_ind_mat[ind, col_y] <- y
      init_ind_mat[ind, col_z] <- z

      # Set status of individual: status=1 => alive
      init_ind_mat[ind, col_status] <- 1

    } else { # no suitable voxels

      # Set status of individual: status=2 =>> dead
      init_ind_mat[ind, col_status] <- 2

      # Set coordinates to 1 (might cause problems in later model if not)
      init_ind_mat[ind, col_x] <- 1
      init_ind_mat[ind, col_y] <- 1
      init_ind_mat[ind, col_z] <- 1
    }
  } # loop individuals

  # Convert to df and save
  init_ind_mat_df <- as.data.frame(init_ind_mat)
  names(init_ind_mat_df) <- col_names

  # Save Initial Epiphyte Matrix
  if (!is.null(path_to_output)) {
    utils::write.csv(init_ind_mat_df, path_to_output, row.names = FALSE)
  } else {
    return(init_ind_mat_df)
  }
}

check_distr_params <- function(distr_params) {

  exptd_params <- c(
    "IndividualsPerSpecies",
    "ScalingPerHa",
    "PercentageMaturePerSpecies",
    "SurfaceBiomassScaling"
  )

  missing_params <- exptd_params[!exptd_params %in% names(distr_params)]
  if (length(missing_params > 0)) {
    stop(err_msg_missing_params("distr_params", missing_params))
  }

  # No NAs, NULL, or NaN!
  is_missing_val <- sapply(distr_params, function(x)  {
    any(is.na(x)) || any(is.nan(x)) || any(is.null(x))
  })
  if (any(is_missing_val)) {
    wrong_params <- exptd_params[is_missing_val]
    stop(paste(c(
      "The following elements of distr_params are NA, NaN, or NULL:",
      wrong_params), rep(" ", length(wrong_params) + 1)))
  }

  if (!is.logical(distr_params$ScalingPerHa)) {
    stop("ScalingPerHa must be TRUE or FALSE.")
  }

  numeric_params <- exptd_params[exptd_params != "ScalingPerHa"]
  is_numeric <- sapply(distr_params[numeric_params], function(x) all(is.numeric(x)))
  if (any(!is_numeric)) {
    wrong_params <- numeric_params[!is_numeric]
    stop(paste(c(
      "The following elements of distr_params must be numeric:",
      wrong_params), rep(" ", length(wrong_params) + 1)))
  }

  if (distr_params$PercentageMaturePerSpecies < 0 ||
      distr_params$PercentageMaturePerSpecies > 100) {
    stop("PercentageMaturePerSpecies must be a number between 0 and 100.")
  }

  if (distr_params$IndividualsPerSpecies < 0) {
    stop("IndividualsPerSpecies must be a positive number.")
  }

  if ((distr_params$IndividualsPerSpecies %% 1) != 0) {
    stop("IndividualsPerSpecies must be an integer.")
  }

  if (distr_params$SurfaceBiomassScaling <= 0) {
    stop("SurfaceBiomassScaling must be a strictly positive number.")
  }

}
