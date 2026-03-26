#' Resolve the dispersal step of the simulation
#'
#' Each year, reproduction (recruitment) and seed dispersal are resolved
#' in a singe step.
#'
#' Fecundity corresponds to the potential average number of recruits per
#' individual (\eqn{n_{RPot}}) and is the product of three elements:
#' * the linear, mass-dependent base fecundity: \eqn{a + b \times M}, where M
#' is the mass, a is `InterceptRecruitment` and b is `SlopeRecruitment`.
#' * the mass-to-reproduction allocation, `RecruitmentInvestmentRel` (a species
#' trait).
#' * a relative increase in mass reproduction allocation as the individual gets
#' closer to the maximum mass:
#' \deqn{1 + Rec_{inc} \frac{M - M_{mat}}{M_{max} - M_{mat}}} where
#' \eqn{Rec_{inc}} is species trait `RecruitmentInc`.
#' This gives the expected number of offspring in a voxel with a 1m square
#' surface area of available substrate.
#'
#' The fecundity is multiplied by the dispersal probability matrix
#' (`prob_disp_matrix`) and the available surface area matrix to obtain the
#' expected number of offspring in each voxel, which is then sampled in a
#' Poisson distribution.
#' Probabilities are set to 0 in all voxels that fall outside of the
#' light niche of the species.
#'
#' @param E a `data.frame` containing the individual epiphytes present in the
#' landscape
#' @param SpeciesPool a `data.frame` containing species-level traits
#' @param Microhabitat the microhabitat matrix, containing surface area, loss
#' and light conditions.
#' @param SurfaceBiomassScaling a strictly positive parameter scaling how much
#' surface area an individual occupies as a function of its mass:
#' \deqn{S = M^{2/3} / g_S}
#' @param InterceptRecruitment a positive number (or zero), the intercept of the
#' relation between mass and fecundity.
#' @param SlopeRecruitment a number between 0 and 1, the slope of the relation
#' between mass and fecundity.
#' @param prob_disp_matrix matrix, the output of [calc_prob_disp_matrix()].
#' Contains the 3D probability distribution of dispersing to all surrounding
#' voxels within a matrix twice as large as the `Microhabitat` matrix.
#' @param centralPoint a numeric vector of length 3 containing the central
#'  X, Y and Z coordinates of `prob_disp_matrix`
#' @param max_id integer, the highest ID among all individuals
#'
#' @returns a list containing the following elements:
#' * `nbIndsBeforeDisp` a 1D array containing the number of individuals alive
#' for each species *before* dispersal and recruitment.
#' * `E` the updated individual `data.frame`
#' * `recruitment_df` a `data.frame` summarising recruitment (expected and
#' realised number of recruits for each species) for the output.
#' * `max_id` the updated maximum individual ID
#'
#' @export
#'
resolve_repro_dispersal <- function(E,
                             Microhabitat,
                             SurfaceBiomassScaling,
                             centralPoint,
                             InterceptRecruitment,
                             SlopeRecruitment,
                             prob_disp_matrix,
                             SpeciesPool,
                             max_id) {

  dimPlot <- dim(Microhabitat)[1:3]

  # Store number of individuals at beginning of time step
  NumberOfSpecies <- nrow(SpeciesPool)
  nbIndsBeforeDisp <- array(rep(0, NumberOfSpecies))
  for (sp in seq_len(NumberOfSpecies)) {
    nbIndsBeforeDisp[sp] <- length(which(E$SpeciesID == sp & E$Status == 1))
  }

  # Habitat variables
  light_mat <- Microhabitat[,,,3]
  dim(light_mat) <- dim(Microhabitat)[1:3]
  # keep same dimensions
  avail_sa_matrix <- Microhabitat[,,,1]
  dim(avail_sa_matrix) <- dim(Microhabitat)[1:3]
  avail_sa_matrix <- get_surf_area_mat(E, avail_sa_matrix, SurfaceBiomassScaling)

  # Initialize potential recruitment dataframe
  unique_species <- unique(E$SpeciesID)
  recruitment_df <- data.frame(
    "species_index" = seq_len(NumberOfSpecies),
    "exptd_nb_recruits" = numeric(NumberOfSpecies),
    "nb_recruits" = numeric(NumberOfSpecies)
  )

  # Loop over all species
  # TODO: this loop could be parallelised
  for (i in seq_len(NumberOfSpecies)) {

    sp <- unique_species[i]
    this_species <- which(SpeciesPool$SpeciesID == sp)

    # Generate initially empty matrix to store the probabilities for recruitment
    exptd_nb_recruits_matrix <- array(
      rep(0, dimPlot[1] * dimPlot[2] * dimPlot[3]),
      dim = c(dimPlot[1], dimPlot[2], dimPlot[3])
    )

    # Matrix containing all mature individuals of one species
    mass_maturity <- SpeciesPool$MassAtMaturity[this_species]
    mature_inds <- E[E$SpeciesID == sp & E$Mass >= mass_maturity, ]

    minLight <- SpeciesPool$MinLight[i]
    maxLight <- SpeciesPool$MaxLight[i]

    # Probability matrix for each species:
    # Depending on the position of each mature individual,
    # the total probability for the species is calculated.
    #
    # The second part of the equation accounts for the actual size of the individual
    # in relation to the maximum size for which the recruitment per individual is defined
    for (j in seq_len(nrow(mature_inds))) {

      dist_to_center <- c(
        centralPoint[1] - mature_inds$X[j],
        centralPoint[2] - mature_inds$Y[j],
        centralPoint[3] - mature_inds$Z[j]
      )
      x_coords <- seq(dist_to_center[1] + 1, dist_to_center[1] + dimPlot[1])
      y_coords <- seq(dist_to_center[2] + 1, dist_to_center[2] + dimPlot[2])
      z_coords <- seq(dist_to_center[3] + 1, dist_to_center[3] + dimPlot[3])

      # Mass-dependent fecundity coefficient
      mass_coeff <- (InterceptRecruitment + SlopeRecruitment * mature_inds$Mass[j]) *
        SpeciesPool$RecruitmentInvestmentRel[this_species] # base mass-to-reproduction allocation

      # Relative mass growth since the individual has reached maturity
      # 0 = just reached maturity
      # 1 = max mass reached
      rel_growth <- (mature_inds$Mass[j] - mass_maturity) /
        (SpeciesPool$MaximumMass[this_species] - mass_maturity)
      # Reproduction allocation increases with relative growth
      incr_alloc_coeff <- 1 + (SpeciesPool$RecruitmentInc[this_species] * rel_growth)

      # Dispersal probability * fecundity = expected nb offspring in each xyz
      exptd_nb_recruits_matrix <- exptd_nb_recruits_matrix +
        prob_disp_matrix[x_coords, y_coords, z_coords, sp] *
        mass_coeff * incr_alloc_coeff
    }

    # Store potential normalized number of recruits
    # We will use this to populate sp_output_mat later
    recruitment_df$exptd_nb_recruits[i] <- sum(exptd_nb_recruits_matrix)

    # Matrix containing all voxel for which the light requirements are fulfilled
    pot_hab_matrix <- ifelse(
      light_mat >= minLight & light_mat <= maxLight,
      1, 0
    )

    # Disable unsuitable cells and scale with surface area
    exptd_nb_recruits_matrix <- exptd_nb_recruits_matrix *
      pot_hab_matrix * avail_sa_matrix # TODO: confirm fecundity scales with SA?

    # Calculate number of recruits based on final probability matrix
    nb_recruits_matrix <- array(
      stats::rpois(length(exptd_nb_recruits_matrix), exptd_nb_recruits_matrix),
      dim = dim(exptd_nb_recruits_matrix)
    )

    # Increment recruit counts
    totalNbRecruits <- sum(nb_recruits_matrix)
    recruitment_df$nb_recruits[sp] <- totalNbRecruits

    if (totalNbRecruits > 0) {

      # Add new recruits to epiphyte matrix
      recruit_coords <- extract_ind_coords(nb_recruits_matrix)
      idx_recruits <- seq(nrow(E) + 1, nrow(E) + totalNbRecruits)
      E[idx_recruits, ] <- NA
      E$X[idx_recruits] <- recruit_coords$x
      E$Y[idx_recruits] <- recruit_coords$y
      E$Z[idx_recruits] <- recruit_coords$z
      E$Mass[idx_recruits] <- 0  # Initial size
      E$Age[idx_recruits] <- 0
      E$SurfaceAreaOccupied[idx_recruits] <- 0
      E$Status[idx_recruits] <- 1  # status 1:alive
      recruits_ids <- seq(max_id + 1, max_id + totalNbRecruits)
      E$IndividualID[idx_recruits] <- recruits_ids
      E$SpeciesID[idx_recruits] <- sp

      max_id <- max_id + totalNbRecruits
    } # if any recruits

  } # species loop

  disp_items <- list(
    "nbIndsBeforeDisp" = nbIndsBeforeDisp,
    "E" = E,
    "recruitment_df" = recruitment_df,
    "max_id" = max_id
  )
  return(disp_items)
}

