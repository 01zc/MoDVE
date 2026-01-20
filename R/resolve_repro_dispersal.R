#' Resolve the dispersal step of the simulation
#'
#' @param E epiphyte data frame
#' @param Microhabitat microhabitat matrix
#' @param SurfaceBiomassScaling numeric
#' @param centralPoint numeric vector of length 3 with coordinates of the central point
#' @param InterceptRecruitment numeric parameter
#' @param SlopeRecruitment numeric parameter
#' @param prob_disp_matrix matrix, the output of `calc_prob_dist_matrix()`
#' @param SpeciesPool data frame containing the species parameters
#' @param max_id integer, the highest ID among all individuals
#'
#' @returns a list
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

  # Deduce surface area
  avail_sa_matrix <- get_surf_area_mat(E, Microhabitat[,,,1], sim_params$SurfaceBiomassScaling)

  # Initialize potential recruitment dataframe
  unique_species <- unique(E$SpeciesID)
  recruitment_df <- data.frame(
    species_index = seq_len(NumberOfSpecies),
    exptd_nb_recruits = numeric(NumberOfSpecies),
    nb_recruits = numeric(NumberOfSpecies)
  )

  # Loop over all species
  # TODO: this loop could be parallelised
  for (i in seq_len(NumberOfSpecies)) {

    sp <- unique_species[i]

    # Generate initially empty matrix to store the probabilities for recruitment
    exptd_nb_recruits_matrix <- array(
      rep(0, dimPlot[1] * dimPlot[2] * dimPlot[3]),
      dim = c(dimPlot[1], dimPlot[2], dimPlot[3])
    )

    # Matrix containing all mature individuals of one species
    mature_inds <- E[E$SpeciesID == sp & E$Mass >= E$MassAtMaturity, ]

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
      mass_coeff <- (InterceptRecruitment + SlopeRecruitment) * mature_inds$Mass[j] *
        mature_inds$RecruitmentInvestmentRel[j] # base mass-to-reproduction allocation

      # Relative mass growth since the individual has reached maturity
      # 0 = just reached maturity
      # 1 = max mass reached
      rel_growth <- (mature_inds$Mass[j] - mature_inds$MassAtMaturity[j]) /
        (mature_inds$MaximumMass[j] - mature_inds$MassAtMaturity[j])
      # Reproduction allocation increases with relative growth
      incr_alloc_coeff <- 1 + (mature_inds$RecruitmentInc[j] * rel_growth)

      # Dispersal probability * fecundity = expected nb offspring in each xyz
      exptd_nb_recruits_matrix <- exptd_nb_recruits_matrix +
        prob_disp_matrix[x_coords, y_coords, z_coords, sp] * mass_coeff * incr_alloc_coeff
    }

    # Store potential normalized number of recruits
    # We will use this to populate sp_output_mat later
    recruitment_df$exptd_nb_recruits[i] <- sum(exptd_nb_recruits_matrix)

    # Matrix containing all voxel for which the light requirements are fulfilled
    pot_hab_matrix <- ifelse(
      Microhabitat[, , , 3] >= minLight &
        Microhabitat[, , , 3] <= maxLight,
      1, 0
    )

    # Disable unsuitable cells and scale with surface area
    exptd_nb_recruits_matrix <- exptd_nb_recruits_matrix *
      pot_hab_matrix * avail_sa_matrix # TODO: confirm fecundity scales with SA?

    # Calculate number of recruits based on final probability matrix
    nb_recruits_matrix <- array(
      rpois(length(exptd_nb_recruits_matrix), exptd_nb_recruits_matrix),
      dim = dim(exptd_nb_recruits_matrix)
    )

    # Increment recruit counts
    totalNbRecruits <- sum(nb_recruits_matrix)
    recruits_tbl$nb_recruits[sp] <- totalNbRecruits

    if (totalNbRecruits > 0) {

      # Add new recruits to epiphyte matrix
      recruit_coords <- extract_ind_coords(nb_recruits_matrix)
      idx_recruits <- seq(nrow(E) + 1, nrow(E) + totalNbRecruits)
      E[idx_recruits, names(SpeciesPool)] <- SpeciesPool[sp, ]
      E$X[idx_recruits] <- recruit_coords$x
      E$Y[idx_recruits] <- recruit_coords$y
      E$Z[idx_recruits] <- recruit_coords$z
      E$Mass[idx_recruits] <- 0  # Initial size
      E$Status[idx_recruits] <- 1  # status 1:alive
      recruits_ids <- seq(max_id + 1, max_id + length(totalNbRecruits))
      E$IndividualID[idx_recruits] <- recruits_ids

      E[is.na(E)] <- 0  # convert all NA to 0 so that the R script matches the Matlab
      # TODO: this is likely to cause bugs

      max_id <- max_id + length(x_recruits)
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

