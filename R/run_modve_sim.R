#' Run the MoDVE simulation
#'
#' Main function of the package, runs the epiphyte community simulation.
#'
#' @param sim_params a list of parameters containing at least the following
#' elements:
#'  * `InitialTimeSteps`: index for the first year
#'  * `TimeSteps`: number of time steps (years) to run the simulation for
#'  * `StopCriterionHa`: limit density of individuals per ha. If the number of
#' individuals exceeds `StopCriterionHa / 10000 * dimX * dimY`, the simulation
#' exits.
#'  * `hasDynamicMicrohabitat`: `TRUE/FALSE`, does the microhabitat matrix change
#' with each time steps? If `TRUE`, `Microhabitat` must be a vector of paths
#' to each of the microhabitat matrices (one per year).
#' * `massDepCompetition`: `TRUE` = larger individuals get priority in
#' voxel attribution, otherwise (`FALSE`) individuals are distributed randomly.
#' * `use_mass_dep_mortality`: if `FALSE` individuals die randomly
#' according to `MortRateRandom`, if `TRUE` mortality is mass-dependent, using
#' `MortRateMass * (mass^MortRateMassScaling)`.
#' * `MortRateRandom`: numeric between 0 and 1, the probability of an
#' individual dying if `use_mass_dep_mortality == FALSE`
#' * `MortRateMass`: numeric, if `use_mass_dep_mortality == TRUE` the
#' coefficient for the effect of mass on the probability of death
#' * `MortRateMassScaling`: numeric between `-Inf` and `0`, if
#' `use_mass_dep_mortality == TRUE` the exponent for the effect of mass on the
#' probability of death. Must be negative or zero.
#' * `SurfaceBiomassScaling`: a strictly positive parameter (\deqn{g_S}) scaling how much
#' surface area an individual occupies as a function of its mass:
#' \deqn{S = M^{2/3} / g_S}
#' * `InterceptRecruitment`: a positive number (or zero), the intercept of the
#' relation between mass and fecundity.
#' * `SlopeRecruitment`: a number between 0 and 1, the slope of the relation
#' between mass and fecundity
#'
#' @param SpeciesPool a `data.frame` containing the species traits, as generated
#' e.g. with [draw_species_traits()]. Contains one row per species and the
#' following columns:
#' * `MaximumMass` positive numeric, the maximum mass an individual can reach.
#'  * `MassAtMaturity` numeric between 0 and and `MaximumMass`,
#'  fraction of `MaximumMass` above which at individual can reproduce.
#'  * `GrowthRate` numeric between 0 an 1, the fraction of remaining growth an
#'  individual gains in a single year (i.e, \eqn{\Delta m = K *
#'  (M_{max} - M)}), under optimal light conditions.
#'  * `DispersalKernel` positive numeric, the dispersal kernel.
#'  * `DispersalKernelAsymmetry` numeric between 0 and 1, the dispersal asymmetry.
#'  \eqn{D_{k_A} = 0.5} corresponds to symmetric dispersal; with \eqn{D_{k_A} = 1} individuals
#'  disperse stricly below themselves; with \eqn{D_{k_A} = 0} individuals never disperse
#'  only above themselves or at their height.
#'  * `RecruitmentInvestmentRel` numeric between 0 and 1, a coefficient scaling
#'  the mass-dependent fecundity coefficient (see [resolve_repro_dispersal()]),
#'  representing the fraction of available biomass invested in fecundity
#'  * `RecruitmentInc` numeric between 0 and 1, scaling how fecundity increases
#'  with the growth stage of the epiphyte. This factor ranges from `1` when
#'  `Mass` = `MassAtMaturity`, and `2 * RecruitmentInc` when
#'  `Mass` = `MaximumMass`.
#'  * `MinLight` positive numeric, minimum light conditions under which this species can survive
#'  * `MaxLight` positive numeric, maximum light conditions under which this species can survive
#'  * `OptimumLight` positive numeric, optimum light conditions under which individuals of this species
#'  grow and reproduce at the maximum rate. It is calculated as the average of `MinLight` and `MaxLight.`
#'  * `LightResponseA` first term of the parabolic light-growth response function.
#'  * `LightResponseB` second term of the parabolic light-growth response function.
#'  * `LightResponseC` third term of the parabolic light-growth response function.
#'
#' @param Microhabitat a 4D matrix where the first three dimensions
#' corresponding to a 3D habitat space, and the last one containing values of:
#' 1. the microhabitat for the available surface area,
#' 2. % of surface area lost in the previous year (if dynamic) and
#' 3. light intensity in each voxel;
#' or a path to a csv file containing such a matrix.
#' If `hasDynamicMicrohabitat` is `TRUE`, `Microhabitat` must be a vector of
#' paths to such matrices`,` with length `Timesteps`.
#'
#' @param InitDist a `data.frame` containing the distribution and initial
#' attributes of individuals at the beginning of the simulation, as generated
#' e.g. with [draw_initial_individuals()]. Contains one row per initial individual,
#' and the following columns:
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
#' * `SpeciesID` which species this individual belongs
#' to.
#'
#' @param path_to_ind_output where to save the individual-level output.
#' Must end with `.csv`.
#' @param path_to_sp_output where to save the species-level output.
#' Must end with `.csv`.
#' @param path_to_comm_output where to save the community-level output.
#' Must end with `.csv`.
#'
#'@export
run_modve_sim <- function(sim_params,
                          SpeciesPool,
                          Microhabitat,
                          InitDist,
                          path_to_ind_output,
                          path_to_sp_output,
                          path_to_comm_output) {

  # Check parameters and read input if necessary
  check_sim_params(sim_params)
  timeSteps <- sim_params$timeSteps

  if (!is.data.frame(SpeciesPool)) {
    # Then it must be a path to the input
    if (!grepl("*.csv$", SpeciesPool)) {
      stop("SpeciesPool should be a data frame or a valid path to a .csv file.")
    }
    if (!file.exists(SpeciesPool)) {
      stop(paste0(SpeciesPool, " doesn't exist.\n"))
    } # error
    else SpeciesPool <- utils::read.csv(SpeciesPool, sep = ",", header = TRUE)
  }
  check_species_df(SpeciesPool)
  NumberOfSpecies <- nrow(SpeciesPool)  # number of species per 25x25m plot

  if (!is.array(Microhabitat)) {
    # Then it must be a path or vector of paths
    for (i in seq_along(Microhabitat)) {
      if (!grepl("*.rds$", Microhabitat[i])) {
        stop("Microhabitat should be an array or a valid path to a .rds file.")
      }
      if (!file.exists(Microhabitat[i])) {
        stop(paste0(Microhabitat[i], " doesn't exist.\n"))
      }
    }

    if (sim_params$hasDynamicMicrohabitat) {
      if (!length(Microhabitat) == timeSteps) {
        stop("For dynamic habitats, Microhabitat should have one element for each time step.\n")
      } else {
        # Stash paths for later time steps
        microhab_files <- Microhabitat
      }
    }
    # If all checks ok, read the first one
    Microhabitat <- readRDS(Microhabitat[1])
    check_microhabitat(Microhabitat)
  }

  dims <- dim(Microhabitat)
  dimX <- dims[1]
  dimY <- dims[2]
  dimZ <- dims[3]

  if (!is.data.frame(InitDist)) {
    # Then it must be a path to the input
    if (!grepl("*.csv$", InitDist)) {
      stop("InitDist should be a data frame or a valid path to a .csv file.")
    }
    if (!file.exists(InitDist)) {
      stop(paste0(InitDist, " doesn't exist.\n"))
    }
    else InitDist <- utils::read.csv(InitDist, sep = ",", header = TRUE)
  }
  check_init_dist(InitDist, SpeciesPool, dims)
  E <- InitDist
  # Add columns to E for additional info
  E[, c("TotalSurfaceInVoxel", "LightInVoxel", "SurfaceLossInVoxel")] <- 0
  max_id <- nrow(E)  # to trace individual IDs

  # Prepare output
  dir_output <- dirname(path_to_ind_output)
  if (!dir.exists(dir_output)) {
    stop(paste0("Output directory ", dir_output, " does not exist."))
  }
  if (!grepl("*.csv$", path_to_ind_output)) {
    stop("path_to_ind_output must be a csv file")
  }

  dir_output <- dirname(path_to_sp_output)
  if (!dir.exists(dir_output)) {
    stop(paste0("Output directory ", dir_output, " does not exist."))
  }
  if (!grepl("*.csv$", path_to_sp_output)) {
    stop("path_to_sp_output must be a csv file")
  }

  sp_output_headers <- species_output_names()
  {
    # Column indices
    col_sp_t <- 1; col_sp_id <- 2; col_nb_inds_begin <- 3; col_nb_inds_end <- 4;
    col_nb_mature_inds <- 5; col_nb_rec <- 6; col_nb_rec <- 7;
    col_nb_dead_branch <- 8; col_nb_dead_light <- 9; col_nb_dead_comp <- 10;
    col_nb_dead_base <- 11; col_growth_rate <- 12; col_growth_log <- 13
    col_birth <- 14; col_death <- 15; col_size <- 16; col_avg_age <- 17
    col_min_light <- 18; col_max_light <- 19; col_mean_light <- 20;
    col_min_height <- 21; col_max_height <- 22; col_mean_height <- 23
    nb_cols_sp_output <- col_mean_height
  }
  # Initialize Matrix where species parameters are saved
  sp_output <- array(
    rep(0, (timeSteps * NumberOfSpecies) * nb_cols_sp_output),
    dim = c(timeSteps * NumberOfSpecies, nb_cols_sp_output)
  )

  dir_output <- dirname(path_to_comm_output)
  if (!dir.exists(dir_output)) {
    stop(paste0("Directory ", dir_output, " does not exist."))
  }
  if (!grepl("*.csv$", path_to_comm_output)) {
    stop("path_to_comm_output must be a csv file")
  }
  comm_output_headers <- comm_output_names()
  comm_output <- data.frame(matrix(
    0.0, nrow = timeSteps, ncol = length(comm_output_headers)
  ))
  colnames(comm_output) <- comm_output_headers

  # Stop criterion: if stop density is exceeded, the simulation ends
  # to cap memory usage
  StopNbInds <- sim_params$StopCriterionHa / 10000 * dimX * dimY

  # Calculate the probability to disperse in surrounding voxels
  # Generate probabilities for a matrix twice as large as microhabitat
  expanded_dims <- dims[1:3] * 2 + 1
  expanded_mat_central_point <- floor(expanded_dims / 2) + 1
  prob_disp_matrix <- calc_prob_disp_matrix(
    expanded_mat_central_point,
    expanded_dims,
    SpeciesPool
  )

  # Year loop
  for (t in seq_len(timeSteps)) {

    year_nb <- sim_params$InitialTimeStep + t - 1 # actual time step

    # Check if the stop criterion is met
    nbIndsAlive <- length(which(E$Status == 1))
    if (nbIndsAlive > StopNbInds) {
      writeLines(paste0(
        "Time ", year_nb, ": population has exceeded max threshold of ", StopNbInds,
        " individuals. Ending simulation."
        ))
      break
    }

    # Update microhabitat if applicable
    if (sim_params$hasDynamicMicrohabitat && t > 1) {
      Microhabitat <- readRDS(microhab_files[t])
      check_microhabitat(Microhabitat)
      if (!all.equal(dim(Microhabitat), dims)) {
        stop(
          paste("Invalid microhabitat matrix at time", year_nb,
                ": number of dimensions must be the same as the first matrix")
          )
      }
    }

    # Update how many species are alive at beginning of the year
    InitialNumberSpecies <- length(unique(E$SpeciesID[E$Status == 1]))
    nbIndsBeforeDispTotal <- length(which(E$Status == 1))

    # Dispersal
    disp_items <- resolve_repro_dispersal(
      E, Microhabitat, sim_params$SurfaceBiomassScaling,
      expanded_mat_central_point, sim_params$InterceptRecruitment,
      sim_params$SlopeRecruitment, prob_disp_matrix,  SpeciesPool, max_id
    )

    # Unwrap dispersal output
    nbIndsBeforeDisp <- disp_items$nbIndsBeforeDisp
    recruitment_df <- disp_items$recruitment_df
    E <- disp_items$E
    max_id <- disp_items$max_id

    # Store potential normalized number of recruits in sp_output
    for (i in seq_len(nrow(recruitment_df))) {
      sp_idx <- recruitment_df$index[i]
      row_idx <- (sp_idx - 1) * timeSteps + t
      sp_output[row_idx, col_nb_rec] <-
        recruitment_df$nb_potential_recruits[i]
    }
    NumberRecruits <- length(which(E$Status == 1)) - nbIndsBeforeDispTotal
    nbRecruitsPerSpecies <- disp_items$recruitment_df$nb_recruits

    # Growth
    E <- resolve_growth(E, SpeciesPool, Microhabitat, sim_params$SurfaceBiomassScaling)

    # Mortality (except from competition)
    E <- resolve_mortality(E, SpeciesPool, Microhabitat, sim_params$use_mass_dep_mortality,
                           sim_params$MortRateRandom, sim_params$MortRateMass,
                           sim_params$MortRateMassScaling)

    # Mortality due to competition for space
    E <- resolve_competition(E, Microhabitat, sim_params$massDepCompetition)

    # Age increment
    E$Age <- E$Age + 1

    # Species-level output
    for (sp in seq_len(NumberOfSpecies)) {

      row_nb <- (sp - 1) * timeSteps + t
      is_sp <- E$SpeciesID == sp

      nb_alive <- sum(E$Status == 1 & is_sp, na.rm = TRUE)
      nb_dead_comp <- sum(E$Status == 2 & is_sp, na.rm = TRUE)
      nb_dead_branch <- sum(E$Status == 3 & is_sp, na.rm = TRUE)
      nb_dead_light <- sum(E$Status == 4 & is_sp, na.rm = TRUE)
      nb_dead_base <- sum(E$Status == 5 & is_sp, na.rm = TRUE)
      nb_inds_begin <- nbIndsBeforeDisp[sp]

      sp_output[row_nb, col_sp_t] <- year_nb
      sp_output[row_nb, col_sp_id] <- sp
      sp_output[row_nb, col_nb_inds_begin] <- nb_inds_begin
      sp_output[row_nb, col_nb_inds_end] <- nb_alive
      sp_output[row_nb, col_nb_mature_inds] <- sum(
        E$Status == 1 & is_sp & E$Mass >= SpeciesPool$MassAtMaturity[sp],
        na.rm = TRUE
      )
      sp_output[row_nb, col_nb_rec] <- nbRecruitsPerSpecies[sp]
      sp_output[row_nb, col_nb_dead_branch] <- nb_dead_branch
      sp_output[row_nb, col_nb_dead_light] <- nb_dead_light
      sp_output[row_nb, col_nb_dead_comp] <- nb_dead_comp
      sp_output[row_nb, col_nb_dead_base] <- nb_dead_base

      if (nb_alive > 0 &&  nb_inds_begin > 0) {
        sp_output[row_nb, col_growth_rate] <- sp_output[row_nb, col_nb_inds_end] /
          sp_output[row_nb, col_nb_inds_begin]
        sp_output[row_nb, col_growth_log] <- log(sp_output[row_nb, col_growth_rate])
        sp_output[row_nb, col_birth] <- nbRecruitsPerSpecies[sp] /
          nb_inds_begin
        sp_output[row_nb, col_death] <- (
          nb_dead_branch + nb_dead_light + nb_dead_comp + nb_dead_base
        ) / nb_inds_begin
        sp_output[row_nb, col_size] <- mean(E$Mass[is_sp])
        sp_output[row_nb, col_avg_age] <- mean(E$Age[is_sp])
        sp_output[row_nb, col_min_light] <- min(E$LightInVoxel[is_sp])
        sp_output[row_nb, col_max_light] <- max(E$LightInVoxel[is_sp])
        sp_output[row_nb, col_mean_light] <- mean(E$LightInVoxel[is_sp])
        sp_output[row_nb, col_min_height] <- min(E$Z[is_sp])
        sp_output[row_nb, col_max_height] <- max(E$Z[is_sp])
        sp_output[row_nb, col_mean_height] <- mean(E$Z[is_sp])
      } else {
        sp_output[row_nb, col_growth_rate:col_mean_height] <- NA
      }

    } # species loop

    # Community-level output
    MortalityCompetition <- length(which(E$Status == 2))
    MortalityBranchFall <- length(which(E$Status == 3))
    MortalityLight <- length(which(E$Status == 4))
    MortalityNatural <- length(which(E$Status == 5))
    comm_output$timeStep[t] <- year_nb
    comm_output$NumberSpeciesBeginning[t] <- InitialNumberSpecies
    comm_output$NumberSpeciesEnd[t] <- length(unique(E$SpeciesID[E$Status == 1]))
    comm_output$NumberIndividualsBeginning[t] <- nbIndsBeforeDispTotal
    comm_output$NumberIndividualsEnd[t] <- length(which(E$Status == 1))
    comm_output$nb_recruits_matrix[t] <- NumberRecruits
    comm_output$MortalityBranchFall[t] <- MortalityBranchFall
    comm_output$MortalityLight[t] <- MortalityLight
    comm_output$MortalityCompetition[t] <- MortalityCompetition
    comm_output$MortalityNatural[t] <- MortalityNatural
    comm_output$BranchSurfaceIndex[t] <- sum(Microhabitat[, , , 1]) /
      (dimX[1] * dimY[2])
    comm_output$EpiphyteFilling[t] <- sum(
      mass_to_surf_area(E$Mass, sim_params$SurfaceBiomassScaling)
      ) / sum(Microhabitat[, , , 1])

    # Command window information
    msg <- "--------------------------------------------"
    msg <- paste_wrap(msg, "Time step: ", year_nb)
    msg <- paste_wrap(msg, "Number of individuals: ", comm_output$NumberIndividualsEnd[t])
    msg <- paste_wrap(msg, "Number of species: ", comm_output$NumberSpeciesEnd[t])
    msg <- paste_wrap(msg, "Number of recruits: ", NumberRecruits)
    msg <- paste_wrap(msg, "MortalityBranchFall: ", MortalityBranchFall)
    msg <- paste_wrap(msg, "MortalityLight: ", MortalityLight)
    msg <- paste_wrap(msg, "MortalityCompetition: ", MortalityCompetition)
    msg <- paste_wrap(msg, "MortalityNatural: ", MortalityNatural)
    msg <- paste_wrap(msg, "Time: ", format(Sys.time(), "%H:%M:%OS3"))
    writeLines(msg)

    # Save Epiphyte matrix for every time step
    ind_output_file <- sub("*.csv$", paste0("_", year_nb, ".csv"), path_to_ind_output)
    utils::write.csv(
      E[, inds_output_names()],
      ind_output_file,
      row.names = FALSE
      )

    # Save sp_output for every time step
    sp_output_df <- as.data.frame(sp_output)
    names(sp_output_df) <- sp_output_headers
    utils::write.csv(
      sp_output_df,
      path_to_sp_output,
      row.names = FALSE
      )

    # Save comm_output for every time step (overwrite old one)
    utils::write.csv(
      comm_output,
      path_to_comm_output,
      append = FALSE,
      row.names = FALSE
      )

    # Remove dead individuals from Epimatrix
    E <- E[E$Status <= 1, ]

  } # time loop

}
