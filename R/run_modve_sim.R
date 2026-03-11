#' Run the MoDVE simulation
#'
#' Main function of the package, runs the epiphyte lifecycle simulation.
#'
#' @param sim_params description
#' @param SpeciesPool descriptioflm
#' @param Microhabitat desc
#' @param InitDist description
#' @param path_to_ind_output description
#' @param path_to_sp_output description
#' @param path_to_comm_output description
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
  exptd_params <- c(
    "InitialTimeStep", "timeSteps", "StopCriterionHa","MicrohabitatType",
    "Imax", "CompetitionMethod", "MortalityMethod", "MortRateMass",
    "MortRateMassScaling", "MortRateRandom", "SurfaceBiomassScaling",
    "SlopeRecruitment", "InterceptRecruitment"
  )
  missing_params <- exptd_params[!exptd_params %in% names(sim_params)]
  if (length(missing_params > 0)) {
    stop(err_msg_missing_params("sim_params", missing_params))
  }
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
  NumberOfSpecies <- nrow(SpeciesPool)  # number of species per 25x25m plot

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
  E <- InitDist
  # Add columns to E for additional info
  E[, c("TotalSurfaceInVoxel", "LightInVoxel", "SurfaceLossInVoxel")] <- 0
  max_id <- nrow(E)  # to trace individual IDs

  isHabitatDynamic <- sim_params$MicrohabitatType == 1
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

    if (isHabitatDynamic) {
      if (!length(Microhabitat) == timeSteps) {
        stop("For dynamic habitats, Microhabitat should have one element for each time step.\n")
      } else {
        # Stash paths for later time steps
        microhab_files <- Microhabitat
      }
    }
    # If all checks ok, read the first one
    Microhabitat <- readRDS(Microhabitat[1])
  }

  #  Convert relative light values to absolute ?mol*m-2*s-1
  Microhabitat[,,,3] <- Microhabitat[,,,3] * sim_params$Imax
  dims <- dim(Microhabitat)
  dimX <- dims[1]
  dimY <- dims[2]
  dimZ <- dims[3]

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

  # Generation loop
  for (t in seq_len(timeSteps)) {

    gen_nb <- sim_params$InitialTimeStep + t - 1 # actual time step

    # Check if the stop criterion is met
    nbIndsAlive <- length(which(E$Status == 1))
    if (nbIndsAlive > StopNbInds) {
      writeLines(paste0(
        "Time ", gen_nb, ": population has exceeded max threshold of ", StopNbInds,
        " individuals. Ending simulation."
        ))
      break
    }

    # Update microhabitat if applicable
    if (isHabitatDynamic && t > 1) {
      Microhabitat <- readRDS(microhab_files[t])
      if (!all.equal(dim(Microhabitat), dims)) {
        stop(
          paste("Invalid microhabitat matrix at time", t,
                ": number of dimensions must be the same as the first matrix")
          )
      }
      Microhabitat[,,,3] <- Microhabitat[,,,3] * sim_params$Imax
    }

    # Update how many species are alive at beginning of generation
    InitialNumberSpecies <- length(unique(E$SpeciesID[E$Status == 1]))
    nbIndsBeforeDispTotal <- length(which(E$Status == 1))

    # Dispersal
    stop("TODO: make sure we use the expanded prob matrix as in old script")
    disp_items <- resolve_repro_dispersal(
      E, Microhabitat, sim_params$SurfaceBiomassScaling,
      expanded_mat_central_point, sim_params$InterceptRecruitment,
      sim_params$SlopeRecruitment,
      prob_disp_matrix,  SpeciesPool, max_id
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
    E <- resolve_mortality(E, SpeciesPool, Microhabitat, sim_params$MortalityMethod,
                           sim_params$MortRateRandom, sim_params$MortRateMass,
                           sim_params$MortRateMassScaling)

    # Mortality due to competition for space
    E <- resolve_competition(E, Microhabitat, sim_params$CompetitionMethod)

    # Age increment
    E$Age <- E$Age + 1

    # Species-level output
    for (sp in seq_len(NumberOfSpecies)) {

      nb_alive <- sum(E$Status == 1 & is_sp, na.rm = TRUE)
      nb_dead_comp <- sum(E$Status == 2 & is_sp, na.rm = TRUE)
      nb_dead_branch <- sum(E$Status == 3 & is_sp, na.rm = TRUE)
      nb_dead_light <- sum(E$Status == 4 & is_sp, na.rm = TRUE)
      nb_dead_base <- sum(E$Status == 5 & is_sp, na.rm = TRUE)
      nb_inds_begin <- nbIndsBeforeDisp[sp]

      row_nb <- (sp - 1) * timeSteps + t
      is_sp <- E$SpeciesID == sp

      sp_output[row_nb, col_sp_t] <- gen_nb
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
    comm_output$timeStep[t] <- gen_nb
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
    comm_output$EpiphyteFilling[t] <- sum(E$Mass^(2/3)) /
      sim_params$SurfaceBiomassScaling / sum(Microhabitat[, , , 1])

    # Command window information
    msg <- "--------------------------------------------"
    msg <- paste_wrap(msg, "Time step: ", gen_nb)
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
    ind_output_file <- sub("*.csv$", paste0("_", t, ".csv"), path_to_ind_output)
    utils::write.csv(E[, inds_output_names()], ind_output_file, row.names = FALSE)

    # Save sp_output for every time step
    sp_output_df <- as.data.frame(sp_output)
    names(sp_output_df) <- sp_output_headers
    utils::write.csv(sp_output_df, path_to_sp_output, row.names = FALSE)

    # Save comm_output for every time step (overwrite old one)
    utils::write.csv(comm_output, path_to_comm_output, append = FALSE, row.names = FALSE)

    # Remove dead individuals from Epimatrix
    E <- E[E$Status <= 1, ]

  } # time loop

}
