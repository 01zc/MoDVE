
#' Functions from previous versions of MoDVE
#'
#' These functions are kept here to test that new versions remain consistent
#' with previous versions


#' Find which voxels are intersected by a branch segment
#'
#' This function implements the solution that was used in the original
#' Matlab implementation of MoDVE. The function uses a simple and quick
#' rule, but the solution is incorrect. It has been replaced with the Bresenham
#' algorithm (`find_intersecting_voxels()`).
#'
#' The original version is kept here to help testing that new versions still
#' produce results consistent with the older versions.
#'
#' @param seg_start a numeric vector of length 3 with the x-y-z coordinates
#' of the start point of the segment
#' @param seg_end a numeric vector of length 3 with the x-y-z coordinates
#' of the end point of the segment
#'
#' @returns a list of length 3 integer vectors, the 3D coordinates of all voxels
#' intersected by the segment

old_find_intersecting_voxels <- function(seg_start, seg_end) {

  intersctd_voxels <- list()

  UniqueX <- c(ceiling(seg_start[1]), ceiling(seg_end[1]))
  UniqueY <- c(ceiling(seg_start[2]), ceiling(seg_end[2]))
  UniqueZ <- c(ceiling(seg_start[3]), ceiling(seg_end[3]))

  numX <- sum((UniqueX[2] - UniqueX[1]) > 0, na.rm=TRUE) + 1
  numY <- sum((UniqueY[2] - UniqueY[1]) > 0, na.rm=TRUE) + 1
  numZ <- sum((UniqueZ[2] - UniqueZ[1]) > 0, na.rm=TRUE) + 1

  for (x in seq_len(numX)) {
    xid <- UniqueX[x]

    for (y in seq_len(numY)) {
      yid <- UniqueY[y]

      for (z in seq_len(numZ)) {
        zid <- UniqueZ[z]
        intersctd_voxels <- c(intersctd_voxels, list(c(xid, yid, zid)))
      }
    }
  }
  return(intersctd_voxels)
}

old_compute_prob_matrix_norm <- function(centralPoint, dimX, dimY, dimZ, NumberOfSpecies, SpeciesPool) {
  # Erzeugen der Distanzmatrix mit allen Distanzen zum
  DistanceMatrix <- array(rep(0, dimX * dimY * dimZ), dim=c(dimX, dimY, dimZ))

  for (i in seq_len(dimX)) {
    for (j in seq_len(dimY)) {
      for (k in seq_len(dimZ)) {
        x1 <- c(i, j, k)
        x2 <- c(centralPoint[1], centralPoint[2], centralPoint[3])
        DistanceMatrix[i, j, k] <- sqrt(sum((x1 - x2)^2))  # call to pdist() in the matlab script
      }
    }
  }

  # Erzeugen der Wahrschienlichkeitsmatrix anhand der Distanzmatrix und dem
  # artspezischien Wert bb aus der Epiphytenmatrix
  # Google translate: Generating the probability matrix based on the distance
  # matrix and the species value from the epiphyte matrix
  # negExp = @(distance,bb) exp(-distance.*bb); %Negative Exponential function
  ProbabilityMatrix <- array(rep(0, dimX * dimY * dimZ * NumberOfSpecies), dim=c(dimX, dimY, dimZ, NumberOfSpecies))
  ProbabilityMatrixNormalized <- array(rep(0, dimX * dimY * dimZ * NumberOfSpecies), dim=c(dimX, dimY, dimZ, NumberOfSpecies))

  for (i in seq_len(NumberOfSpecies)) {
    exponentE <- SpeciesPool$DispersalKernel[i]
    dispersalAsymmetry <- SpeciesPool$DispersalKernelAsymmetry[i]

    ProbabilityMatrix[, , , i] <- exp(-DistanceMatrix * exponentE)  # call to negExp(DistanceMatrix(:,:,:),exponentE) in matlab

    # Apply dispersal asymmetry (probability to disperse downwards higher than upwards dispersal)
    # WARNING: The index i in the 4th dimension was missing from the Matlab script and so we were
    # getting "incorrect number of dimensions" in R. Added it here but need to ask if this is what
    # it was supposed to do.
    id3_1 <- int_seq(from=centralPoint[3], to=dimZ, by=1)
    id3_2 <- int_seq(from=1, to=centralPoint[3] - 1, by=1)
    ProbabilityMatrix[, , id3_1, i] <- ProbabilityMatrix[, , id3_1, i] * ((1 - dispersalAsymmetry) / 0.5)
    ProbabilityMatrix[, , id3_2, i] <- ProbabilityMatrix[, , id3_2, i] * (dispersalAsymmetry / 0.5)

    ProbabilityMatrixNormalized[, , , i] <- ProbabilityMatrix[, , , i] / sum(ProbabilityMatrix[, , , i])  # sum(sum(sum(ProbabilityMatrix(:,:,:,i)))) in matlab
  }

  return(ProbabilityMatrixNormalized)
}

old_dispersal <- function(NumberOfSpecies,
                      E,
                      Microhabitat,
                      SurfaceBiomassScaling,
                      dimPlot,
                      centralPoint,
                      InterceptRecruitment,
                      SlopeRecruitment,
                      ProbabilityMatrixNormalized,
                      SpeciesPool,
                      MaxIndividualID) {
  # Store number of individuals at beginning of time step
  IntialNumberIndividuals <- array(rep(0, NumberOfSpecies))
  for (g in seq_len(NumberOfSpecies)) {
    # Count indices where SpeciesID is g and Status is 1
    IntialNumberIndividuals[g] <- length(which(E$SpeciesID == g & E$Status == 1))
  }
  IntialNumberIndividualsTotal <- length(which(E$Status == 1))
  InitialNumberSpecies <- length(unique(E$SpeciesID[E$Status == 1]))
  NumberRecruitsPerSpecies <- array(rep(0, NumberOfSpecies))

  # Calculate free surface area per voxel
  AvailableSurfaceArea <- Microhabitat[, , , 1]
  for (i in seq_len(nrow(E))) {
    SurfaceAreaNeededInVoxel <- E$Mass[i]^(2/3) / SurfaceBiomassScaling
    AvailableSurfaceArea[E$X[i], E$Y[i], E$Z[i]] <- max(0, AvailableSurfaceArea[E$X[i], E$Y[i], E$Z[i]] - SurfaceAreaNeededInVoxel)
  }

  # Check if there are species left (~isempty(E) in matlab)
  if (nrow(E) > 0) {
    unique_species <- unique(E$SpeciesID)  # list with species IDs of all present species

    # Initialize potential recruitment dataframe
    PotentialRecruitment <- data.frame(matrix(0, nrow=length(unique_species), ncol=2))
    colnames(PotentialRecruitment) <- c("index", "potential_recruit")

    # loop over all species
    for (i in seq_len(length(unique_species))) {
      # Generate initially empty matrix to store the probabilities for recruitment
      ProbabilityMatrixPerSpecies <- array(rep(0, dimPlot[1] * dimPlot[2] * dimPlot[3]), dim=c(dimPlot[1], dimPlot[2], dimPlot[3]))

      # Matrix containing all mature individuals of one species
      MatureIndividulsPerSpecies <- E[E$SpeciesID == unique_species[i] & E$Mass >= E$MassAtMaturity, ]

      # ~isempty(MatureIndividulsPerSpecies) in matlab
      if (nrow(MatureIndividulsPerSpecies) > 0) {

        # Probability matrix for each species: Depending on the position of each mature individual,
        # the total probability for the species is calculated.
        # The second part of the equation accounts for the actual size of the individual
        # in relation to the maximum size for which the recruitment per individual is defined
        for (j in seq_len(nrow(MatureIndividulsPerSpecies))) {
          idx1 <- seq(from=centralPoint[1] - MatureIndividulsPerSpecies$X[j] + 1, to=centralPoint[1] - MatureIndividulsPerSpecies$X[j] + dimPlot[1], by=1)
          idx2 <- seq(from=centralPoint[2] - MatureIndividulsPerSpecies$Y[j] + 1, to=centralPoint[2] - MatureIndividulsPerSpecies$Y[j] + dimPlot[2], by=1)
          idx3 <- seq(from=centralPoint[3] - MatureIndividulsPerSpecies$Z[j] + 1, to=centralPoint[3] - MatureIndividulsPerSpecies$Z[j] + dimPlot[3], by=1)
          idx4 <- MatureIndividulsPerSpecies$SpeciesID[j]

          factor1 <- (InterceptRecruitment + (SlopeRecruitment * MatureIndividulsPerSpecies$Mass[j])) * MatureIndividulsPerSpecies$RecruitmentInvestmentRel[j]
          factor2 <- (MatureIndividulsPerSpecies$Mass[j] - MatureIndividulsPerSpecies$MassAtMaturity[j]) / (MatureIndividulsPerSpecies$MaximumMass[j] - MatureIndividulsPerSpecies$MassAtMaturity[j])
          factor3 <- 1 + (MatureIndividulsPerSpecies$RecruitmentInc[j] * factor2)

          ProbabilityMatrixPerSpecies <- ProbabilityMatrixPerSpecies +
            ProbabilityMatrixNormalized[idx1, idx2, idx3, idx4] *
            factor1 * factor3

        }

        # Store potential normalized number of recruits. We will use this to populate SummaryMatrixSpecies later
        PotentialRecruitment$index[i] <- i
        PotentialRecruitment$potential_recruit[i] <- sum(ProbabilityMatrixPerSpecies)  # potential recruitment / sum(sum(sum(ProbabilityMatrixPerSpecies))) in matlab

        # Matix containing all voxel for which the light requirements are fulfilled
        # We use the first row from MatureIndividulsPerSpecies. Since its elements have the same SpeciesID
        # then the MinLight and MaxLight is the same for all rows.
        pot_habitat <- ifelse((Microhabitat[, , , 3] >= MatureIndividulsPerSpecies$MinLight[1]) & (Microhabitat[, , , 3] <= MatureIndividulsPerSpecies$MaxLight[1]), 1, 0)

        # Final probabiliy matrix for new recruits
        probability_recruits <- ProbabilityMatrixPerSpecies *
          pot_habitat * AvailableSurfaceArea

        # Calculate number of recuits based on final probability matrix
        Recruits <- array(
          rpois(length(probability_recruits), probability_recruits),
          dim=dim(probability_recruits)
          )  # poissrnd(probability_recruits) in matlab

        # Add new recruits to epiphyte matrix
        num_recruits <- sum(Recruits)  # sum(sum(sum(Recruits))) in matlab
        NumberRecruitsPerSpecies[unique_species[i]] <- num_recruits

        if (num_recruits > 0) {
          ids <- arrayInd(which(Recruits > 0), dim(Recruits))
          xInd <- ids[, 1]
          yInd <- ids[, 2]
          zInd <- ids[, 3]

          while (num_recruits > length(xInd)) {
            tmp_ids <- arrayInd(which(Recruits > 0), dim(Recruits))
            Recruits[tmp_ids] = Recruits[tmp_ids] - 1

            tmp_ids <- arrayInd(which(Recruits > 0), dim(Recruits))
            xInd <- append(xInd, tmp_ids[, 1])
            yInd <- append(yInd, tmp_ids[, 2])
            zInd <- append(zInd, tmp_ids[, 3])
          }
          vec_recruits <- seq(from=nrow(E) + 1, to=nrow(E) + length(xInd), by=1)

          # Copy species information to Epiphyte matrix
          E[vec_recruits, names(SpeciesPool)] <- SpeciesPool[unique_species[i], ]
          E$X[vec_recruits] <- xInd
          E$Y[vec_recruits] <- yInd
          E$Z[vec_recruits] <- zInd
          E$Mass[vec_recruits] <- 0  # Initial size
          E$Status[vec_recruits] <- 1  # status 1:alive
          E$IndividualID[vec_recruits] <- seq(from=MaxIndividualID + 1, to=MaxIndividualID + length(xInd), by=1)  # individual ID

          E[is.na(E)] <- 0  # convert all NA to 0 so that the R script matches the Matlab

          MaxIndividualID <- MaxIndividualID + length(xInd)
        }
      }
    }
  } else {
    # Initialize empty potential recruitment dataframe
    PotentialRecruitment <- data.frame(matrix(0, nrow=0, ncol=2))
    colnames(PotentialRecruitment) <- c("index", "potential_recruit")
  }

  disp_items <- list(
    "IntialNumberIndividuals" = IntialNumberIndividuals,
    "NumberRecruitsPerSpecies" = NumberRecruitsPerSpecies,
    "InitialNumberSpecies" = InitialNumberSpecies,
    "IntialNumberIndividualsTotal" = IntialNumberIndividualsTotal,
    "E" = E,
    "PotentialRecruitment" = PotentialRecruitment,
    "MaxIndividualID" = MaxIndividualID
  )
  return(disp_items)
}


