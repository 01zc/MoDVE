test_that("Growth meets expectations", {

  # Initialise state
  dimensions <- sample(2:10, 3, replace = TRUE)
  SurfaceBiomassScaling <- runif(1, 0, 1)
  nb_species <- 3
  repeat {
    SpeciesPool <- create_rnd_species_df(nb_species)
    Microhabitat <- create_rnd_microhabitat(SpeciesPool, dimensions, SurfaceBiomassScaling)
    distr_params <- list(
      "IndividualsPerSpecies" = sample(1:10, 1),
      "ScalingPerHa" = FALSE,
      "PercentageMaturePerSpecies" = runif(1, 0, 1),
      "SurfaceBiomassScaling" = SurfaceBiomassScaling
    )
    E <- draw_initial_individuals(distr_params, SpeciesPool, Microhabitat) |>
      dplyr::filter(Status == 1) # dead individuals don't grow

    # We need at least one live individual of species 1 for the last test
    if (nrow(E[E$SpeciesID == 1,]) > 0) break;
  }

  # Add mass species traits to individual table for reference in tests
  E$max_mass <- E$growth_rate <- 0
  for (i in 1:nrow(E)) {
    which_sp <- which(SpeciesPool$SpeciesID == E$SpeciesID[i])
    E$max_mass[i] <- SpeciesPool$MaximumMass[which_sp]
    E$growth_rate[i] <- SpeciesPool$GrowthRate[which_sp]
  }

  # Save initial table for resets
  E_init <- E

  # Suitability matrix (all perfectly suitable for now)
  SuitabilityMat <- array(1, dim = c(dimensions, nb_species))

  # Growth never decreases mass, or exceed max mass
  mass_never_decreases <- TRUE
  mass_never_above_max <- TRUE

  for (i in 1:1000) {
    mass_before <- E$Mass
    E <- resolve_growth(E,
                        SpeciesPool,
                        Microhabitat,
                        SuitabilityMat,
                        SurfaceBiomassScaling
                        )
    mass_diff <- E$Mass - mass_before
    if (any(mass_diff < -1e-06)) mass_never_decreases <- FALSE
    if (any(E$Mass - E$max_mass > 1e-06)) mass_never_above_max <- FALSE
  }
  expect_true(mass_never_decreases)
  expect_true(mass_never_above_max)

  # Growth scales with suitability
  # Perfect suitability test is implicitly part of above test
  # Unsuitable = no growth
  E <- E_init # reset
  SuitabilityMat[,,,] <- 0
  E <- resolve_growth(
    E, SpeciesPool, Microhabitat, SuitabilityMat, SurfaceBiomassScaling
  )
  expect_equal(E$Mass, E_init$Mass)

  # Intermediate suitability
  E <- E_init[E_init$SpeciesID == 1,]
  rnd_suit <- runif(1)
  SuitabilityMat[,,,] <- rnd_suit
  Microhabitat[,,,3] <- SpeciesPool$OptimumLight[1]
  exptd_mass <- E$Mass + E$growth_rate * (E$max_mass - E$Mass) * rnd_suit
  E <- resolve_growth(
    E, SpeciesPool, Microhabitat, SuitabilityMat, SurfaceBiomassScaling
  )
  expect_equal(E$Mass, exptd_mass)
})

