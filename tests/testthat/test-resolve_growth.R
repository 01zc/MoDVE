source("../test-utils.R")

test_that("Growth meets expectations", {

  dimensions <- sample(1:10, 3, replace = TRUE)
  SurfaceBiomassScaling <- runif(1, 0, 1)
  SpeciesPool <- create_rnd_species_df(5)
  Microhabitat <- create_rnd_microhabitat(SpeciesPool, dimensions, SurfaceBiomassScaling)
  distr_params <- list(
    "IndividualsPerSpecies" = sample(1:10, 1),
    "ScalingPerHa" = FALSE,
    "PercentageMaturePerSpecies" = runif(1, 0, 1),
    "SurfaceBiomassScaling" = SurfaceBiomassScaling
  )
  E <- draw_initial_individuals(distr_params, SpeciesPool, Microhabitat)
  E$max_mass <- E$growth_rate <- 0
  for (i in 1:nrow(E)) {
    which_sp <- which(SpeciesPool$SpeciesID == E$SpeciesID[i])
    E$max_mass[i] <- SpeciesPool$MaximumMass[which_sp]
    E$growth_rate[i] <- SpeciesPool$GrowthRate[which_sp]
  }
  E_init <- E

  mass_never_decreases <- TRUE
  mass_never_abov_max <- TRUE
  for (i in 1:1000) {
    mass_before <- E$Mass
    E <- resolve_growth(
      E, SpeciesPool, Microhabitat, SurfaceBiomassScaling
    )
    mass_diff <- E$Mass - mass_before
    if (any(mass_diff < 0)) mass_never_decreases <- FALSE
    if(any(E$Mass > E$max_mass)) mass_never_abov_max <- FALSE
  }
  expect_true(mass_never_decreases)
  expect_true(mass_never_abov_max)

  # No growth outside of light niche
  E <- E_init # reset
  if (any(SpeciesPool$MinLight < 1))
  SpeciesPool$MinLight <- pmax(1, SpeciesPool$MinLight)
  Microhabitat[,,,3] <- min(SpeciesPool$MinLight) / 2
  E <- resolve_growth(
    E, SpeciesPool, Microhabitat, SurfaceBiomassScaling
  )
  expect_equal(E$Mass, E_init$Mass)

  # At optimum light, growth is equal to k * (max_mass - mass)
  Microhabitat[,,,3] <- SpeciesPool$OptimumLight[1]
  E <- E[E$SpeciesID == 1,]
  exptd_mass <- E$Mass + E$growth_rate * (E$max_mass - E$Mass)
  E <- resolve_growth(
    E, SpeciesPool, Microhabitat, SurfaceBiomassScaling
  )
  expect_equal(E$Mass, exptd_mass)
})


