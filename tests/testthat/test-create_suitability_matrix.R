test_that("light niche", {

  # Piece of UT salvaged from the growth test, use this as a starting point?

  # Create
  E <- E_init # reset
  Microhabitat[,,,3] <- min(SpeciesPool$MinLight) / 2
  E <- resolve_growth(
    E, SpeciesPool, Microhabitat, SuitabilityMat, SurfaceBiomassScaling
  )
  expect_equal(E$Mass, E_init$Mass)

  # At optimum light, growth is equal to k * (max_mass - mass)
  # We carry this test only for the first species
  E <- E_init[E_init$SpeciesID == 1,]
  Microhabitat[,,,3] <- SpeciesPool$OptimumLight[1]
  exptd_mass <- E$Mass + E$growth_rate * (E$max_mass - E$Mass)
  E <- resolve_growth(
    E, SpeciesPool, Microhabitat, SuitabilityMat, SurfaceBiomassScaling
  )
  expect_equal(E$Mass, exptd_mass)
})
