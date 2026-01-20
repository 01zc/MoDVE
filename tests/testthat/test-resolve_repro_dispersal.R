source("../old_funcs.R")

species_params <- parse_config("../config_a2.toml")

test_that("consistent with old version", {

  nb_species <- 3

  SpeciesPool <- 1:nb_species |>
    purrr::imap(function(i) { draw_species_traits(species_params)}) |>
    purrr::map(as.data.frame) |>
    purrr::list_rbind()

  surface_biomass_scaling <- runif(1, 0, 100)

  # Initialise a random 3D grid with individuals
  dimensions <- sample(1:20, 3)
  surface_area_mat <- array(10, dim = dimensions)
  nb_inds <- sample(1:dimensions[1], 1)

  centralPoint <- find_central_point(dimensions)

  nb_inds <- 20
  #create_rnd_epiphyte_tbl <- function() {}

  E <- tibble::tibble(
    SpeciesID = 1,
    IndividualID = seq_len(nb_inds),
    Status = 1,
    # Initialise individuals randomly
    X = rep(sample(dimensions[1], nb_inds, replace = TRUE)),
    Y = rep(sample(dimensions[2], nb_inds, replace = TRUE)),
    Z = rep(sample(dimensions[3], nb_inds, replace = TRUE)),
    Mass = runif(nb_inds, 0, 1),
    MassAtMaturity = Mass,
    MaxMass = Mass,
    RecruitmentInvestmentRel = runif(nb_inds),
    # Compute expected surface area
    expected_sa = Mass^(2/3) / surface_biomass_scaling
  )

  idx_recruits <- 1:3
  sp <- 1
  E[idx_recruits, names(SpeciesPool)] <- SpeciesPool[sp, ]
  E$X[idx_recruits] <- recruit_coords$x
  E$Y[idx_recruits] <- recruit_coords$y
  E$Z[idx_recruits] <- recruit_coords$z
  E$Mass[idx_recruits] <- 0  # Initial size
  E$Status[idx_recruits] <- 1  # status 1:alive
  recruits_ids <- seq(max_id + 1, max_id + length(totalNbRecruits))
  E$IndividualID[idx_recruits] <- recruits_ids

  Microhabitat <- array(0, c(dimensions, 3))
  Microhabitat[,,,1] <- 0 # surface area
  Microhabiat[,,,3] <- 0 # light

  prob_disp_matrix <- calc_prob_disp_matrix(
    centralPoint, dimensions[1], dimensions[2], dimensions[3],
    SpeciesPool
  )


  disp_list <- old_dispersal(
    nb_species,
    E,
    Microhabitat,
    SurfaceBiomassScaling,
    dimensions,
    centralPoint,
    InterceptRecruitment,
    SlopeRecruitment,
    prob_disp_matrix,
    SpeciesPool,
    max_id
  )

  res_obs <- resolve_repro_dispersal(
    E,
    Microhabitat,
    SurfaceBiomassScaling,
    centralPoint,
    InterceptRecruitment,
    SlopeRecruitment,
    prob_disp_matrix,
    SpeciesPool,
    max_id
  )

})
