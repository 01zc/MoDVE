test_that("reads input", {
  skip() # wip
  # check that input matrices / dfs are read correctly when a path is passed

  sim_params <- parse_config(test_path("configs", "config_a4.toml"))

    # InitDist

  run_modve_sim(
    sim_params,
    SpeciesPool,
    Microhabitat,
    InitDist,
    path_to_ind_output = NULL,
    path_to_sp_output = NULL,
    path_to_comm_output = NULL
  )

})
