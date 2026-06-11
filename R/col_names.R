species_trait_names <- function(microclimate_opts) {
  sptr_names <- c("MaximumMass", "MassAtMaturity", "GrowthRate",
                  "DispersalKernel", "DispersalKernelAsymmetry",
                  "RecruitmentInvestmentRel", "RecruitmentInc",
                  "MinLight", "MaxLight", "OptimumLight",
                  "LightResponseA", "LightResponseB", "LightResponseC")
  if (microclimate_opts$use_humidity)
    sptr_names <- c(sptr_names, "MinHum", "MaxHum", "OptimumHum")
  if (microclimate_opts$use_temperature)
    sptr_names <- c(sptr_names, "MinTemp", "MaxTemp", "OptimumTemp")
  if (microclimate_opts$use_wind)
    sptr_names <- c(sptr_names, "MinWind", "MaxWind", "OptimumWind",
                    "DispersalKernelWindEffect")
  return(sptr_names)
}

inds_input_names <- function() {
  return(c("X", "Y", "Z", "Mass", "Status", "IndividualID",
           "SurfaceAreaOccupied", "Age", "SpeciesID"))
}

inds_output_names <- function() {
  return(c(
    "SpeciesID", "IndividualID", "Status", "Mass", "Age", "X", "Y", "Z",
    "TotalSurfaceInVoxel", "SurfaceLossInVoxel", "LightInVoxel", "HumInVoxel",
    "TempInVoxel", "WindInVoxel"
    ))
}

species_output_names <- function(microclimate_opts) {
  sp_out_names <- c(
    "TimeStep", "SpeciesID", "NumberIndividualsBeginning",
    "NumberIndividualsEnd", "NumberMatureIndividuals", "NumberRecruits",
    "NumberRecruitsPotential", "NumberMortalityBranchFall",
    "NumberMortalityLight", "NumberMortalityCompetition",
    "NumberMortalityNatural", "PopulationGrowthRate",
    "PopulationGrowthRateLog", "BirthRate", "DeathRate", "AverageMass",
    "AverageAge", "MinLight", "MaxLight", "MeanLight"
  )
  if (microclimate_opts$use_humidity)
    sp_out_names <- c(sp_out_names,  "NumberMortalityHum", "MinHum", "MaxHum", "MeanHum")
  if (microclimate_opts$use_temperature)
    sp_out_names <- c(sp_out_names, "NumberMortalityTemp", "MinTemp", "MaxTemp", "MeanTemp")
  if (microclimate_opts$use_wind)
    sp_out_names <- c(sp_out_names, "NumberMortalityWind", "MinWind", "MaxWind", "MeanWind")
  return(sp_out_names)
}

comm_output_names <- function(microclimate_opts) {
  comm_out_names <- c(
    "timeStep", "NumberSpeciesBeginning", "NumberSpeciesEnd",
    "NumberIndividualsBeginning", "NumberIndividualsEnd", "nb_recruits_matrix",
    "MortalityBranchFall", "MortalityLight", "MortalityCompetition",
    "MortalityNatural", "MortalityHum", "MortalityTemp",
    "MortalityWind", "BranchSurfaceIndex", "EpiphyteFilling"
  )
  if (microclimate_opts$use_humidity)
    comm_out_names <- c(comm_out_names,  "MortalityHum")
  if (microclimate_opts$use_temperature)
    comm_out_names <- c(comm_out_names, "MortalityTemp")
  if (microclimate_opts$use_wind)
    comm_out_names <- c(comm_out_names, "MortalityWind")
  return(comm_out_names)
}
