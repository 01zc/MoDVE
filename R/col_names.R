species_trait_names <- function() {
  return(c("MaximumMass", "MassAtMaturity", "GrowthRate",
           "DispersalKernel", "DispersalKernelAsymmetry",
           "RecruitmentInvestmentRel", "RecruitmentInc",
           "MinLight", "MaxLight", "OptimumLight",
           "LightResponseA", "LightResponseB", "LightResponseC"))
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

species_output_names <- function() {
  return(c(
    "TimeStep", "SpeciesID", "NumberIndividualsBeginning",
    "NumberIndividualsEnd", "NumberMatureIndividuals", "NumberRecruits",
    "NumberRecruitsPotential", "NumberMortalityBranchFall",
    "NumberMortalityLight", "NumberMortalityCompetition",
    "NumberMortalityNatural", "PopulationGrowthRate",
    "PopulationGrowthRateLog", "BirthRate", "DeathRate", "AverageMass",
    "AverageAge", "MinLight", "MaxLight", "MeanLight",
    "NumberMortalityHum", "NumberMortalityTemp", "NumberMortalityWind",
    "MinHum", "MaxHum", "MeanHum",
    "MinTemp", "MaxTemp", "MeanTemp",
    "MinWind", "MaxWind", "MeanWind"
  ))
}

comm_output_names <- function() {
  return(c(
    "timeStep", "NumberSpeciesBeginning", "NumberSpeciesEnd",
    "NumberIndividualsBeginning", "NumberIndividualsEnd", "nb_recruits_matrix",
    "MortalityBranchFall", "MortalityLight", "MortalityCompetition",
    "MortalityNatural", "MortalityHum", "MortalityTemp",
    "MortalityWind", "BranchSurfaceIndex", "EpiphyteFilling"
  ))
}
