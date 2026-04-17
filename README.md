# Modeling the Dynamics of Vascular Epiphytes

<!-- badges: start -->
[![R-CMD-check](https://github.com/01zc/MoDVE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/01zc/MoDVE/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

MoDVE runs an individual-based simulation of vascular (i.e., ferns and 
angiosperms) epiphytes growing on exposed tree surfaces in a 3D forest 
environment.

The simulation is an `R` (with `Rcpp` components) implementation of the epiphyte
model introduced in [*Petter et al., 2020*](https://onlinelibrary.wiley.com/doi/10.1002/ece3.7255)
and originally developed in `Java`.

3D tree surfaces are typically generated with a functional-structural forest stand model
[MoF3D](https://github.com/julianoscabral/MoF3D) ([Petter et al. 2021](https://onlinelibrary.wiley.com/doi/10.1002/ece3.7255)), although users may
import surfaces generated from other sources.

Each year, the epiphytes produce and disperse seedlings in the habitat matrix,
undergo growth, mortality and competition, based on user-defined species traits
and local environmental conditions.

# Installation

The package can be installed from this repository:

```r
remotes::install_github("https://github.com/01zc/MoDVE", build_vignettes = TRUE)
```

# Usage

Running the simulation requires three types of input data:
- A microhabitat matrix containing the available surface area for epiphytes to
grow on and local light conditions.
- A table containing species-level traits that defines the community.
- A table containing the position and trait values of initial individuals.

The preparation of each input type is handled by a dedicated function, which 
usage is covered in its corresponding vignette:


| Input | Description | Vignette |
|------------------|------------------|------------------|
| `Microhabitat` | Generation of the microhabitat matrix from `MoF3D` output | `vignette("generate_microhabitat")` | 
| `SpeciesPool` | Generation of the species traits | `vignette("create_species_pool")` |
| `InitDist` | Generation of the initial individuals and their distribution | `vignette("create_initial_distributions")` |

Once the different inputs have been assmbled, running the simulation is straightforward:

```r
run_modve_sim(
    sim_params = config, # simulation controls
    SpeciesPool = SpeciesPool, # Species trait table or path to table
    Microhabitat = Microhabitat, # Habitat matrix, or path to matrix
    InitDist = InitDist, # Initial individuals table, or path to table
    path_to_ind_output = path_to_ind_output, 
    path_to_sp_output = path_to_sp_output,
    path_to_comm_output = path_to_comm_output
  )
```

See the corresponding vignette fora full walkthrough:

```r
vignette("run_model", package = "MoDVE")
```

# Help and bug report

To report a bug, or suggest a feature, please open an issue.

# License

`MoDVE` is made publicly available under the terms of the GNU General Public License v3.

# References

Petter, G., Kreft, H., Ong, Y., Zotz, G., Sarmento, J. (2020). Modeling the 
long-term dynamics of tropical forests: from leaf traits to whole-tree growth 
patterns.
[https://www.sciencedirect.com/science/article/pii/S0304380021002866?via%3Dihub](https://www.sciencedirect.com/science/article/pii/S0304380021002866?via%3Dihub).

Petter, G.; Zotz, G.; Kreft, H.; Sarmento Cabral, J. (2021). Agent-based 
modelling of the effects of forest dynamics, selective logging, and fragment 
size on epiphyte communities. Ecology and Evolution, 11, 2937–2951. 
[https://doi.org/10.1002/ece3.7255](https://doi.org/10.1002/ece3.7255)

