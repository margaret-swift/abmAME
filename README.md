# abmAME

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6951937.svg)](thttps://doi.org/10.5281/zenodo.6951937)

--------------------------------------------------------------------------------

abmAME, or "abmAnimalMovementExtension" extends the *[abmAnimalMovement](https://github.com/BenMMarshall/abmAnimalMovement)* package functionality to include:
- Introduction of linear barriers (e.g. fences, roads, rivers) with a perceived permeability for the species/sex/individual being modeled.
- Extension of matrix-based landscapes to include real-world georeferenced rasters in simulated agent step choices.
- [TODO] Hierarchical expansion of shelter use to define both local and global movement restrictions.
- [TODO] Seasonal-based selection of landscape rasters and shelter attraction strength in agent movement choices.

--------------------------------------------------------------------------------

The *abmAnimalMovement* package simulates animal movement use a discrete time agent-based model, programmed in C++ via the Rcpp package. The simulations include a number of key internal and external movement influences, as well as parameters for navigation and mobility capacity of the animal.

A more complete description of the package, alongside a demonstration can be found at DOI: TBC. Or a draft version of that manuscript within the package Github [here](https://github.com/BenMMarshall/abmAnimalMovement/blob/main/notebook/manuscript/Agent-based_model_walkthrough.pdf).

--------------------------------------------------------------------------------

## Installation

**Install the fencing extension to the abmAnimalMovement package**

To install *this package*:

```
install.packages("devtools")
devtools::install_github("margaret-swift/abmAME")
```



## Installation of abmAnimalMovement

**Install from CRAN**

To install [abmAnimalMovement](https://CRAN.R-project.org/package=abmAnimalMovement) from CRAN:

```
install.packages("abmAnimalMovement")

```

**Install with GitHub**

To install the development versions of the abmAnimalMovement package from GitHub, use the `install_github` function from the `devtools` library.

```
install.packages("devtools")
devtools::install_github("BenMMarshall/abmAnimalMovement")
```



## Core simulation function

The `abm_simulate()` function is the main purpose of the package. Guidance on how to parametrise the simulation can be found in accompanying documentation; however, a more detailed walk-through can be found in the preprint *COMING SOON*. <!-- add link here when ready -->
