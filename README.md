# ADLP-sandbox
Sandbox repository for ADLP - ensemble reserving package

## Introduction 

We present ADLP (Accident and Development period adjusted Linear Pools), a tailored ensemble technique for general insurance loss reserving. ADLP seeks to combine various loss reserving models, leveraging their strengths, with combination weights optimised to enhance the ensemble's distributional forecasting performance.

This package originates from the paper "Ensemble distributional forecasting for insurance loss reserving," while also offering users ample flexibility to choose or create component models for the ensemble, and to employ data partitioning for calibrating either the component models or the combination weights, aligning with their experiences. 

## Package Overview

This section provides an overview of the folders and files located in this repository; their purposes will also be briefly introduced. 

* R: stores the sources of R codes used in constructing the package functions.
    * `train_val_split.R`
        * Defines the functions for partitioning the claims triangle into training and validation sets.
    * `components.R`
        * Defines the functions for storing the component models used in the ensemble, and functions for calculating the density, mean, and cumulative distribution of the component models. The simulation function for component models is also contained. 
    * `custom_model.R`:
        * Defines the functions to build customised models.
    * `partitions.R`:
        * Defines the functions to partition the data used for calibrating the ensemble weights.
    * `mm_optim.R`:
        * Defines the functions to optimise the ensemble weights based on the Minorisation-Maximisation (MM) algorithm.
    * `adlp.R`:
        * Defines the functions to calibrate an ADLP ensemble, and functions to calculate the density, Log Score and CRPS of the fitted ADLP objects. The simulation function for ADLP ensembles is also contained.
    * `S3_methods.R`:
        * Contains miscellaneous functions used for predictions and results printing.
 * DemoData: stores demomonstration datasets used for the ADLP package demonstration (`ADLP-demo.Rmd`).
 * vignettes: contains the demonstration file for the ADLP package (`ADLP-demo.Rmd`).




# Reference

For a full description of ADLP's structure and modelling details, readers should refer to:

Avanzi, B., Li, Y., Wong, B., & Xian, A. (2022). Ensemble distributional forecasting for insurance loss reserving. [arXiv preprint arXiv:2206.08541](https://doi.org/10.48550/arXiv.2206.08541).

To cite this package in publications, please use:

`citation("ADLP")`

# Install Package

To install the development version of the package from this GitHub repository, do

```
if (!require(remotes)) install.packages("remotes")
remotes::install_github("agi-lab/ADLP", build_vignettes = TRUE)
```

After the installation, run:

`library(ADLP)`

as you would normally do will load the package. View a full demonstration of the package by running

`vignette("ADLP-demo", package = "ADLP")`
