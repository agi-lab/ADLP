# ADLP
Repository for ADLP - ensemble reserving package

## Introduction 

We present ADLP (Accident and Development period adjusted Linear Pools), a tailored ensemble technique for general insurance loss reserving. ADLP seeks to combine various loss reserving models, leveraging their strengths, with combination weights optimised to enhance the ensemble's distributional forecasting performance.

This package originates from the paper "Ensemble distributional forecasting for insurance loss reserving," while also offering users ample flexibility to choose or create component models for the ensemble, and to employ data partitioning for calibrating either the component models or the combination weights, aligning with their experiences. 

# Reference

For a full description of ADLP's structure and modelling details, readers should refer to:

Avanzi, B., Li, Y., Wong, B., & Xian, A. (2022). Ensemble distributional forecasting for insurance loss reserving. [arXiv preprint arXiv:2206.08541](https://doi.org/10.48550/arXiv.2206.08541).

To cite this package in publications, please use:

`citation("ADLP")`

# Install Package

To install the development version of the package from this GitHub repository, do

```
if (!require(remotes)) install.packages("remotes")
remotes::install_github("agi-lab/ADLP/ADLP-package", build_vignettes = TRUE)
```

After the installation, run:

`library(ADLP)`

as you would normally do will load the package. View a full demonstration of the package by running

`vignette("ADLP-demo", package = "ADLP")`
