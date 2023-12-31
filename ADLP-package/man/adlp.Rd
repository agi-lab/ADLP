% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S3_methods.R, R/adlp.R
\name{print.adlp}
\alias{print.adlp}
\alias{adlp}
\title{Accident and Development period Adjusted Linear Pools (ADLP) Models}
\usage{
\method{print}{adlp}(object, ...)

adlp(components_lst, newdata, partition_func, param_tol = 1e-06, ...)
}
\arguments{
\item{...}{Other named parameters passed onto \link[MM_optim]{MM_optim}}

\item{components_lst}{List of \code{adlp_components}}

\item{newdata}{Validation data to fit the ADLP paritions on}

\item{partition_func}{Partition function used to subset the data. ADLP weights
will be generated for each partition. To specify partition preferences,
set the parameter to \code{adlp_partition_none} if no partitioning is required.
For partitioning the claims triangle by accident periods with predetermined weights,
use \code{adlp_partition_ap}. Alternatively, users can create a custom partition
function by defining the cut-off accident period for each subset manually.}

\item{param_tol}{Tolerance for weights. Any value less than tolerance in
magnitude is assumed zero.}
}
\value{
Object of class \code{adlp}
}
\description{
Class to estimate an ADLP model fitted by Minorization-Maximisation.
}
\details{
See \link[adlp_component]{adlp_component} and \link[adlp_components]{adlp_components}
objects for more information on required format for inputs.

See \link[adlp_partition]{adlp_partition} for information on valid partition
functions.

For an understanding of how partitions affect the performance of the ADLP ensemble,
one might refer to Avanzi, Li, Wong and Xian (2022)
}
\references{
Avanzi, B., Li, Y., Wong, B., & Xian, A. (2022). Ensemble distributional forecasting for insurance loss reserving. arXiv preprint arXiv:2206.08541.
}
