#############################################################################
##                    Datasets Included in This Package                    ##
#############################################################################

#' Claims Data in data.frame Format
#'
#' A data.frame of claims, with the corresponding Accident Period (`origin`)
#' and Development Period (`dev`). A `calendar` column is also included (as the
#' sum of `dev` and `origin`. This format is required for the ADLP package
#'
#' @format A `data.frame` with 4 components:
#' \describe{
#'   \item{origin}{Accident Period}
#'   \item{dev}{Development Period}
#'   \item{calendar}{Accident Period + Development Period}
#'   \item{claims}{Claim amount}
#' }
#' @examples
#' test_claims_dataset$claims
"test_claims_dataset"


#' Test ADLP Component
#'
#' A `adlp_component` object created for examples.
#'
#' @format A `adlp_component` format, see \link[ADLP]{adlp_component}.
#'
#' @examples
#' test_adlp_component
"test_adlp_component"
