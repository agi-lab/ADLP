

#' Custom Model Wrapper
#'
#' Function to define basic functionality needed for a custom model that does
#' not fit the general framework of models that align with \link[ADLP]{adlp_component}
#'
#' @param formula Formula needed that defines all variables required for the model
#' @param data Initial training data for model
#'
#' @return An object of class `custom_model`. `custom_model` is a list that
#' stores the required formula to update the model and the data used to update
#' the model.
#'
#' @details
#' Custom model should support the S3 method `formula` and `update`.
#'
#' @examples
#' data("test_claims_dataset")
#' custom_model <- custom_model(claims~., data=test_claims_dataset)
#'
#'
#' @export
custom_model <- function(formula, data, ...) {
    z <- list(
        formula = stats::formula(formula),
        data = data,
        ...
    )
    class(z) <- "custom_model"
    z
}

#' @rdname custom_model
#' @param object Object of type `custom model`
#' @param data data to update custom model
#' @param ... Additional variables for update
#' @export
update.custom_model <- function(object, data, ...) {
    z <- object
    z$data = data
    z
}
