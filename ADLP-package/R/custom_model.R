

#' Custom Model Wrapper
#'
#' Function to define basic functionality needed for a custom model that does
#' not fit the general framework of models that align with \link[adlp_component]{adlp_component}
#'
#' @param formula Formula needed that defines all variables required for the model
#' @param data Initial training data for model
#'
#' @details
#' Custom model should support the S3 method `formula` and `update`.
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
#' @export
update.custom_model <- function(x, data) {
    z <- x
    z$data = data
    z
}
