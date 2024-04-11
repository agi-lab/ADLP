
#' @rdname adlp
#' @param x Object of class `adlp`
#' @export
print.adlp <- function(x, ...) {
    print(x$components_lst)
    print("ADLP weights:")
    print(x$model_weights)
}

#' @rdname adlp_func
#'
#' @param object Object of class `adlp`
#' @param newdata new data for prediction. Defaults to NULL
#' @param ... Other parameters to pass onto predict
#'
#' @details
#' Predicts the central estimates based on the ADLP component models and weights.
#'
#' @export
predict.adlp <- function(object, newdata = NULL, ...) {

    if (is.null(newdata)) {
        newdata <- object$data
    }

    component_mu = calc_adlp_component_lst(
        object$components_lst, newdata, "full", "mu"
    )
    mu_partitions <- object$partition_func(component_mu)
    n.partitions <- length(mu_partitions)
    ensemble_w <- object$model_weights
    mu <- c()
    mu_ij <- c()
    data_ij <- paste(newdata$origin, newdata$dev, sep = "-")

    for (j in 1:n.partitions) {
        meta_partition <- mu_partitions[[j]][, -c(1, 2)]
        mu_predict <- as.matrix(meta_partition) %*% ensemble_w[[j]]
        mu <- c(mu, mu_predict)
        mu_ij <- c(mu_ij, paste(mu_partitions[[j]]$origin, "-", mu_partitions[[j]]$dev, sep = ""))
    }

    mu_index <- match(data_ij, mu_ij)
    ensemble_mu <- mu[mu_index]
    ensemble_mu <- cbind(newdata[, 1:2], ensemble_mu)
    ensemble_mu

}

#' @rdname adlp_component
#' @param x Object of class `adlp_component`
#' @export
print.adlp_component <- function(x, ...) {
    print(paste("ADLP Component: ", collapse = ""))
    print(x$model_train$call)
}

#' @rdname adlp_components
#' @param x Object of class `adlp_components`
#' @export
print.adlp_components <- function(x, ...) {
    print("ADLP Components:")
    print(names(x))
}

