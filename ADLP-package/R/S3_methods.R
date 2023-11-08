
#' @rdname adlp
#' @export
print.adlp <- function(object, ...) {
    print(object$components_lst)
    print("ADLP weights:")
    print(object$model_weights)
}

#' @rdname adlp_func
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

    mu_index <- match(mu_ij, data_ij)
    ensemble_mu <- mu[mu_index]
    ensemble_mu <- cbind(newdata[, 1:2], ensemble_mu)
    ensemble_mu

}

#' @rdname adlp_component
#' @export
print.adlp_component <- function(object, ...) {
    print(paste("ADLP Component: ", collapse = ""))
    print(object$model_train$call)
}

#' @rdname adlp_components
#' @export
print.adlp_components <- function(object, ...) {
    print("ADLP Components:")
    print(names(object))
}

