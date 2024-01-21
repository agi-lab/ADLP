
#' Accident and Development period Adjusted Linear Pools (ADLP) Models
#'
#' Class to estimate an ADLP model fitted by Minorization-Maximisation.
#'
#' @param components_lst List of `adlp_components`
#' @param newdata Validation data to fit the ADLP paritions on
#' @param partition_func Partition function used to subset the data. ADLP weights
#' will be generated for each partition. To specify partition preferences,
#' set the parameter to `adlp_partition_none` if no partitioning is required.
#' For partitioning the claims triangle by accident periods with predetermined weights,
#' use `adlp_partition_ap`. Alternatively, users can create a custom partition
#' function by defining the cut-off accident period for each subset manually.
#' @param param_tol Tolerance for weights. Any value less than tolerance in
#' magnitude is assumed zero.
#' @param ... Other named parameters passed onto \link[MM_optim]{MM_optim}
#'
#' @return Object of class `adlp`
#'
#' @details
#' See \link[adlp_component]{adlp_component} and \link[adlp_components]{adlp_components}
#' objects for more information on required format for inputs.
#'
#' See \link[adlp_partition]{adlp_partition} for information on valid partition
#' functions.
#'
#' For an understanding of how partitions affect the performance of the ADLP ensemble,
#' one might refer to Avanzi, Li, Wong and Xian (2022)
#'
#' @references Avanzi, B., Li, Y., Wong, B., & Xian, A. (2022). Ensemble distributional forecasting for insurance loss reserving. arXiv preprint arXiv:2206.08541.
#' @export
adlp <- function(
    components_lst, newdata, partition_func, param_tol = 1e-16, ...
) {

    # Calculate densities
    component_dens = calc_adlp_component_lst(
        components_lst, newdata, model = "train", calc = "pdf"
    )

    n.components <- length(components_lst)
    w_init<-rep(1/n.components, n.components)
    model_weights<-list()
    optim_MM_par <- list()

    valid_partitions <- partition_func(component_dens, ...)
    n.partitions <- length(valid_partitions)

    for (i in 1:n.partitions) {

        # Transform the partitions into a data.frame structure for processing.
        partition_in <- do.call('rbind', valid_partitions[1:i])
        dens_in <- partition_in[, 3:(n.components+2)]

        #Train the model weights using the MM Algorithm
        optim_MM <- MM_optim(
            w_init = w_init,
            dat = as.matrix(dens_in),
            ...
        )
        finalw_MM <- optim_MM$finalparams
        finalw_MM <- ifelse(abs(finalw_MM) < param_tol, 0, finalw_MM)
        finalw_MM <- finalw_MM/sum(finalw_MM)

        #Calculate the predictive density by the ensemble
        model_weights[[i]] <- finalw_MM
        optim_MM_par[[i]] <- optim_MM
    }

    z <- list(
        components_lst = components_lst,
        model_weights = model_weights,
        partition_func = partition_func,
        optim_MM = optim_MM,
        newdata = newdata
    )
    class(z) <- "adlp"

    return (z)
}

#' Accident and Development period Adjusted Linear Pools (ADLP) Functions
#'
#' @description
#' Family of functions used to support ADLP inference and prediction.
#'
#' @return `data.frame` of results, where the first and second columns correspond
#' to the `$origin` and `$dev` columns from the triangles. An index column for
#' `simulation #` is also included when simulating ADLP.
#'
#' @name adlp_func
NULL

#' @rdname adlp_func
#'
#' @param adlp Object of class `adlp`
#' @param newdata Data to perform the function on
#' @param model Whether the `train` or `full` model should be used in function
#'
#' @details
#' Calculates the probability density ad each point, given `newdata`.
#'
#' @export
adlp_dens <- function(adlp, newdata, model = c("train", "full")) {

    model <- match.arg(model)
    component_dens = calc_adlp_component_lst(adlp$components_lst, newdata, model, "pdf")

    dens_partitions <- adlp$partition_func(component_dens)
    n.partitions <- length(dens_partitions)
    ensemble_w <- adlp$model_weights
    dens <- c()
    dens_ij <- c()
    data_ij <- paste(newdata$origin, newdata$dev, sep = "-")

    for (j in 1:n.partitions) {

        meta_partition <- dens_partitions[[j]][, -c(1, 2)]
        dens_predict <- as.matrix(meta_partition) %*% ensemble_w[[j]]
        dens <- c(dens, dens_predict)
        dens_ij <- c(dens_ij, paste(dens_partitions[[j]]$origin, "-", dens_partitions[[j]]$dev, sep = ""))
    }

    dens_index <- match(data_ij, dens_ij)
    ensemble_dens <- dens[dens_index]
    ensemble_dens <- cbind(newdata[, 1:2], ensemble_dens)
    ensemble_dens
}

#' @rdname adlp_func
#'
#' @param adlp Object of class `adlp`
#' @param newdata Data to perform the function on
#' @param model Whether the `train` or `full` model should be used in function
#' @param epsilon Offset added to the density before calculating the log
#'
#' @details
#' Calculates the log score, which is the log of the probability density, with
#' an offset `epsilon` to handle zero densities.
#'  Log Score is a strictly proper scoring rule.
#' For full discussion of the mathematical details and
#' advantages of Log Score, one might refer to Gneiting and Raftery (2007)
#'
#' @export
adlp_logS <- function(adlp, newdata, model = c("train", "full"), epsilon = 1e-6) {
    model <- match.arg(model)
    component_dens <- adlp_dens(adlp, newdata, model)
    dens_index <- component_dens[, 1:2]
    dens_val <- component_dens[, -(1:2)]
    dens_val <- log(dens_val + epsilon)
    z <- cbind(dens_index, dens_val)
    z
}

#' @rdname adlp_func
#'
#' @param adlp Object of class `adlp`
#' @param newdata Data to perform the function on
#' @param response_name The column name of the response variable; in string format
#' @param model Whether the `train` or `full` model should be used in function
#' @param lower The lower limit to calculate CRPS; the default value is set to be 1
#' @param upper The upper limit to calculate CRPS; the default value is set to be
#' twice the maximum value of the response variable in the dataset
#' @param sample_n The number of evenly spaced values to sample between lower and upper range of
#' numeric integration used to calculate CRPS. This sample function is designed to
#' constrain memory usage during the computation of CRPS,
#' particularly when dealing with large response variables.
#'
#' @details
#' Continuously Ranked Probability Score (CRPS) is calculated for each data point.
#' `lower` and `upper` are used as limits when approximating the integral.
#' CRPS is a strictly proper scoring rule.
#' For full discussion of the mathematical details and
#' advantages of CRPS, one might refer to Gneiting and Raftery (2007).
#' The CRPS function has been discretized in this context to ensure
#' adaptability to various distributions.
#' For details, one might refer to
#' Gneiting and Ranjan (2011)
#'
#' @references
#' Gneiting, T., Raftery, A. E., 2007. Strictly proper scoring rules, prediction, and estimation. Journal of the American Statistical Association 102 (477), 359–378.
#'
#' Gneiting, T., Ranjan, R., 2011. Comparing density forecasts using threshold-and quantile-weighted scoring rules. Journal of Business & Economic Statistics 29 (3), 411–422.
#'
#' @export
adlp_CRPS <- function(adlp, newdata, response_name, model = c("train", "full"), lower = 1, upper=NULL, sample_n = 2000) {

    model <- match.arg(model)

    # Default upper is 2x the largest response in model.frame

    response_y <- newdata[, response_name]


    if (is.null(upper)) {
         upper <- round(2*max(response_y),0)
    }

    # Sample evenly spaced values between the lower and upper bound by using the quantile function:
    z <- floor(quantile(lower:upper, seq(0, 1, by = 1/sample_n)))[-1]

    z <- sort(z)

    # Calculate CDF as at lower:upper for each data point
    component_cdfs <- calc_adlp_component_lst(
        adlp$components_lst,
        newdata,
        model,
        "cdf",
        y = rep(z, each = nrow(newdata))
    )

    component_cdfs_lst <- split(component_cdfs, rep(z, each = nrow(newdata)))

    # Calculate CRPS for adlp
    I<-function(y, z) ifelse(y<=z,1,0)

    out_partitions <- adlp$partition_func(newdata)
    n.partitions <- length(out_partitions)

    ensemble_w <- adlp$model_weights
    crps <- c()
    crps_ij <- c()
    data_ij <- paste(newdata$origin, newdata$dev, sep = "-")

    for (k in 1:n.partitions) {
        partition_ind <- rownames(newdata) %in% rownames(out_partitions[[k]])
        w <- ensemble_w[[k]]


        y = out_partitions[[k]][, response_name]

        pred_CDF_ensemble <- lapply(
            component_cdfs_lst,
            function (x) {
                as.matrix(x[partition_ind, -(1:2)]) %*% w
            }
        )

        pred_vs_obs_cdf <- unlist(pred_CDF_ensemble) -
            I(rep(y, length(z)), rep(z, each = length(y)))

        crps_ensemble <- unlist(lapply(
            split(pred_vs_obs_cdf^2, rep(1:length(y), length(z))),
            sum
        ))

        crps <- c(crps, crps_ensemble)
        crps_ij <- c(crps_ij, paste(out_partitions[[k]]$origin, "-", out_partitions[[k]]$dev, sep = ""))
    }

    crps_index <- match(data_ij, crps_ij)
    ensemble_crps <- crps[crps_index]
    ensemble_crps <- cbind(newdata[, 1:2], ensemble_crps)
    ensemble_crps
}

#' @rdname adlp_func
#'
#' @param n number of simulations
#'
#' @details
#' Simulations of ADLP predictions, given component models and ADLP weights.
#'
#' @export
adlp_simulate <- function(n, adlp, newdata = NULL) {

    if (is.null(newdata)) {
        newdata <- adlp$data
    }

    U <- stats::runif(n)
    all_sims <- c()

    for (i in 1:n) {

        component_sim = calc_adlp_component_lst(adlp$components_lst, newdata, model = "full", calc = "sim")
        sim_partitions <- adlp$partition_func(component_sim)
        n.partitions <- length(sim_partitions)
        ensemble_w <- adlp$model_weights

        sims <- c()
        sim_ij <- c()
        data_ij <- paste(newdata$origin, newdata$dev, sep = "-")

        for (j in 1:n.partitions) {
            mix_ind_subset <- findInterval(U[i], cumsum(unlist(ensemble_w[[j]])))+1
            sim <- sim_partitions[[j]][, -c(1, 2)][,mix_ind_subset]

            sims <- c(sims, sim)
            sim_ij <- c(sim_ij, paste(sim_partitions[[j]]$origin, "-", sim_partitions[[j]]$dev, sep = ""))
        }

        sim_index <- match(data_ij, sim_ij)
        simulation <- sims[sim_index]
        simulation <- cbind(list(sim=i), newdata[, 1:2], simulation)

        all_sims <- rbind(all_sims, simulation)
    }

    all_sims
}
