
#' Accident and Development period Adjusted Linear Pools Component Models
#'
#' Class to store component models and related functions required for ADLP
#' estimation, prediction and goodness of fit.
#'
#' @param model_train Model trained on training data
#' @param model_full Model trained on all in-sample data
#' @param calc_dens function to calculate the pdf of each point
#' @param calc_mu function to calculate the estimated mean of each point
#' @param calc_cdf function to calculate the cdf of each point
#' @param sim_fun function to simulate new from
#' @param ... Other named parameters required for the model or any of its
#'  related functions to run.
#'
#' @return Object of class `adlp_component`
#'
#' @details
#' Component models `model_train` and `model_full` are designed to be objects of
#' class `glm`, `lm`, or similar. The models would desirably have a S3 method for
#' `formula. For models that do not fit under this umbrella,
#' see \link[custom_model]{custom_model}. For a potential list of candidate models,
#' one might refer to Avanzi, Li, Wong and Xian (2022).
#'
#' Functions as assumed to have the following parameter naming convention:
#' \itemize{
#'  \item{`y` as the response variable}
#'  \item{`model` as the modeling object `model_train` or `model_full`}
#'  \item{`newdata` to designate new data}
#' }
#' Other inputs not in this list will need to be intialised with the `adlp_component`
#'
#' @references Avanzi, B., Li, Y., Wong, B., & Xian, A. (2022). Ensemble distributional forecasting for insurance loss reserving. arXiv preprint arXiv:2206.08541.
#' @export
adlp_component <- function(
        model_train,
        model_full,
        calc_dens,
        calc_mu,
        calc_cdf,
        sim_fun,
        ...
) {
    z = list(
        model_train = model_train,
        model_full = model_full,
        calc_dens = calc_dens,
        calc_mu = calc_mu,
        calc_cdf = calc_cdf,
        sim_fun = sim_fun,
        ...
    )

    class(z) <- "adlp_component"
    return (z)
}

#' Accident and Development period Adjusted Linear Pools Component Models
#'
#' @param ... Individual `adlp_components`
#'
#' @details
#' Class to structure a list of \link[adlp_component]{adlp_components}.
#'
#' @name adlp_components
#' @export
adlp_components <- function(...) {
    z <- list(...)

    for (i in names(z)) {
        if (class(z[[i]]) != "adlp_component") {
            stop(paste("Argument `", i, "` is not an adlp_component", collapse = ""))
        }
    }

    class(z) <- "adlp_components"
    return (z)
}

#' @name component_extract_model
#' Extracts the desired model from component depending on input
#' @noRd
component_extract_model <- function(component, model = c("train", "full")) {
    model <- match.arg(model)
    if (model == "train") {
        component_model = component$model_train
    } else if (model == "full") {
        component_model = component$model_full
    } else {
        stop ("Incorrect argument `model`.")
    }
    return (component_model)
}

#' Accident and Development period Adjusted Linear Pools Component Models
#'
#' @param component Object of class `adlp_component`
#' @param newdata Claims Triangle and other information. `data.frame` format of
#'  claims and related information for each cell. Dataframe will have columns
#'   `origin` and `dev` as columns 1 and 2 respectively.
#' @param y Optional vector of `y` to be used in pdf or cdf calculations. Will
#' default to the response fetched by `model.frame`.
#' @param model Whether the training component model or the full component model
#' should be used
#' @param calc Type of calculation to perform
#'
#' @details
#' Calls the specified function for an object of class `adlp_component`.
#'
#' @export
calc_adlp_component <- function(component, newdata,
                                y = NULL,
                                model = c("train", "full"),
                                calc = c("pdf", "cdf", "mu", "sim")) {

    calc <-  match.arg(calc)
    model <- match.arg(model)
    component_model <- component_extract_model(component, model)

    if (is.null(y)) {
        y <- stats::model.frame(stats::formula(component_model), newdata)[, 1]
    }

    if (calc == "pdf") {
        func <- component$calc_dens
    } else if (calc == "cdf") {
        func <- component$calc_cdf
    } else if (calc == "mu") {
        func <- component$calc_mu
    } else if (calc == "sim") {
        func <- component$sim_fun
    } else {
        stop ("Incorrect argument `calc`.")
    }

    dens_args <- names(formals(func))
    args <- list()
    if ("y" %in% dens_args) {args <- c(args, list(y=y))}
    if ("model" %in% dens_args) {args <- c(args, list(model=component_model))}
    if ("newdata" %in% dens_args) {args <- c(args, list(newdata=newdata))}

    # Only select necessary/required arguments from component to be used in dens
    if (any(names(component) %in% dens_args)) {
        args <- c(
            args,
            component[names(component)[names(component) %in% dens_args]]
        )
    }

    return (
        do.call(
            what = func,
            args = args
        )
    )
}

#' @rdname calc_adlp_component
#'
#' @details
#' `calc_adlp_component_lst` is a wrapper for `calc_adlp_component` for each
#' component in the list `components_lst`. This wrapper also contains functionality
#' to signal the component that causes an error if it is occuring downstream.
#'
#' @export
calc_adlp_component_lst <- function(
        components_lst, newdata,
        model = c("train", "full"),
        calc = c("pdf", "cdf", "mu", "sim"), ...) {

    model = match.arg(model)
    calc = match.arg(calc)
    component_dens <- do.call(
        "cbind",
        lapply(names(components_lst), function(x) {
            tryCatch(
                expr = {
                    calc_adlp_component(
                        component = components_lst[[x]],
                        newdata = newdata,
                        model = model,
                        calc = calc,
                        ...
                    )
                },
                error = function(err) {
                    err$message <- paste("Error occured in component function:", x, "\n", err)
                    stop(err)
                }
            )
        })
    )
    colnames(component_dens) <- names(components_lst)
    if (nrow(component_dens) >= nrow(newdata)) {
        origin <- rep_len(newdata[,1], nrow(component_dens))
        dev <- rep_len(newdata[,2], nrow(component_dens))
        dens <- cbind(origin, dev, component_dens)
    } else {
        dens <- cbind(newdata[,1:2], component_dens)
    }
    as.data.frame(dens)
}

adlp_component_lst_logS <- function(components_lst, newdata, model = c("train", "full"), epsilon = 1e-8) {
    model <- match.arg(model)
    component_dens <- calc_adlp_component_lst(components_lst, newdata, model, "pdf")
    dens_index <- component_dens[, 1:2]
    dens_val <- component_dens[, -(1:2)]
    dens_val <- log(dens_val + epsilon)
    z <- cbind(dens_index, dens_val)
    z
}
