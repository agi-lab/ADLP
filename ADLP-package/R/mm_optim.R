
#' @rdname MM_optim
#' @noRd
MM_func <- function(w,dat){
    wd <- dat%*%w
    dens <- dat / as.vector(wd)
    m <- apply(dens, MARGIN=2, FUN=mean)
    return(m)
}

#' Minorization-Maximisation Algorithm performed to fit the ADLPs
#' @description The Minorization-Maximization algorithm aims to optimize a
#' surrogate objective function that approximates the Log Score. This approach
#'  typically results in fast and stable convergence, while ensuring that
#'   combination weights adhere to the constraints of being non-negative and
#'   summing to one. For detailed description of the algorithm, one might refer
#'    to: Conflitti, De Mol, and Giannone (2015)
#'
#' @param w_init initial weights for each ADLP
#' @param dat matrix of densities for each ADLP
#' @param niter maximum number of iterations. Defaults to 500
#'
#' @references Conflitti, Cristina, Christine De Mol, and Domenico Giannone. "Optimal combination of survey forecasts." International Journal of Forecasting 31.4 (2015): 1096-1103.
##' @name MM_optim
##'
#' @return An object of class `mm_optim`. `mm_optim` is a list that stores the
#' results of the MM algorithm performed, including the final parameters, the
#' final loss and numer of iterations.
#'
#' @examples
#' w_init <- rep(1/3, 3)
#' set.seed(1)
#' density_data <- matrix(runif(9), nrow = 3, ncol = 3)
#' MM_optim(w_init, density_data, niter = 500)
#'
#' @export
MM_optim <- function(w_init, dat, niter = 500){
    # To avoid taking logs and division of zero densities
    dat <- dat + 1e-16
    w <- w_init
    i <- 0
    prev_loss <- -mean(log(dat%*%w_init))
    while (i < niter) {
        w <- w*MM_func(w, dat)
        current_loss <- -mean(log(dat%*%w))
        if (abs(current_loss - prev_loss) < 1e-16) {
            break
        }
        prev_loss <- current_loss
        i <- i + 1
    }

    z <- list(
        finalparams = w,
        finalNLL = current_loss,
        n_iters = i
    )
    class(z) <- "mm_optim"
    return(z)
}


