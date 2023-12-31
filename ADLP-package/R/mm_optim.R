
#' @rdname MM_optim
#' @noRd
MM_func <- function(w,dat){
    wd <- dat%*%w
    dens <- dat / as.vector(wd)
    m <- apply(dens, MARGIN=2, FUN=mean)
    return(m)
}

#' Minorization-Maximisation Algorithm performed to fit the ADLPs
#' @description The Minorization-Maximization algorithm aims to optimize a surrogate objective function that approximates the Log Score. This approach typically results in fast and stable convergence, while ensuring that combination weights adhere to the constraints of being non-negative and summing to one. For detailed description of the algorithm, one might refer to: Conflitti, De Mol, and Giannone (2015)
#' @references Conflitti, Cristina, Christine De Mol, and Domenico Giannone. "Optimal combination of survey forecasts." International Journal of Forecasting 31.4 (2015): 1096-1103.
#' @name MM_optim
#' @export
MM_optim <- function(w_init, dat, niter = 500, ...){
    # To avoid taking logs and division of zero densities
    #dat <- dat + 1e-16
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


