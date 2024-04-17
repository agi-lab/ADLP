
check_df <- function(df) {

    err_msg <- "Inputted dataframe does not have the
    correctly specified columns, "

    if (!methods::is(df, "data.frame")) {stop(err_msg)}

    if (is.null(df$dev) | is.null(df$origin)) {stop(err_msg)}
}

#' Train-Validation Split of Claims Triangle
#'
#' General framework for any user-defined training/validation/testing split of
#' claims triangle data. The ADLP package contains three default splitting algorithms.
#'
#' @usage train_val_split(df, ...)
#' @param df Claims Triangle and other information. `data.frame` format of
#'  claims and related information for each cell. Dataframe will have columns
#'   `origin` and `dev` as columns 1 and 2 respectively.
#' @param ... Other parameters used to calculate train/test splitting.
#'
#' @return List containing `$train`, `$valid`, `$test`, which should partition
#'  the input `df`.
#'
#' @seealso \link[ADLP]{train_val_split_by_AP},
#'  \link[ADLP]{train_val_split_method1},
#'  \link[ADLP]{train_val_split_method2}
#' @name train_val_split
train_val_split <- function(df, ...) {
    return (
        list(
            train = df,
            valid = df,
            test = df
        )
    )
}

#' Train-Validation Split by Accident Period
#'
#' Function for training/validation splitting.
#'
#' @param df Claims Triangle and other information. `data.frame` format of
#'  claims and related information for each cell. Dataframe will have columns
#'   `origin` and `dev` as columns 1 and 2 respectively.
#' @param accident_periods Vector of accident periods. Will be equivalent to
#'  `1:Triangle_Size`
#' @param max_dev_periods Vector of development periods
#' @param test Returns the test set if `TRUE`
#'
#' @return List containing `$train`, `$valid`, `$test`, which should partition
#'  the input `df`.
#'
#' @details
#' Assigns training set defined by a maximum development period for each
#' accident period: \eqn{(x_{ij} <= MaxDP(i))}.
#'
#' Validation set is therefore cells outside of this period but within the
#' upper triangle. The test set is all observations in the lower triangle.
#'
#' @examples
#' data("test_claims_dataset")
#'
#' train_val <- train_val_split_by_AP(
#'     df = test_claims_dataset,
#'     accident_periods = 1:40,
#'     max_dev_periods = 40:1,
#'     test = TRUE
#' )
#'
#' @seealso \link[ADLP]{train_val_split}
#' @export
train_val_split_by_AP <- function(
        df,
        accident_periods,
        max_dev_periods,
        test = FALSE
) {
    check_df(df)
    stopifnot({
        length(accident_periods) == length(max_dev_periods)
        all(sort(accident_periods) == min(accident_periods):max(accident_periods))
    })

    train_idx <- rep(FALSE, nrow(df))
    for (i in 1:length(accident_periods)) {
        ap <- accident_periods[i]
        max_dp <- max_dev_periods[i]

        train_idx <- train_idx | (df$origin == ap & df$dev <= max_dp)
    }

    val_idx <- (df$origin + df$dev <= (max(accident_periods) + 1)) & !train_idx

    if (test) {
        z <- list(
            train = df[train_idx, ],
            valid = df[val_idx, ],
            test = df[!(train_idx | val_idx), ]
        )
    } else {
        z <- list(
            train = df[train_idx, ],
            valid = df[val_idx, ]
        )
    }

    return (z)
}

#' Train-Validation Split by Accident Period Method 1
#'
#' Function for training/validation splitting.
#'
#' @param df Claims Triangle and other information. `data.frame` format of
#'  claims and related information for each cell. Dataframe will have columns
#'   `origin` and `dev` as columns 1 and 2 respectively.
#' @param tri.size Triangle size.
#' @param val_ratio Value between 0 and 1 as the approximate size of validation
#'  set.
#' @param test Returns the test set if `TRUE` .
#'
#' @return List containing `$train`, `$valid`, `$test`, which should partition
#'  the input `df`.
#'
#' @details
#' Approximates the validation set by taking the `n` most recent calendar years
#' as validation to best fit `val_ratio`.
#'
#' Validation set is therefore cells outside of this period but within the
#' upper triangle. The test set is all observations in the lower triangle.
#'
#' Note that accident period 1 and development period 1 will always be within
#' the training set.
#'
#' @examples
#'
#' data("test_claims_dataset")
#'
#' train_val <- train_val_split_method1(
#'     df = test_claims_dataset,
#'     tri.size = 40,
#'     val_ratio = 0.3,
#'     test = TRUE
#' )
#'
#' @seealso \link[ADLP]{train_val_split}
#' @export
train_val_split_method1 <- function(df, tri.size, val_ratio, test = FALSE) {

    check_df(df)
    stopifnot({
        tri.size >= max(df$origin, df$dev);
        val_ratio >= 0;
        val_ratio <= 1;
    })

    # Number of validation cells per accident period (approx)
    # = val_ratio * tri.size * (tri.size + 1) / 2 /tri.size  = val_ratio * (tri.size + 1) / 2
    # Note: the total number of data points of the upper triangle is: tri.size * (tri.size + 1) / 2
    n = floor(val_ratio * (tri.size + 1) / 2)
    train_AP <- 1:tri.size
    train_DP <- rev(train_AP) - 1 - n
    # Keep all observations from AP1 and DP1
    train_DP[1] = tri.size
    train_DP = pmax(train_DP, 1)

    z <- train_val_split_by_AP(df, train_AP, train_DP, test)
    return (z)
}

#' Train-Validation Split by Accident Period Method 2
#'
#' Function for training/validation splitting.
#'
#' @param df Claims Triangle and other information. `data.frame` format of
#'  claims and related information for each cell. Dataframe will have columns
#'   `origin` and `dev` as columns 1 and 2 respectively.
#' @param tri.size Triangle size.
#' @param val_ratio Value between 0 and 1 as the approximate size of validaiton
#'  set.
#' @param test Returns the test set if `TRUE` .
#'
#' @return List containing `$train`, `$valid`, `$test`, which should partition
#'  the input `df`.
#'
#' @details
#' Approximates the validation set by defining the training set as the cells
#' below the function \eqn{((b^{1/a} - x^{1/a})^a)}. Where \eqn{b} is equal to
#' the triangle size and \eqn{a} is optimised to best fit `val_ratio`.
#'
#' The training set is therefore cells outside of this period but within the
#' upper triangle. The test set is all observations in the lower triangle.
#'
#' Note that accident period 1 and development period 1 will always be within
#' the training set.
#'
#' @examples
#'
#' data("test_claims_dataset")
#'
#' train_val <- train_val_split_method1(
#'     df = test_claims_dataset,
#'     tri.size = 40,
#'     val_ratio = 0.3,
#'     test = TRUE
#' )
#'
#' @seealso \link[ADLP]{train_val_split}
#' @export
train_val_split_method2 <- function(df, tri.size, val_ratio, test = FALSE) {

    # Checks
    check_df(df)
    stopifnot({
        tri.size >= max(df$origin, df$dev);
        val_ratio >= 0;
        val_ratio <= 1;
    })

    # Defines Splitting Fucntion
    split_fun <- function(x, a) {
        (tri.size^(1/a) - x^(1/a))^a
    }

    alpha <- stats::uniroot(
        function(a) a*beta(a, a+1)*tri.size^2 / (tri.size*tri.size/2) - (1 - val_ratio),
        interval = c(1, 5)
    )$root

    train_AP <- 1:tri.size
    train_DP <- round(split_fun(train_AP-1, alpha))
    # Keep all observations from AP1 and DP1
    train_DP[1] = tri.size
    train_DP = pmax(train_DP, 1)

    z <- train_val_split_by_AP(df, train_AP, train_DP, test)
    return (z)
}



