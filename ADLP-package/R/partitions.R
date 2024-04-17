
#' Accident and Development period Adjusted Linear Pools partition function
#'
#' General framework for any user-defined function for partitions of
#' claims triangle data.
#'
#' @usage adlp_partition(df, ...)
#' @param df data.frame format of claims and related information for each cell.
#'  Dataframe will have columns `origin` and `dev` as columns 1 and 2 respectively.
#' @param ... Other parameters used to calculate ADLP partitions
#'
#' @return List containing the `df` as a result of the partitions.
#'
#' @seealso \link[ADLP]{adlp_partition_none},
#'  \link[ADLP]{adlp_partition_ap}
#' @name adlp_partition
adlp_partition <- function(df, ...) {
    return (adlp_partition_none(df))
}

#' @rdname adlp_partition
#'
#' @details
#' `adlp_partition_none` is the default functionality with no partitions. This is
#' equivalent to the standard linear pooling.
#'
#' @examples
#' data("test_claims_dataset")
#' adlp_partition_none(test_claims_dataset)
#'
#' @export
adlp_partition_none <- function(df) {
    return (list(
        df
    ))
}

#' @rdname adlp_partition
#'
#' @param tri.size Triangle size in claims
#' @param size Number of partitions
#' @param weights a vector of weights for the size of each partition.
#'
#' @details
#' `adlp_partition_ap` will partition the claims triangle by accident period,
#' Where the choice of accident period to partition will be determined to most
#' closely resemble the desired `weights`.
#'
#' The choice of accident period relies on a greedy algorithm that aims to find the
#' accident period that provides the amount of cells that is larger or equal to the
#' desired split.
#'
#' @examples
#' data("test_claims_dataset")
#' adlp_partition_ap(test_claims_dataset, tri.size = 40, size = 3)
#'
#' @export
adlp_partition_ap <- function(df, tri.size, size = 1, weights = rep(1, size)) {
    if (size == 1) {
        return (adlp_partition_none(df))
    }

    weights <- weights/sum(weights)
    partition_sizes <- weights * tri.size * (tri.size + 1) / 2

    # Note that below is a greedy algorithm (i.e. will try to find the closest)
    #   fitting triangle from the latest accident years and moves upwards
    accident_periods <- c(-1)
    for (i in size:2) {
        lower_triangle <- choose(accident_periods[size + 1 - i] + 1, 2)
        k <- floor(0.5*(1 + sqrt(8*(partition_sizes[i] + lower_triangle) + 1)))
        accident_periods <- c(accident_periods, k)
    }

    ap_partitions <- (rev(tri.size - accident_periods[2:size]))
    ap_partitions <- c(0, ap_partitions, tri.size)

    partitions <- list()
    for (i in 1:(size)) {
        partitions[[i]] <- df[
            (df$origin > ap_partitions[i]) & (df$origin <= ap_partitions[i+1]),
        ]
    }
    return (partitions)
}
