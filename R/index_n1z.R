#' A subdata of size n1 that mimics the proportion of different levels of delta and Z as in the original data.
##' @title A sub sample of size n1 with similar proportion of different levels of delta and Z as in the original data.
##' @return An index of selected samples of size n1
##' @author Atika Urmi
##' @param delta a vector of censoring indicator in the original data.
##' @param Z a vector of conditional covariate (categorized)
##' @param n1 desired sample size in the sub data
#' @export
index_n1z <- function(delta, n1, z) {
    n <- length(delta)
    ind <- NULL
    ud <- unique(delta)
    uz <- unique(z)
    for (j in 1:length(ud)) {
        for (i in 1:length(uz)) {
            ind1 <- sample(which(delta == ud[j] & z == uz[i]), n1 * (sum(delta == ud[j] & z == uz[i])/n))
            ind <- c(ind, ind1)
        }
    }
    ind
}
