#' CKIDS-Conditional kernel-based independence dual screening
##' @title CKIDS
##' @return a list of a matrix of correlation(omega) and the index of covariates based on the highest correlations with outcome conditioning on Z.
##' @author Atika Urmi
##' @param x a matrix of covariates
##' @param y observed time points
##' @param z conditional covariate (categorized)
##' @param delta the censoring indicator where 1 indicates the event, 0 indicates censored.
##' @param swap  indicates whether inverse regression method should be considered. If swap='TRUE',it will estimate Y given X. otherwise, it will estimate X|Y.
#' @export


CKIDS <- function(x, y, delta, z, swap) {
    # complete data with z being the factorized conditional variable variable
    z <- as.integer(factor(z))
    if (min(z) == 0) {
        z <- z + 1
    }
    total <- nrow(x)
    omega1 <- matrix(0, ncol(x), max(z))
    omega2 <- matrix(0, ncol(x), max(z))
    wz <- NULL
    for (i in 1:max(z)) {
        wz[i] <- sum(z == i)/total
    }
    for (i in 1:max(z)) {
        x1 <- x[z == i, ]
        y1 <- y[z == i]
        delta1 <- delta[z == i]
        dual_screen <- KIDS(x1, y1, delta1, swap)
        omega1[, i] <- wz[i] * dual_screen$omega[, 1]
        omega2[, i] <- wz[i] * dual_screen$omega[, 2]
    }
    omega1con <- rowSums(omega1)
    omega2con <- rowSums(omega2)
    r1 <- rank(-omega1con)
    r2 <- rank(-omega2con)
    rank_index <- order(pmin(r1, r2), pmax(r1, r2))  #returns the index of x for lowest to highest rank
    # if returns 4, then x4 has rank1 colnames(omega) <- paste('omega', 1:2, sep = '')
    list(omega = cbind(omega1con, omega2con), top_covariates = rank_index)
}
