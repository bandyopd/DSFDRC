#' KIDS-kernel-based independence dual screening
##' @title KIDS
##' @return a list of a matrix of correlation(omega) and the order of covariates based on the highest correlations.
##' @author Chenlu Ke, Atika Urmi
##' @param x a matrix of covariates
##' @param y observed time points
##' @param delta the censoring indicator where 1 indicates the event, 0 indicates censored.
##' @param swap if 'TRUE', interchanges two continuous variable (X and Y) while calculating the smoothing Kernel G().
##' default is 'FALSE'.
#' @export
KIDS <- function(x, y, delta, swap = F) {
    # x=covariate matrix/data

    n <- nrow(x)
    n0 <- sum(delta == 0)
    n1 <- n - n0
    w0 <- n0/n  #weights for delta==0
    w1 <- n1/n
    p = ncol(x)
    y0 <- y[delta == 0]
    y1 <- y[delta == 1]
    omega <- matrix(0, p, 2)

    if (swap == F) {
        # ir=T

        for (i in 1:p) {
            x0 <- x[, i][delta == 0]
            x1 <- x[, i][delta == 1]
            sigma2x = 0.5 * median(dist(x[, i], diag = T, upper = T)^2)

            mK = mean(K(x[, i], sigma2x))
            mK0 = mean(K(x0, sigma2x))
            mK1 = mean(K(x1, sigma2x))
            wmK = w0 * mK0 + w1 * mK1
            omega1 <- (wmK - mK)/(1 - mK)


            a <- sum(G(y0, n0) * K(x0, sigma2x))/n0 - mK0
            b <- sum(G(y1, n1) * K(x1, sigma2x))/n1 - mK1
            c <- w0 * mK0 + w1 * mK1
            omega21 <- (w0 * a + w1 * b)/(1 - c)

            omega[i, ] <- c(omega1, omega21)
        }

    } else {

        for (i in 1:p) {
            x0 <- x[, i][delta == 0]
            x1 <- x[, i][delta == 1]
            sigma2x = 0.5 * median(dist(x[, i], diag = T, upper = T)^2)
            mK = mean(K(x[, i], sigma2x))
            mK0 = mean(K(x0, sigma2x))
            mK1 = mean(K(x1, sigma2x))
            wmK = w0 * mK0 + w1 * mK1
            omega1 <- (wmK - mK)/(1 - mK)

            sigma2y = 0.5 * median(dist(y, diag = T, upper = T)^2)
            mKy0 <- mean(K(y0, sigma2y))
            mKy1 <- mean(K(y1, sigma2y))
            mKy <- mean(K(y, sigma2y))
            a <- sum(G(x0, n0) * K(y0, sigma2y))/n0 - mKy0
            b <- sum(G(x1, n1) * K(y1, sigma2y))/n1 - mKy1
            c <- w0 * mKy0 + w1 * mKy1
            omega22 <- (w0 * a + w1 * b)/(1 - c)

            omega[i, ] <- c(omega1, omega22)
        }
    }

    omega
    r1 <- rank(-omega[, 1])
    r2 <- rank(-omega[, 2])
    rank_index <- order(pmin(r1, r2), pmax(r1, r2))  #returns the index of x for lowest to highest rank
    # if returns 4, then x4 has rank1
    colnames(omega) <- paste("omega", 1:2, sep = "")
    list(omega = omega, top_covariates = rank_index)
}
