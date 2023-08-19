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
            Kx0 <- K(x0, sigma2x)
            Kx1 <- K(x1, sigma2x)

            mK = mean(K(x[, i], sigma2x))
            mK0 = mean(Kx0)
            mK1 = mean(Kx1)
            wmK = w0 * mK0 + w1 * mK1
            omega1 <- (wmK - mK)/(1 - mK)

            Gy0 <- G(y0, n0)
            Gy1 <- G(y1, n1)
            a <- sum(Gy0 * Kx0)/n0 - mK0
            b <- sum(Gy1 * Kx1)/n1 - mK1
            omega21 <- (w0 * a + w1 * b)/(1 - wmK)

            omega[i, ] <- c(omega1, omega21)
        }

    } else {

        for (i in 1:p) {
            x0 <- x[, i][delta == 0]
            x1 <- x[, i][delta == 1]
            sigma2y = 0.5 * median(dist(y, diag = T, upper = T)^2)
            Ky0 <- K(y0, sigma2y)
            Ky1 <- K(y1, sigma2y)
            mKy0 = mean(Ky0)
            mKy1 = mean(Ky1)
            wmKy = mKy0 * w0 + mKy1 * w1

            sigma2x = 0.5 * median(dist(x[, i], diag = T, upper = T)^2)
            Kx = K(x[, i], sigma2x)
            Kx0 = K(x0, sigma2x)
            Kx1 = K(x1, sigma2x)
            mKx = mean(Kx)
            mKx0 = mean(Kx0)
            mKx1 = mean(Kx1)
            wmKx = mKx0 * w0 + mKx1 * w1
            omega1 = (wmKx - mKx)/(1 - mKx)

            Gx0 = G(x0, n0)
            Gx1 = G(x1, n0)
            a = sum(Gx0 * Ky0)/n0 - mKy0
            b = sum(Gx1 * Ky1)/n1 - mKy1
            omega22 = (a * w0 + b * w1)/(1 - wmKy)

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

