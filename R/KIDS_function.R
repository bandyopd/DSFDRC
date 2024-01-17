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

    n <- nrow(x)
    p <- ncol(x)
    id0 <- delta == 0
    n0 <- sum(id0)
    n1 <- n - n0
    w0 <- n0/n
    w1 <- 1 - w0
    ohat <- matrix(0, p, 2)

    if (swap == F) {
        x0 <- x[id0, ]
        x1 <- x[!id0, ]
        # bwy <- 1.06*sd(y)*n^(-1/5)
        tGGy0 <- NW_RBF(y[id0], n0)
        tGGy1 <- NW_RBF(y[!id0], n1)

        Omega <- function(u, u0, u1) {

            distu = dist(u, diag = T, upper = T)^2
            sigma2u = 0.5 * median(distu)
            if (sigma2u == 0) {
                sigma2u = 0.001
            }
            K = exp(-as.matrix(distu)/2/sigma2u)
            K0 = exp(-as.matrix(dist(u0, diag = T, upper = T)^2)/2/sigma2u)
            K1 = exp(-as.matrix(dist(u1, diag = T, upper = T)^2)/2/sigma2u)
            mK = mean(K)
            mK0 = mean(K0)
            mK1 = mean(K1)
            wmK = mK0 * w0 + mK1 * w1
            omega1 = (wmK - mK)/(1 - mK)

            Hy0 = sum(tGGy0 * K0)/n0 - mK0
            Hy1 = sum(tGGy1 * K1)/n1 - mK1
            omega2 = (Hy0 * w0 + Hy1 * w1)/(1 - wmK)

            return(c(omega1, omega2))
        }

        for (j in 1:p) {
            ohat[j, ] = Omega(x[, j], x0[, j], x1[, j])
        }

    } else {
        x0 <- x[id0, ]
        x1 <- x[!id0, ]
        disty = dist(y, diag = T, upper = T)^2
        sigma2y = 0.5 * median(disty)
        Ky0 <- exp(-as.matrix(dist(y[id0], diag = T, upper = T)^2)/2/sigma2y)
        Ky1 <- exp(-as.matrix(dist(y[!id0], diag = T, upper = T)^2)/2/sigma2y)
        mKy0 = mean(Ky0)
        mKy1 = mean(Ky1)
        wmKy = mKy0 * w0 + mKy1 * w1

        Omega <- function(u, u0, u1) {

            distu = dist(u, diag = T, upper = T)^2
            sigma2u = 0.5 * median(distu)
            if (sigma2u == 0) {
                sigma2u = 0.001
            }
            K = exp(-as.matrix(distu)/2/sigma2u)
            K0 = exp(-as.matrix(dist(u0, diag = T, upper = T)^2)/2/sigma2u)
            K1 = exp(-as.matrix(dist(u1, diag = T, upper = T)^2)/2/sigma2u)
            mK = mean(K)
            mK0 = mean(K0)
            mK1 = mean(K1)
            wmK = mK0 * w0 + mK1 * w1
            omega1 = (wmK - mK)/(1 - mK)

            # bwx <- 1.06*sd(u)*n^(-1/5)
            tGGx0 = NW_RBF(u0, n0)
            tGGx1 = NW_RBF(u1, n1)
            Hy0 = sum(tGGx0 * Ky0)/n0 - mKy0
            Hy1 = sum(tGGx1 * Ky1)/n1 - mKy1
            omega2 = (Hy0 * w0 + Hy1 * w1)/(1 - wmKy)

            return(c(omega1, omega2))
        }

        for (j in 1:p) {
            ohat[j, ] = Omega(x[, j], x0[, j], x1[, j])
        }

    }

    r1 <- rank(-ohat[, 1])
    r2 <- rank(-ohat[, 2])
    rank_index <- order(pmin(r1, r2), pmax(r1, r2))
    colnames(ohat) <- paste("ohat", 1:2, sep = "")
    list(omega = ohat, top_covariates = rank_index)


}

