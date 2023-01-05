#' This function calculates the ECCFIC based correlation between two random vector X and Y.
#' X is always continious here, Y can be either continuous or discrete
##' @title ECCFIC- correlation
##' @return ECCFIC correlation between X and Y
##' @author Chenlu Ke
##' @param x a numeric vector/matrix
##' @param y a numeric vector
##' @return ECCFIC correlation between X and Y
#' @param est 1 for continuous y, 2 for categorical y
#' @param kernel a kernel to use for x, 'gaussian' or 'distance'
#' @param alpha exponent on distance in (0,2], if a 'distance' kernel is used
#' @param sigma bandwidth, if a 'gaussian' kernel is used,
#'              default is heuristic median pairwise distances of x
#' @param bw bandwidth of the gaussian smoothing kernel applied on y,
#'           bandwidths suggested by Silverman (1986) are used unless otherwise specified.
ECH2 <- function(x, y, est = 1, kernel = "gaussian", sigma = "default", alpha = 1, bw = "default") {

    x = as.matrix(x)
    n = nrow(x)
    if (n == 1) {
        return(0)
    }

    if (kernel == "gaussian") {
        if (sigma == "default") {
            sigma = sqrt(0.5 * median(dist(x)^2))
            if (sigma == 0) {
                sigma = 0.001
            }
        }
    }

    if (identical(x, as.matrix(y))) {
        if (kernel == "gaussian") {
            K = dnorm(as.matrix(dist(x, diag = T, upper = T)), mean = 0, sd = sigma)
            return((dnorm(0, 0, sigma) - mean(K)) * sqrt(2 * pi) * sigma)
        } else {
            if (kernel == "distance") {
                K = 0.5 * (as.matrix(dist(x, diag = T, upper = T))^alpha)
                return(mean(K))
            }
        }
    }

    if (est == 1) {
        y = as.matrix(y)
        p = ncol(y)
        if (bw == "default") {
            if (p == 1) {
                bw = 1.06 * sd(y) * n^(-1/5)
            } else {
                bw = (4/n/(p + 2))^(1/(p + 4)) * sum(diag(cov(y)))/p
            }
        }
        H = diag(n) - matrix(1, n, n)/n
        if (kernel == "gaussian") {
            K = dnorm(as.matrix(dist(x, diag = T, upper = T)), mean = 0, sd = sigma) * sqrt(2 * pi) * sigma
        } else {
            if (kernel == "distance") {
                normx = matrix(apply(x * x, 1, sum)^(alpha/2), n, n)
                K = 0.5 * (normx + t(normx) - as.matrix(dist(x, diag = T, upper = T))^alpha)
            }
        }
        G = dnorm(as.matrix(dist(y, diag = T, upper = T)), mean = 0, sd = bw)
        Gstar = t(G/rowSums(G))
        return(sum(diag(K %*% H %*% Gstar %*% t(Gstar) %*% H))/n)
    } else {
        L = ifelse(as.matrix(dist(y, diag = T, upper = T)) == 0, 1, 0)
        L = n * L/rowSums(L) - 1
        if (kernel == "gaussian") {
            K = dnorm(as.matrix(dist(x, diag = T, upper = T)), mean = 0, sd = sigma)
            return((sum(K * L)/(n^2)) * sqrt(2 * pi) * sigma)
        } else {
            if (kernel == "distance") {
                K = -0.5 * (as.matrix(dist(x, diag = T, upper = T))^alpha)
                return(sum(K * L)/(n^2))
            }
        }
    }

}

rho_ech2 <- function(x, y, est = 1, ...) {
    return(ECH2(x, y, est, ...)/ECH2(x, x, est, ...))
}

