#' These kernel functions returnf Gaussian kernel as a reproducing (K()) or smoothing (G()) kernel
##' @title Reproducing Gaussian Kernel K(.)
##' @return Reproducing kernel (K) for a vector
##' @author Chenlu Ke
##' @param x the vector to calculate the kernel
##' @param sigma2 the bandwidth of the reproducing kernel K(), default is the heuristic median pairwise distance
#' @export

K <- function(x, sigma2 = "default") {
    dist = dist(x, diag = T, upper = T)^2
    if (sigma2 == "default") {
        sigma2 = 0.5 * median(dist)
    }
    if (sigma2 == 0) {
        sigma2 = 0.001
    }
    return(exp(-as.matrix(dist)/2/sigma2))
}

##' @title Smoothing Kernel G(.)
##' @return Smoothing Kernel (G) for a vector
##' @author Chenlu Ke
##' @param x the vector to calculate the kernel
##' @param n length of x
#' @export

G <- function(y, n) {
    bw = 1.06 * sd(y) * n^(-1/5)
    dist = dist(y, diag = T, upper = T)^2
    A = exp(-as.matrix(dist)/2/bw^2)
    Astar = A/rowSums(A)
    tAA = t(Astar) %*% Astar
    return(tAA)
}
