
##' @title Nadara-Watson estimates of smoothing kernel G
##' @return returns a (nxn) matrix of the estimated probability density function (PDF) of a random vector 
##' @author Chenlu Ke
##' @param y the one dimensional vector for which the density is to be estimated
##' @param n length of y
#' @export

NW_RBF <- function(y, n) {
    bw = 1.06 * sd(y) * n^(-1/5)
    dist = dist(y, diag = T, upper = T)^2
    G = exp(-as.matrix(dist)/2/bw^2)
    Gstar = G/rowSums(G)
    tGG = t(Gstar) %*% Gstar
    return(tGG)
}
