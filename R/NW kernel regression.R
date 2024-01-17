
##' @title Nadaraya-Watson kernel regression
##' @return returns a (nxn) matrix of kernel regression estimate using Nadaraya-Watson kernel regression
##' @author Chenlu Ke
##' @param y the vector to calculate the kernel
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
