
#' For a data with size n and c% censorig rate, this function returns the index of a sub data of desired size n1 with c% censoring
##' @title Randomly selected sub sample of size n1 for prespecified censoring rate c
##' @return an index of the sub sample of size n1 with c% censoring rate
##' @author Atika Urmi
##' @param delta a vector of censoring indicator in the original data
##' @param n1 desired sample size in the sub data
#' @export

index_n1 <- function(delta, n1) {
    n = length(delta)
    p1 <- round(sum(delta == 1) * (n1/n), 0)  #percent of event in the full data
    p2 <- round(sum(delta == 0) * (n1/n), 0)  #percent of censored in the full data
    ind <- sort(c(sample(which(delta == 1), p1), sample(which(delta == 0), p2)))
    ind
}
