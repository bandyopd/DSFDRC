
#'  A subdata of size n1 that mimics the proportion of censored observations in the original data.
##' @title A sub sample of size n1 with censoring rate similar to original data
##' @return An index of selected samples of size n1.
##' @author Atika Urmi
##' @param delta a vector of censoring indicator in the original data. This is required to calculate the censoring rate in original data.
##' @param n1 desired sample size in the sub data
#' @export

index_n1 <- function(delta, n1) {
    n = length(delta)
    p1 <- round(sum(delta == 1) * (n1/n), 0)  #percent of event in the full data
    p2 <- round(sum(delta == 0) * (n1/n), 0)  #percent of censored in the full data
    ind <- sort(c(sample(which(delta == 1), p1), sample(which(delta == 0), p2)))
    ind
}
