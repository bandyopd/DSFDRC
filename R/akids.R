#' a_KIDS-alpha-controlled kernel-based independence dual screening
##' @title a_KIDS
##' @return lists containing covariates selected in the screening step (d1) and final step (d2) for a given FDR level
##' @author Atika Farzana Urmi, Chenlu Ke
##' @param x a matrix of covariates
##' @param y observed time points
##' @param delta the censoring indicator where 1 indicates the event, 0 indicates censored.
##' @param n1 the desired sample size for the screening step
##' Note: n1+n2=n and n1<n2
##' @param d number of predictors to be selected in the screening step
##' @param fdr a single/vector of prespecified type I error rate.
##' @param swap if 'TRUE', interchanges two continuous variable (X and Y) while calculating the smoothing Kernel G().
##' default is 'FALSE'.
#' @export
a_KIDS <- function(x, y, delta, n1, d, fdr, swap = F) {
    z1 <- index_n1(delta, n1)
    step1 <- KIDS(x[z1, ], y[z1], delta[z1], swap = swap)
    d1 <- unname(unlist(step1[2])[1:d])
    # knockoff construction
    xo <- x[, d1]  #original d1 covariates (n1+n2)
    xk <- x[, d1]
    xk[-z1, ] <- create.second_order(as.matrix(xo[-z1, ]), method = "sdp")  # d1covariates(original n1+knockoff n2)
    # step2

    W1 <- KIDS(xo, y, delta, swap = swap)$omega[, 1] - KIDS(xk, y, delta, swap = swap)$omega[, 1]
    W2 <- KIDS(xo, y, delta, swap = swap)$omega[, 2] - KIDS(xk, y, delta, swap = swap)$omega[, 2]
    z2 <- threshold(W1, W2, fdr)
    d2 <- vector("list", length(fdr))
    for (i in 1:nrow(z2)) {
        d2[[i]] <- unique(as.vector(d1[which(W1 >= z2[i, 1] | W2 >= z2[i, 2])]))
        ifelse(length(d2[[i]]) == 0, d2[[i]] <- 0, d2[[i]])
    }
    list(d1, d2)
}





