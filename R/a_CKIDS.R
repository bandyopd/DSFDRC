
#' aCKIDS-alpha controlled Conditional kernel-based independence dual screening
##' @title a_CKIDS
##' @return lists containing covariates selected in the screening step (d1) and final step (d2) conditioning on Z with a given FDR level
##' @author Atika Urmi
##' @param x a matrix of covariates
##' @param y observed time points
##' @param z conditional covariate (categorized)
##' @param delta the censoring indicator where 1 indicates the event, 0 indicates censored.
##' @param n1 the desired sample size for the screening step
##' Note: n1+n2=n and n1<n2
##' @param d number of predictors to be selected in the screening step
##' @param fdr a single/vector of prespecified type I error rate.
##' @param swap  indicates whether inverse regression method should be considered. If swap='TRUE',it will estimate Y given X. otherwise, it will estimate X|Y.
#' @export

a_CKIDS <- function(x, y, delta, n1, d, z, fdr, swap) {
    # z=binary factorized conditional variable
    x <- data.matrix(x)
    z1 <- index_n1z(delta, n1, z)
    step1 <- CKIDS(x[z1, ], y[z1], delta[z1], z[z1], swap)
    d1 <- unname(unlist(step1[2])[1:d])
    # knockoff construction
    xo <- x[, d1]  #original d1 covariates (n1+n2)
    xk <- x[, d1]
    catz <- unique(z)
    for (i in 1:length(catz)) {
        tind <- setdiff(which(z == catz[i]), z1)
        xk[tind, ] <- create.second_order(xo[tind, ], method = "sdp")
    }

    # step2

    W1 <- CKIDS(xo, y, delta, z, swap)$omega[, 1] - CKIDS(xk, y, delta, z, swap)$omega[, 1]
    W2 <- CKIDS(xo, y, delta, z, swap)$omega[, 2] - CKIDS(xk, y, delta, z, swap)$omega[, 2]
    z2 <- threshold(W1, W2, fdr)
    d2 <- vector("list", length(fdr))
    for (i in 1:nrow(z2)) {
        d2[[i]] <- unique(as.vector(d1[which(W1 >= z2[i, 1] | W2 >= z2[i, 2])]))
        ifelse(length(d2[[i]]) == 0, d2[[i]] <- 0, d2[[i]])
    }
    list(d1, d2)
}

## a_KIDS
a_KIDS <- function(x, y, delta, n1, d, fdr, swap) {
    z1 <- index_n1(delta, n1)
    # z1<-sample(1:length(delta),n1,F)
    step1 <- KIDS(x[z1, ], y[z1], delta[z1], swap)
    d1 <- unname(unlist(step1[2])[1:d])
    # knockoff construction
    xo <- x[, d1]  #original d1 covariates (n1+n2)
    xk <- x[, d1]
    xk[-z1, ] <- create.second_order(as.matrix(xo[-z1, ]), method = "sdp")  # d1covariates(original n1+knockoff n2)
    # step2

    W1 <- KIDS(xo, y, delta, swap)$omega[, 1] - KIDS(xk, y, delta, swap)$omega[, 1]
    W2 <- KIDS(xo, y, delta, swap)$omega[, 2] - KIDS(xk, y, delta, swap)$omega[, 2]
    z2 <- threshold(W1, W2, fdr)
    d2 <- vector("list", length(fdr))
    for (i in 1:nrow(z2)) {
        d2[[i]] <- unique(as.vector(d1[which(W1 >= z2[i, 1] | W2 >= z2[i, 2])]))
        ifelse(length(d2[[i]]) == 0, d2[[i]] <- 0, d2[[i]])
    }
    list(d1, d2)
}





