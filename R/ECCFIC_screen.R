#' ECCFIC_screen ranks the predictors based on their ECCFIC-correlation with survival time and censoring indicator
##' @title ECCFIC-screen
##' @return ECCFIC correlation between X and Y
##' @author Atika Farzana Urmi
##' @param x_mat a matrix/dataframe of continious covariates
##' @param time a numeric vector of survival time
##' @param delta a numeric vector of censoring indicator
##' @return the corrleation between X and time given delta,correlation between x and delta and the rank of the
##' covariates based on importance
#' @param kernel a kernel to use for x, 'gaussian' or 'distance',default gaussian.
#' if for 'gaussian' kernel, bandwidth is defined as the heuristic median pairwise distances of x
#' @export
ECCFIC_screen <- function(time, delta, x_mat, kernel = "gaussian") {
    x_mat <- data.frame(x_mat)
    data <- cbind(time, delta, x_mat)
    n = nrow(data)
    del0 <- data[which(data$delta == 0), c(3:ncol(data))]
    del1 <- data[which(data$delta == 1), c(3:ncol(data))]
    time0 <- data[which(data$delta == 0), c(1)]
    time1 <- data[which(data$delta == 1), c(1)]
    whole <- data[, c(3:ncol(data))]
    del_whole <- data[, c(2)]

    ECCFIC <- foreach(i = 1:ncol(del0), .combine = c) %do% {
        x <- del0[, c(i)]
        y <- time0
        rho_ech2(x, y, est = 1)
    }
    ECCFIC <- (nrow(del0)/n) * ECCFIC  #weighting correlation
    ECCFIC[!is.finite(ECCFIC)] <- 0
    v <- colnames(del0)
    final1 <- data.frame(v, ECCFIC)

    colnames(final1)[1:2] <- c("variable", "cor")
    # del1
    ECCFIC <- foreach(i = 1:ncol(del1), .combine = c) %do% {

        x <- del1[, c(i)]
        y <- time1
        rho_ech2(x, y, est = 1)
    }
    ECCFIC <- (nrow(del1)/n) * ECCFIC  #weighting correlation
    ECCFIC[!is.finite(ECCFIC)] <- 0
    v <- colnames(del1)

    final2 <- data.frame(v, ECCFIC)

    colnames(final2)[1:2] <- c("variable", "cor")

    # for whole
    ECCFIC <- foreach(i = 1:ncol(whole), .combine = c) %do% {

        x <- whole[, c(i)]
        y <- del_whole
        rho_ech2(x, y, est = 2)

    }
    ECCFIC <- (nrow(whole)/n) * ECCFIC  #weighting correlation

    ECCFIC[!is.finite(ECCFIC)] <- 0

    v <- colnames(whole)
    final3 <- data.frame(v, ECCFIC)
    colnames(final3)[1:2] <- c("variable", "cor")

    mar_step1 <- cbind(colnames(whole), final1[, c(2)] + final2[, c(2)], final3[, c(2)])
    mar_step1 <- data.frame(mar_step1)
    colnames(mar_step1)[1:3] <- c("X", "(time,X)|delta", "(delta,X)")
    # stopCluster(cl)


    mar_step1$`(time,X)|delta` <- as.numeric(mar_step1$`(time,X)|delta`)
    mar_step1$`(delta,X)` <- as.numeric(mar_step1$`(delta,X)`)
    mar_step1$rank_col2 <- rank(-mar_step1$`(time,X)|delta`)
    mar_step1$rank_col3 <- rank(-mar_step1$`(delta,X)`)
    mar_step1$min_rank <- apply(mar_step1[, c(4, 5)], 1, FUN = min)
    mar_step1$max_rank <- apply(mar_step1[, c(4, 5)], 1, FUN = max)
    mar_step1$`final rank` <- data.table::frank(list(mar_step1$min_rank, mar_step1$max_rank), ties.method = "dense")
    mar_step1 <- mar_step1[order(mar_step1$`final rank`, decreasing = FALSE), ]
    d_our <- mar_step1[, c(1, 2, 3, 8)]
    d_our
}
