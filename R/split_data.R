#' ECCFIC_screen ranks the predictors based on their ECCFIC-correlation with survival time and censoring indicator
##' @title Split_data
##' @return Splits the data into two subsets (n1 & n2 with n1<n2) to perform screening with FDR control
##' @author Atika Farzana Urmi
##' @param data the whole dataset to be splitted into two subsets. The 1st column has to have the name 'time'
##' and the 2nd column has to have the name 'delta'. Rest will be the covariates.
##' @param n1  desired number of rows in 1st subset.Final n1 might differ slightly after adjusting for censoring rate in both data.
##' @param n2  desired number of rows in 2nd subset.Final n2 might differ slightly after adjusting for censoring rate in both data.
##' @param rand_num setting a random number to reproduce the result
##' @return A list od two data sets
#' @export
split_data <- function(data, n1, n2, rand_num) {
    # need to keep censoring rate similar in both data
    data$id <- seq(1:nrow(data))
    data <- data[, c(ncol(data), 1:(ncol(data) - 1))]
    event <- data[data$delta == 1, ]
    censored <- data[data$delta == 0, ]
    ##############################################################################
    per1 <- round(n1/(nrow(data)), 2)
    set.seed(rand_num)
    ################################ keeping censoring rate similar in both data#########
    p1 <- round(nrow(event) * per1, 0)
    p2 <- round(nrow(censored) * per1, 0)
    z <- sample(event$id, p1)
    z2 <- sample(censored$id, p2)
    z1 <- sort(c(z, z2))


    data_n1 <- data[data$id %in% z1, -c(1)]
    data_n2 <- data[!data$id %in% z1, -c(1)]
    return(list(data_n1, data_n2))
}

