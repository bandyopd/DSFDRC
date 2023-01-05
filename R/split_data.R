#' ECCFIC_screen ranks the predictors based on their ECCFIC-correlation with survival time and censoring indicator
##' @title Split_data
##' @return Splits the data into two subsets (n1 & n2 with n1<n2) to perform screening with FDR control
##' @author Atika Farzana Urmi
##' @param data the whole dataset to be splitted into two subsets
##' @param n1 number of rows in 1st subset
##' @param n2 number of rows in 2nd subset
##' @param rand_num setting a random number to reproduce the result
##' @return A list od two data sets
#' @export
split_data <- function(data, n1, n2, rand_num) {
    # need to keep censoring rate similar in both data
    data$id <- seq(1:nrow(data))
    data <- data[, c(ncol(data), 1:(ncol(data) - 1))]
    event <- data[data$delta == 1, ]  #220 row
    censored <- data[data$delta == 0, ]  #296 row
    ############################################################################## 


    # rand_num<- floor(runif(1, min=100, max=9999))#4452
    rand_num <- 7519  #8409#6882#2660#1548#2208 #9568#2925 #8409
    ## split data into n1 and n2.n1<n2. (48% and 52%, larger n1,larger d1 as d1=n1/log(n1))
    per1 <- round(n1/(nrow(data)), 2)
    set.seed(rand_num)
    #################################################################################### 33 sample observations for n1 data z1<-sample(data$id,248) #516*.48=248
    p1 <- round(nrow(event) * per1, 0)
    p2 <- round(nrow(censored) * per1, 0)
    z <- sample(event$id, p1)  #220*.48=106
    z2 <- sample(censored$id, p2)  #296*.48=142
    z1 <- sort(c(z, z2))


    data_n1 <- data[data$id %in% z1, -c(1)]
    data_n2 <- data[!data$id %in% z1, -c(1)]
    return(list(data_n1, data_n2))
}

