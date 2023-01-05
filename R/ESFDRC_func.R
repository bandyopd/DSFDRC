#' ECCFIC-Screening with False Discovery Rate control
##' @title ESFDRC
##' @return The final set of important covariates with a specfic false discovery rate
##' @author Atika Farzana Urmi
##' @param data_n1 a data set with column 1= 'time',column2='delta' and rest are 'covariates' (can be obtained from the function split_data())
##' @param data_n2 a data set with column 1= 'time',column2='delta' and rest are 'covariates' (can be obtained from the function split_data())
##' Note: n1<n2
##' @param rand_num a random seed to reproduce the result
##' @param q a prespecified false discovery rate (usually .05 or .10)
##' @param s the number of covariates to be screened in the 1st step. Default is (n/log(n))
##' following Fan&Lv(2008) where 'n' i the number of rows in 'data_n1'
#' @export
ESFDRC_func <- function(data_n1, data_n2, rand_num, q, s = round(nrow(data_n1)/log(nrow(data_n1)), 0)) {
    p = rand_num
    # step1: screening for d covairates with n1 data
    n1_screen <- ECCFIC_screen(data_n1[, c(1)], data_n1[, c(2)], x_mat = data_n1[, c(3:ncol(data_n1))])
    # screen covariates
    d <- n1_screen[1:s, 1]

    ## intermediate step: create knockoff copies of 'd' covariates with n2 data

    copy_x <- copy_d_n2(data_n2[, which(colnames(data_n2) %in% d)], p)
    data_n2_copy <- cbind(data_n2[, c(1, 2)], copy_x)

    #### step3: Perform ECCFIC on original and knockoff copies of 'd' covarites
    data_original <- rbind(data_n1[, c(1, 2, which(colnames(data_n1) %in% d))], data_n2[, c(1, 2, which(colnames(data_n2) %in% 
        d))])
    data_knockoff <- rbind(data_n1[, c(1, 2, which(colnames(data_n1) %in% d))], data_n2_copy)

    new <- ECCFIC_screen(data_original[, c(1)], data_original[, c(2)], x_mat = data_original[, c(3:ncol(data_original))])
    new <- new[order(new$X), ]
    new2 <- ECCFIC_screen(data_knockoff[, c(1)], data_knockoff[, c(2)], x_mat = data_knockoff[, c(3:ncol(data_knockoff))])
    new2 <- new2[order(new2$X), ]
    #### final set of covariates based on important statistics

    compare <- cbind(new, new2)
    compare <- compare[, -c(4, 5, 8)]
    compare[, c(2, 3, 4, 5)] <- sapply(compare[, c(2, 3, 4, 5)], as.numeric)
    head(compare)  #check
    W1 <- compare[, c(2)] - compare[, c(4)]
    W2 <- compare[, c(3)] - compare[, c(5)]
    t1 <- knockoff.threshold(W1, fdr = q, offset = 0)  #10% FDR control, offset=0
    t2 <- knockoff.threshold(W2, fdr = q, offset = 0)


    d2 <- sort(union(which(W1 >= t1), which(W2 >= t2)))

    compare$X[c(d2)]  #final d2 covariates



    var <- as.vector(compare$X[c(d2)])
    var

}
