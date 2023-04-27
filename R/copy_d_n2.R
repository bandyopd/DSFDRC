# creating knockoff copies of covairates screened in 1st step.Here we put the bigger data (n2) with only the covariates
# selected in the screening step
##' @title copy_d_n2
##' @return a data frame of time,delta and knockoff copies of original covariates
##' @author Atika Farzana Urmi
##' @param data_d_n2 A data frame of time,delta and original 'd' covariates
##' @param rand_num a random number which was used in the data splitting step
##' @noRd
copy_d_n2 <- function(data_d_n2, rand_num) {
    xo <- as.matrix(data_d_n2[, -c(1, 2)])  #original x
    # normalizing constant
    norm_x <- apply(xo, 2, function(x) sqrt(sum(x^2)))

    # normalizing original X
    xo <- sweep(xo, MARGIN = 2, STATS = norm_x, FUN = "/")  #now diag(t(xo)%*%(xo))=1

    # creating knockoff

    set.seed(rand_num)
    xk <- create.fixed(xo)$Xk  ##diag(t(xk)%*%(xk))= 1
    # getting back to original scale
    xk <- sweep(xk, MARGIN = 2, STATS = norm_x, FUN = "*")  #daigonals are not one any more
    xk <- data.frame(xk)
    copy_data <- data.frame(cbind(data_d_n2[, c(1, 2)], xk))
    copy_data
}
