#' This function returns the threshold value for the FDR control part based on the given important statistics (W1,W2)
##' @title threshold for aKIDS
##' @return a value for a prespecified FDR to seperate active and inactive set of predictors in aKIDS procedure
##' @author Chenlu Ke, Atika Urmi
##' @param W1 difference in omega1 between x (original) and xk (kncokoff) for the d covariates selected in screening step
##' @param W2 difference in omega2 between x (original) and xk (kncokoff) for the d covariates selected in screening step
##' Note: n1<n2
##' @param fdr a single or vector of prespecified false discovery rate for which the threshold to be calculated
#' @export
threshold <- function(W1, W2, fdr) {
    df = expand.grid(t1 = c(abs(W1), Inf), t2 = c(abs(W2), Inf))  #all possible combinations of two vectors
    df <- df[order(df$t1, df$t2), ]
    rW <- (rank(-W1) > rank(-W2)) + 1  #lower value, higher rank, coding values to 1/2 instead of 0/1 using +1
    th <- matrix(0, length(fdr), 2)
    offset = 1
    ratio <- sapply(1:nrow(df), function(i) (offset + sum(W1 <= -df[i, 1] | W2 <= -df[i, 2]))/sum(W1 >= df[i, 1] | W2 >= df[i, 2]))
    for (i in 1:length(fdr)) {
        ok <- which(ratio <= fdr[i])
        if (length(ok) == 0) {
            th[i, ] <- NA
        } else {
            # df.ok <- df[ok,,drop=F]
            df.ok <- df[ok, ]
            df.ok1 <- do.call(rbind, by(df.ok, list(df.ok[, 1]), head, n = 1))  #df is sorted from low to high,for each set, #n=1 gives one value, head=top
            df.ok1s <- df.ok1[order(df.ok1[, 2], df.ok1[, 1]), ]  #for each t2, we get t2, we get lowest to highest t1
            df.ok2 <- do.call(rbind, by(df.ok1s, list(df.ok1s[, 2]), head, n = 1))  #similarly for each set of sorted t2, we get 1(lowest) value for t1
            mes <- df.ok2  # minimal elements
            if (nrow(df.ok2) > 1) {
                rm <- NULL
                for (j in 2:nrow(df.ok2)) {
                  if (df.ok2[j, 1] > min(df.ok2[1:(j - 1), 1])) {
                    rm = c(rm, j)
                  }
                }
                if (!is.null(rm)) {
                  mes = mes[-rm, ]
                }
            }
            row.names(mes) <- NULL
            dt.tmp <- apply(mes, 1, function(ti) {
                d2.tmp <- which(W1 >= ti[1] | W2 >= ti[2])
                W <- cbind(W1, W2)
                mean(W[cbind(d2.tmp, rW[d2.tmp])])  #d2.tmp is the row number for which W1 or W2 > mes.
                # rw[d2.tmp] gives 1/2. if 1, W1<W2 and W1 column is selected. otherwise W2<W1 and W2 is selected
            })
            th[i, ] <- c(mes[which.max(dt.tmp), ][, 1], mes[which.max(dt.tmp), ][, 2])  #this is if mes have more than one row
            # d2[[i]] <- unique(as.vector(d1[which(W1>=th$t1|W2>=th$t2)]))
        }
        # d2[[i]] <- unique(as.vector(unlist(sapply(ok, function(i) d1[which(W1>=df[i,1]|W2>=df[i,2])]))))
    }

    return(th)

}

