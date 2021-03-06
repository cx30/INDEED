## Functions defined to help the functions.R
# This file contains several very important functions to run select_sig() in functions.R


#' @title Compute the correlation
#'
#' @description Compute either pearson or spearman correlation coefficient.
#' @param data_group_2 a n*p matrix.
#' @param data_group_1 a n*p matrix
#' @param type_of_cor if NULL, pearson correlation coefficient will be calculated.
#'   Otherwise, a character string "spearman" to calculate spearman correlation
#'   coefficient.
#'
#' @return A list of correlation matrix for both group 1 and group 2
#'


# Compute Pearson correlation or Spearman correlation
compute_cor <- function(data_group_2, data_group_1, type_of_cor) {
    if (is.null(type_of_cor) || type_of_cor == "pearson") {
        cor_group_2 <- cor(data_group_2, method = "pearson")
        cor_group_1 <- cor(data_group_1, method = "pearson")

    } else if (type_of_cor == "spearman") {
        cor_group_2 <- cor(data_group_2, method = "spearman")
        cor_group_1 <- cor(data_group_1, method = "spearman")
    }
    cor <- list("Group2" = cor_group_2, "Group1" = cor_group_1)
}



#' @title Compute the partial correlation
#'
#' @description Compute the partial correlation coefficient.
#' @param pre_inv an inverse covariance matrix.
#'
#' @return An \eqn{n by n} partial correlation matrix
#'

## Compute partial correlation
compute_par <- function(pre_inv) {
    p <- nrow(pre_inv)
    pc <- matrix(0, p, p) # partial correlation
    for (i in 1:(p - 1)) {
        for (j in (i + 1) : p) {
            pc[i, j] = -pre_inv[i, j] / sqrt(pre_inv[i, i] * pre_inv[j, j])
            pc[j, i] = pc[i, j]
        }
    }
    return(pc)
}





#' @title Permutations to build differential network using correlation
#'
#' @description A permutation test that randomly permutes the sample labels in distinct
#'     biological groups for each biomolecule. The difference in each paired biomolecule
#'     is considered significant if it falls into the 2.5% tails on either end of the empirical
#'     distribution curve.
#' @param m number of permutations.
#' @param p number of biomarker candidates.
#' @param n_group_1 number of subjects in group 1.
#' @param n_group_2 number of subjects in group 2.
#' @param data_group_1 a \eqn{n*p} matrix or data.frame containing group 1 data.
#' @param data_group_2 a \eqn{n*p} matrix of data.frame containing group 2 data.
#' @param type_of_cor if NULL, pearson correlation coefficient will be calculated.
#'     Otherwise, a character string "spearman" to calculate spearman correlation
#'     coefficient.
#'
#' @return A multi-dimensional matrix that contains the permutation results

# Permutation to build differential network using correlation
permutation_cor <- function(m, p, n_group_1, n_group_2, data_group_1, data_group_2, type_of_cor) {
    diff_p <- array(0, dim = c(m, p, p))
    pb <- txtProgressBar(min = 0, max = m, style = 3)
    for (t in 1 : m) {
        Sys.sleep(0.1)
        data_group_1_p <- matrix(0, n_group_1, p)
        for (i in 1 : p) {
            data_group_1_p[, i] <- data_group_1[sample(n_group_1), i]
        }
        data_group_2_p <- matrix(0, n_group_2, p)
        for (i in 1 : p) {
            data_group_2_p[, i] <- data_group_1[sample(n_group_2), i]
        }

    if (is.null(type_of_cor)) {
        cor_group_2_p <- cor(data_group_2_p, method = "pearson")
        cor_group_1_p <- cor(data_group_1_p, method = "pearson")
    } else {
        cor_group_2_p <- cor(data_group_2_p, method = "spearman")
        cor_group_1_p <- cor(data_group_1_p, method = "spearman")
    }
        diff_p[t, , ] <- cor_group_2_p - cor_group_1_p

        # update progress bar
        setTxtProgressBar(pb, t)
    }
    close(pb)
    return(diff_p)
}


#' @title Permutations to build differential network using partial correlation
#'
#' @description A permutation test that randomly permutes the sample labels in distinct
#'     biological groups for each biomolecule. The difference in paired partial correlation
#'     is considered significant if it falls into the 2.5% tails on either end of the empirical
#'     distribution curve.
#' @param m number of permutations.
#' @param p number of biomarker candidates.
#' @param n_group_1 number of subjects in group 1.
#' @param n_group_2 number of subjects in group 2.
#' @param data_group_1 a \eqn{n*p} matrix or data.frame containing group 1 data.
#' @param data_group_2 a \eqn{n*p} matrix of data.frame containing group 2 data.
#' @param rho_group_1_opt optimal tuning parameter to sparse the differential network for group 1
#' @param rho_group_2_opt optimal tuning parameter to sparse the differential network for group 2
#'
#' @return A multi-dimensional matrix that contains the permutation results

## Permutation to build differential network using partial correlation
permutation_pc <- function(m, p, n_group_1, n_group_2, data_group_1, data_group_2, rho_group_1_opt, rho_group_2_opt) {
    diff_p <- array(0, dim = c(m, p, p))
    pb <- txtProgressBar(min = 0, max = m, style = 3)
    for(t in 1 : m) {
        Sys.sleep(0.1)
        data_group_1_p <- matrix(0, n_group_1, p)
        for(i in 1 : p) {
            data_group_1_p[, i] <- data_group_1[sample(n_group_1), i]
        }
        data_group_2_p <- matrix(0, n_group_2, p)
        for(i in 1 : p) {
            data_group_2_p[, i] <- data_group_2[sample(n_group_2), i]
        }
        per_group_1 <- glasso(var(data_group_1_p), rho = rho_group_1_opt)
        per_group_2 <- glasso(var(data_group_2_p), rho = rho_group_2_opt)
        pc_group_1_p <- compute_par(per_group_1$wi)
        pc_group_2_p <- compute_par(per_group_2$wi)
        diff_p[t, , ] <- pc_group_2_p - pc_group_1_p
        # update progress bar
        setTxtProgressBar(pb, t)
    }
    close(pb)
    return(diff_p)
}




#' @title Calculate the positive and negative threshold based on the permutation result
#'
#' @description Calculate the positive and negative threshold based on the permutation result.
#'
#' @param thres_left 2.5 percent left tails.
#' @param thres_right 2.5 percent right tails.
#' @param p number of biomarker candidates.
#' @param diff_p permutation results.
#'
#' @return A list of positive and negative threshold

# Calculate the positive and negative threshold based on the permutation result
permutation_thres <- function(thres_left, thres_right, p, diff_p) {
    significant_thres_p <- matrix(0, p, p)
    significant_thres_n <- matrix(0, p, p)
    for (i in 1 : (p-1)) {
        for (j in (i + 1) : p) {
            significant_thres_n[i, j] <- quantile(diff_p[, i, j], probs = thres_left)
            significant_thres_n[j, i] <- significant_thres_n[i, j]
            significant_thres_p[i, j] <- quantile(diff_p[, i, j], probs = thres_right)
            significant_thres_p[j, i] <- significant_thres_p[i, j]
        }
    }
    significant_thres <- list("positive" = significant_thres_p, "negative" = significant_thres_n)
    return(significant_thres)
}




#' @title Calculate differential network score
#'
#' @description Calculate differential network score.
#'
#' @param binary_link binary correlation matrix with 1 indicating positive correlation and -1
#'     indicating negative correlation for each biomolecular pair.
#' @param z_score converted from p-value.
#'
#' @return An activity score associated with each biomarker candidate

# Calculate differential network score
compute_dns <- function(binary_link, z_score) {
    # get adjacent matrix
    diff_d <- abs(binary_link)
    # set diagonal elements to 1
    diag(diff_d) <- 1
    # compute differential network score for each row
    dns <- apply(diff_d, 1, function(x, y = z_score) sum(y[which(x == 1)]))
    return(dns)
}




#' @title Obtain p-values using logistic regression
#'
#' @description Calculate p-values using logistic regression.
#'
#' @param x a data frame consists of data from group 1 and group 2.
#' @param class_label a binary array indicating 0: group 1; 1: group 2.
#' @param Met_name an array of ID.
#'
#' @return p-values

pvalue_logit <- function(x, class_label, Met_name) {
    data_tp <- as.data.frame(t(x))    # n*p
    class_label_tp <- as.data.frame(t(class_label))
    # attach metabolites ID and class label in the data set
    X_df <- cbind(data_tp, class_label_tp)
    colnames(X_df)[1:(ncol(X_df)-1)] <- Met_name
    colnames(X_df)[ncol(X_df)] <- "Class"
    glm.fit <- glm(Class ~. , family = "binomial", data = X_df)
    ## Sort metabolites based on their p-values
    pvalue <- summary(glm.fit)$coefficients[,4][2:ncol(X_df)]
    pvalue_df <- data.frame("ID" = Met_name, "p.value" = pvalue)
    return(pvalue_df)
}



#' @title Create log likelihood error function
#'
#' @description Calculate log likelihood error function.
#'
#' @param data a matrix or data.frame.
#' @param theta a precision matrix.
#'
#' @return log likelihood error function

## Create log likelihood error function
loglik_ave <- function(data, theta){
    loglik <- c()
    loglik <- log(det(theta)) - sum(diag(var(data) %*% theta))
    return(-loglik)
}





#' @title Draw error curve
#'
#' @description Draw error curve using cross-validation.
#'
#' @param data a matrix.
#' @param n_fold specify n to n-fold cross_validation.
#' @param rho multiple regularization parameter values to be evalueated in terms of errors.
#'
#' @return a list of errors and their corresponding \eqn{log(rho)}

## Draw error curve
choose_rho <- function(data, n_fold, rho) {
    # randomly shuffle the data
    Data <- data[sample(nrow(data)), ]
    # create n_fold equally size folds
    folds <- cut(seq(1, nrow(Data)), breaks = n_fold, labels = FALSE)
    # tune parameters
    d <- ncol(Data)
    loglik_cv <- c()
    loglik_rho <- c()
    pb <- txtProgressBar(min = 0, max = length(rho), style = 3) # create progress bar
    for (i in 1 : length(rho)) {
        Sys.sleep(0.1)
        # perform n_fold cross validation
        loglik <- c()
        for (j in 1 : n_fold) {
            # segement your data by fold using the which() function
            testIndexes <- which(folds == j, arr.ind = TRUE)
            testData <- Data[testIndexes, ]
            trainData <- Data[-testIndexes, ]
            # use test and train data partitions however you desire...
            cov <- var(trainData) # compute the covariance matrix
            pre<- glasso(cov, rho = rho[i])
            loglik <- c(loglik, loglik_ave(testData, pre$wi))
        }
    loglik_cv <- c(loglik_cv, sum(loglik) / n_fold)
    loglik_rho <- c(loglik_rho, sd(loglik) / sqrt(n_fold))
    setTxtProgressBar(pb, i) # update progress bar
    }
    close(pb)
    plot(rho, loglik_cv, xlab = expression(lambda), ylab = "Error")
    lines(rho, loglik_cv)
    error <- list("log.cv" = loglik_cv, "log.rho" = loglik_rho)
    return(error)
}






