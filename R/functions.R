#' @title Select biomarker candidates
#'
#' @description A method that integrates differential expression (DE) analysis
#'   and differential network (DN) analysis to select biomarker candidates for
#'   survival time prediction.
#' @param x a data frame consists of data from group 1 and group 2.
#' @param class_label a binary array with 0: group 1; 1: group 2.
#' @param id an array of biomolecule ID.
#' @param partial logical. If TRUE the sparse differential network will be obtained
#'    by using partial correlation. If FALSE, correlation.
#' @param method a character string indicating which correlation coefficient is
#'    to be computed. One of "pearson" (default) or "spearman".
#' @param p_val a path to either a data frame or matrix that contains p-values
#'    or NULL.
#'
#' @inheritParams compute_cor
#' @inheritParams compute_par
#' @inheritParams compute_dns
#' @inheritParams permutation_cor
#' @inheritParams permutation_pc
#' @inheritParams permutation_thres
#' @inheritParams pvalue_logit
#' @inheritParams loglik_ave
#' @inheritParams choose_rho
#'
#'
#'
#' @return A .csv file containing the p-value, node degree, and activity score
#'    for each biomarker candidate
#'
#' @examples
#' select_sig(x = Met_GU, class_label = Met_Group_GU, id = Met_name_GU,
#'                                      partial = NULL, method = NULL, p_val = NULL)
#'
#' @importFrom glasso glasso
#' @importFrom utils write.csv read.table write.table txtProgressBar setTxtProgressBar
#' @importFrom stats qnorm cor quantile var sd glm
#' @importFrom graphics abline title plot lines
#'
#' @export


select_sig <- function(x = NULL, class_label = NULL, id = NULL,
    partial = NULL, method = NULL,  p_val = NULL) {


    data_bind <- rbind(x, class_label)
    raw_group_1 <- data_bind[,data_bind[nrow(data_bind),] == 0][1:(nrow(data_bind) - 1),]  # Group 1: p*n1
    raw_group_2 <- data_bind[,data_bind[nrow(data_bind),] == 1][1:(nrow(data_bind) - 1),]  # Group 2: p*n2

    p <- nrow(raw_group_2)
    n_group_2 <- ncol(raw_group_2)
    n_group_1 <- ncol(raw_group_1)

    # Z-transform the data for group-specific normalization
    data_group_1 <- scale(t(raw_group_1)) # Group 1: n1*p
    data_group_2 <- scale(t(raw_group_2)) # Group 2: n2*p
    cov_group_1 <- var(data_group_1)
    cov_group_2 <- var(data_group_2)


    # If the users want the sparse network using partial correlation
    if (partial == TRUE) {
        set.seed(100)
        ## Apply grapihcal LASSO given data_group_1 and data_group_2
        #  Group 1 first
        n_fold <- 5 # number of folds
        rho <- exp(seq(log(0.6), log(0.01), length.out = 20))
        # draw error curve
        error <- choose_rho(data_group_1, n_fold, rho)
        title(main = "Group 1: Error curve using corss validation")
        # chosse optimal rho
        rho[error$log.cv == min(error$log.cv)] # rho based on minimum rule
        abline(v = rho[error$log.cv == min(error$log.cv)], col = "red", lty = 3)
        # one standard error rule
        abline(h = min(error$log.cv) + error$log.rho[error$log.cv == min(error$log.cv)], col = "blue")
        # rhos that are under blue line
        rho_under_blue <- rho[error$log.cv < min(error$log.cv) + error$log.rho[error$log.cv == min(error$log.cv)]]
        # rhos that are on the right of the red line
        rho_right_red <- rho[rho > rho[error$log.cv == min(error$log.cv)]]
        #rhos that are under blue line and to the right of the red line
        rho_one_standard <- intersect(rho_under_blue, rho_right_red )
        ##############################################################
        # User interaction in console for selecting rho for group 1
        ##############################################################
        min_rule <- rho[error$log.cv == min(error$log.cv)]   # rho based on minumum rule

        # Provide users with an opportunity to interact within the function
        print("The list of rhos for group 1:")
        print(rev(rho))

        rho <- readline("Choose your own regularization parameter rho for group 1? [y/n]: ")
        if (rho == "y") {
            your_rho <- readline(prompt = "Enter your own choice of rho: ")
            rho_group_1_opt <- as.numeric(your_rho)
            print(rho_group_1_opt)
        } else if (rho == "n") {
        rho_based_on_rule <- readline(prompt = "rho based on minimum rule? [y/n]: ")
            if (rho_based_on_rule == "y") {
            rho_group_1_opt <- min_rule
            print(rho_group_1_opt)
            } else {                # select rho based on one standard rule
            rho_group_1_opt <- max(rho_one_standard)
            print(rho_group_1_opt)
            }
        }

        # perform gLASSO
        pre_group_1 <- glasso(cov_group_1, rho = rho_group_1_opt)
        rm(n_fold, rho_under_blue, rho_right_red, rho, error, rho_one_standard, min_rule)
        # examine the precision matrix
        thres <- 1e-3
        sum(abs(pre_group_1$wi) > thres)
        pre_group_1$wi[1:10, 1:10]
        # compute partial correlation
        pc_group_1 <- compute_par(pre_group_1$wi)
        # examine the partial correlation matrix
        sum(abs(pc_group_1) > thres)
        pc_group_1[1:10, 1:10]
        rm(pre_group_1, thres)

        ## Group 2 second
        n_fold <- 5 # number of folds
        rho <- exp(seq(log(0.6), log(0.01), length.out = 20))
        # draw error curve
        error <- choose_rho(data_group_2, n_fold, rho)
        title(main = "Group 2: Error curve using corss validation")
        # chosse optimal rho
        rho[error$log.cv == min(error$log.cv)] # rho based on minimum rule
        abline(v = rho[error$log.cv == min(error$log.cv)], col = "red", lty = 3)
        # one standard error rule
        abline(h = min(error$log.cv) + error$log.rho[error$log.cv == min(error$log.cv)], col = "blue")
        # rhos that are under blue line
        rho_under_blue <- rho[error$log.cv < min(error$log.cv) + error$log.rho[error$log.cv == min(error$log.cv)]]
        # rhos that are on the right of the red line
        rho_right_red <- rho[rho > rho[error$log.cv == min(error$log.cv)]]
        #rhos that are under blue line and to the right of the red line
        rho_one_standard <- intersect(rho_under_blue, rho_right_red )
        #############################################################
        # User interaction in console for selecting rho for group 2
        #############################################################
        min_rule <- rho[error$log.cv == min(error$log.cv)]   # rho based on minumum rule

        # Provide users with an opportunity to interact within the function
        print("The list of rhos for group 2:")
        print(rev(rho))

        rho <- readline("Choose your own regularization parameter rho for group 2? [y/n]: ")
        if (rho == "y") {
            your_rho <- readline(prompt = "Enter your own choice of rho: ")
            rho_group_2_opt <- as.numeric(your_rho)
            print(rho_group_2_opt)
        } else if (rho == "n") {
            rho_based_on_rule <- readline(prompt = "rho based on minimum rule? [y/n]: ")
            if (rho_based_on_rule == "y") {
                rho_group_2_opt <- min_rule
                print(rho_group_2_opt)
            } else {
                rho_group_2_opt <- max(rho_one_standard)
                print(rho_group_2_opt)
            }
        }

        # perform gLASSO
        pre_group_2 <- glasso(cov_group_2, rho = rho_group_2_opt)
        rm(n_fold, rho_under_blue, rho_right_red, rho, error, rho_one_standard, min_rule)
        # examine the precision matrix
        thres <- 1e-3
        sum(abs(pre_group_2$wi) > thres)
        pre_group_2$wi[1:10,1:10]
        # compute partial correlation
        pc_group_2 <- compute_par(pre_group_2$wi)
        # examine the partial correlation matrix
        sum(abs(pc_group_2) > thres)
        pc_group_2[1:10,1:10]
        rm(pre_group_2, thres)

        ## Build differential partial correlation networks
        diff <- pc_group_2 - pc_group_1  # from group 1 to group 2
        thres = 1e-3
        sum(abs(diff) > thres)
        diff[1:10, 1:10]

        ## Permutation test using partial correlation
        num_of_permutations_pc <- readline(prompt = "Enter your desired number of permutations to build differential network using partial correlation [Default: 1000]: ")
        m <- as.numeric(num_of_permutations_pc)
        diff_p <- permutation_pc(m, p, n_group_1, n_group_2, data_group_1, data_group_2, rho_group_1_opt, rho_group_2_opt)
        rm(m, thres, rho_group_1_opt, rho_group_2_opt, num_of_permutations_pc)

    } else {       # Obtain the sparse network using correlation


        # Compute correlaton matrix for each group
        cor <- compute_cor(data_group_2, data_group_1, type_of_cor = method)    # default is pearson correlation

        # Get the correlation matrix
        cor_group_2 <- cor$Group2
        cor_group_1 <- cor$Group1

        # examine the correlation matrix
        thres <- 1e-3
        sum(abs(cor_group_2) > thres)
        cor_group_2[1:10, 1:10]
        sum(abs(cor_group_1) > thres)
        cor_group_1[1:10, 1:10]
        rm(thres)

        # Build differential correlation networks
        diff <- cor_group_2 - cor_group_1 # from group 1 to group 2
        thres = 1e-3
        sum(abs(diff) > thres)
        diff[1:10, 1:10]

        # Permutation test
        num_of_permutations_c <- readline(prompt = "Enter your desired number of permutations to build differential network using correlation [Default: 1000]: ")
        m <- as.numeric(num_of_permutations_c)
        diff_p <- permutation_cor(m, p, n_group_1, n_group_2, data_group_1, data_group_2, type_of_cor = method)
        rm(m, thres, num_of_permutations_c)
    }




    # calculate differential network connections
    thres_left <- 0.025
    thres_right <- 0.975
    significant_thres <- permutation_thres(thres_left, thres_right, p, diff_p)
    rm(thres_left, thres_right)

    # get binary matrix
    significant_thres_p <- significant_thres$positive
    significant_thres_n <- significant_thres$negative
    binary_link <- matrix(0, p, p) # binary connection
    binary_link[diff < significant_thres_n] <- -1
    binary_link[diff > significant_thres_p] <- 1
    weight_link <- matrix(0, p, p) # weight connection
    weight_link[diff < significant_thres_n] <- diff[diff < significant_thres_n]
    weight_link[diff > significant_thres_p] <- diff[diff > significant_thres_p]
    sum(diff < significant_thres_n)
    sum(diff > significant_thres_p)
    binary_link[1:10, 1:10]
    weight_link[1:10, 1:10]
    rowSums(abs(binary_link)) # node degree for differential networks
    rm(diff_p)

    # Convert adjacent matrix into edge list
    edge <- matrix(0, (sum(diff < significant_thres_n) + sum(diff > significant_thres_p)) / 2, 4)
    k <- 1
    for (i in 1:(nrow(binary_link) - 1)) {
        for (j in (i + 1) : nrow(binary_link)) {
            if(binary_link[i, j] != 0) {
                edge[k, 1] <- i
                edge[k, 2] <- j
                edge[k, 3] <- binary_link[i, j]
                edge[k, 4] <- binary_link[i, j]
                k <- k + 1
            }
        }
    }
    edge_dn <- data.frame("Met1" = edge[, 1], "Met2" = edge[, 2], "Binary" = edge[, 3], "Weight" = edge[, 4])
    write.csv(edge_dn, file = "Met_dn.csv", quote = FALSE, row.names = FALSE)
    rm(k, edge, edge_dn , significant_thres)

    # If the p-value table is not provided by users
    if (is.null(p_val) == TRUE) {
        # Calculate p-values using logistic regression if p-values are not provided by users
        pvalue <- pvalue_logit(x, class_label, id)
        p.value <- pvalue$p.value
    } else {     # If the p-value matrix is provided
        pvalue <- read.table(p_val, sep = ",", header = TRUE)
        p.value <- pvalue$p.value           # Extract p-values from the table provided
    }


    # trasfer p-value to z-score
    z_score <- abs(qnorm(1 - p.value/2))
    # calculate differntial network score
    dn_score <- compute_dns(binary_link, z_score)


    indeed_df <- cbind(pvalue, rowSums(abs(binary_link)), dn_score )


    # save the p-value, node degree and activity score of each biomarker candidate
    write.table(indeed_df, file = "INDEED_result.csv", sep=",", quote = FALSE,
        row.names = FALSE, col.names = c("ID", "P-value", "Node Degree", "Activity Score"))
}




