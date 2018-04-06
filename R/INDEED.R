#' INDEED: A package for biomarker candidate prioritization.
#'
#' The INDEED package provides two categories of important functions:
#' select_sig and load_sample_data.
#'
#' @section select_sig function:
#' The select_sig function provides the users with few options. By setting \bold{partial} = TRUE,
#' the differential network is obtained using partial correlation. This function is also interactive
#' in the sense that it provides users with few options in the console for the selection of regularization
#' parameter rho as well number of permutations. With \bold{partial} = FALSE, correlation will be used to get the differential network.
#' In addtion, the default method for the correlation part is "pearson". Otherwise, use \bold{method} = "spearman"
#' to obtain Spearman correlation coefficient. The \bold{p_val}: gives users more flexibility by allowing
#' the users to provide their own p-value tables. If NULL, p-values will be calculated for the users
#' using logistic regression.
#'
#' @section load_sample_data function:
#' The load_sample_data function will load the sample data sets to run the select_sig function.
#'
#' @docType package
#' @name INDEED
NULL
