## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=FALSE, results='hide', message=FALSE, fig.height=5, fig.width=6.5----
## Draw error curve
devtools::load_all()
set.seed(100)
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
  for(i in 1 : length(rho)){
    Sys.sleep(0.1)
    # perform n_fold cross validation
    loglik <- c()
    for(j in 1 : n_fold){
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
data_bind <- rbind(Met_GU, Met_Group_GU)
raw_group_1 <- data_bind[,data_bind[nrow(data_bind),] == 0][1:(nrow(data_bind) - 1),]  # Group 1: p*n1
raw_group_2 <- data_bind[,data_bind[nrow(data_bind),] == 1][1:(nrow(data_bind) - 1),]  # Group 2: p*n2
# Z-transform the data for group-specific normalization
data_group_1 <- scale(t(raw_group_1)) # Group 1: n1*p
data_group_2 <- scale(t(raw_group_2)) # Group 2: n2*p
n_fold <- 5 # number of folds
rho = exp(seq(log(0.6), log(0.01), length.out = 20))
# draw error curve
error <- choose_rho(data_group_1, n_fold, rho)
title(main = "Group 1: Error curve using corss validation")
# chosse optimal rho
rho[error$log.cv == min(error$log.cv)] # rho based on minimum rule
abline(v = rho[error$log.cv == min(error$log.cv)], col = "red", lty = 3)
# one standard error rule
abline(h = min(error$log.cv) + error$log.rho[error$log.cv == min(error$log.cv)], col = "blue")

## ---- echo=FALSE, results='hide', message=FALSE, fig.height=5, fig.width=6.5----
## Draw error curve
set.seed(100)
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
  for(i in 1 : length(rho)){
    Sys.sleep(0.1)
    # perform n_fold cross validation
    loglik <- c()
    for(j in 1 : n_fold){
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

# draw error curve
error <- choose_rho(data_group_2, n_fold, rho)
title(main = "Group 2: Error curve using corss validation")
# chosse optimal rho
rho[error$log.cv == min(error$log.cv)] # rho based on minimum rule
abline(v = rho[error$log.cv == min(error$log.cv)], col = "red", lty = 3)
# one standard error rule
abline(h = min(error$log.cv) + error$log.rho[error$log.cv == min(error$log.cv)], col = "blue")

## ---- echo=FALSE, results='asis'-----------------------------------------
x <- read.csv("../inst/INDEED_result.csv", header = FALSE)
names(x) <- NULL
knitr::kable(head(x, 11))

## ---- echo=FALSE, results='asis'-----------------------------------------

Met_GU <- get("Met_GU")[1:10, 1:10]   # Displying the fisrt 10*10 data in Met_GU
knitr::kable(Met_GU, caption = "**Met_GU**")

## ---- echo=FALSE, results='asis'-----------------------------------------
p_val <- read.table("../inst/pvalue_M_GU.csv", sep = ",", header = TRUE)
knitr::kable(head(p_val, 10), caption = "**P-values**")     # display only the first 10 metabolites associated with their p-values

## ---- echo=FALSE, results='asis'-----------------------------------------


Met_Group_GU <- get("Met_Group_GU")[1:10]
knitr::kable(Met_Group_GU)   # Displying the fisrt 10 data in Met_Group_GU

## ---- echo=FALSE, results='asis'-----------------------------------------

Met_name_GU <- get("Met_name_GU")[1:10]
knitr::kable(Met_name_GU)     # Displying the fisrt 10 data in Met_name_GU

## ---- echo=FALSE, results='asis'-----------------------------------------
x <- read.csv("../inst/INDEED_result_pearson.csv", header = FALSE)
names(x) <- NULL
knitr::kable(head(x, 11))

## ---- echo=FALSE, results='asis'-----------------------------------------
x <- read.csv("../inst/INDEED_result_spearman.csv", header = FALSE)
names(x) <- NULL
knitr::kable(head(x, 11))

