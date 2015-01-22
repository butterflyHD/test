  # ST758 Homework 9
  # Meng Yang

  # clear up
  rm(list = ls())
  
  library(reshape)
  library(MASS)
  library(mc2d)
  library(compiler)
  library(knitr)
  library(ggplot2)
  # Question 1:
  ddirmult <- function(x, alpha, log = FALSE) {
    
    # check x is count data
    if (any(x < 0))
      stop("Error: x should be nonnegative count data")
    # check positivity of alpha
    if (any(alpha < 0))
      stop("Error: Dirichlet-multinomial parameter alpha should be nonnegative")
    
    # check dimensions of input arguments
    if (is.vector(alpha) && length(alpha) > 1) {
      if (is.vector(x) && length(x) > 1) {
        if (length(x) != length(alpha)) {
          stop("Error: sizes of x and alpha do not match.")
        } else {
          # expand x to a matrix of matching size with alpha
          x <- matrix(x, 1, length(x)) 
        }
      } else if (is.vector(x) && length(x) <= 1) {
        stop("Error: x can not be a scalar")
      }
      # expand alpha to a matrix of matching size with x
      alpha <- matrix(alpha, nrow = nrow(x), ncol = length(alpha), byrow = TRUE)
    }
    if (any(dim(alpha) != dim(x)))
      stop("Error: dimensions of alpha and x do not match")
    
    # compute log-densities
    alphaRowSums <- rowSums(alpha)
    xRowSums <- rowSums(x)
    # lgamma(0 + 0) - lgamma(0) will produce NaNs.
    # This fixes x_{ij} = alpha_{ij} = 0 cases.
    # Is there better way to deal with this?
    alpha[(x == 0) & (alpha == 0)] <- 1
    # assemble log-likelihood
    # lgamma(0) throws a lot of warnings but they are valid
    logl <- suppressWarnings(
      lfactorial(xRowSums) + rowSums(lgamma(x + alpha)) + 
        lgamma(alphaRowSums) - (rowSums(lfactorial(x)) + 
                                  rowSums(lgamma(alpha)) + 
                                  lgamma(alphaRowSums + xRowSums))
    )
    # Deal with alphaRowSums == 0 cases
    # fix the lgamma(0 + 0) - lgamma(0) produces NaN here
    logl[(xRowSums == 0) & (alphaRowSums == 0)] <- 0
    logl[(xRowSums > 0) & (alphaRowSums == 0)] <- - Inf
    
    # output
    if (log)
      return(logl)
    else
      return(exp(logl))
  }
  
  dirmultfit <- function(X, weights = NULL, alpha0 = NULL, 
                         tolfun = 1e-6, maxiters = 100, display = FALSE) {
    
    # remove data points with batch size 0
    rsum <- rowSums(X)
    if (any(rsum == 0)) { 
      rmv <- sum(rsum == 0)
      message(paste( "Warning: ", rmv,
                     " rows are removed because the row sums are 0"))
    }
    
    # remove the bins with no observations
    csum <- colSums(X)
    if (any(csum == 0)) { 
      rmv <- sum(csum == 0)
      message(paste("Warning: ", rmv,
                    " columns are removed because the column sums are 0"))
    }
    
    # cleaned up data
    data <- X[rsum != 0, csum != 0]
    N <- nrow(data)       ## Sample size
    d <- ncol(data)           ## Number of parameters
    m <- rowSums(data)  ## batch sizes
    
    # set starting points
    if (is.null(alpha0)) {
      # method of moment estimate
      rho <- sum(colSums((data / m)^2) / (colSums(data / m)))
      alpha0 <- as.vector(colSums(data / m) * (d - rho) / (rho - 1) / N)
      alpha0[alpha0 <= 0] = 1e-6
    } else {
      alpha0 <- alpha0[csum != 0]
      if (!is.vector(alpha0) && length(alpha0) != d) {
        stop("Error: dimension of alpha0 does not match data")
      } else if (any(alpha0 <= 0)) {
        # user provided starting values
        stop("Error: starting values should be positive")
      }
    }
    
    # prepare some variables before looping
    alphaHat <- alpha0
    alphaSum <- sum(alphaHat)
    loglIter <- sum(ddirmult(data, alpha0, log = TRUE))
    # backtrack max iterations
    backtrackMaxiters <- 10
    
    ##----------------------------------------##
    ## The Newton loop
    if (maxiters == 1) iter <- 1
    else { 
      for (iter in 2:maxiters) {
        
        # score vector
        alphaMat <- matrix(alphaHat, nrow = N, ncol = d, byrow = TRUE)
        score <- colSums(digamma(data + alphaMat)) - 
          N * (digamma(alphaHat) - digamma(alphaSum)) -
          sum(digamma(alphaSum + m))
        # observed info. matrix = diag(obsinfoDvec) - obsinfoC
        obsinfoDvec <- N * trigamma(alphaHat) - 
          colSums(trigamma(data + alphaMat))
        obsinfoDvecInv <- 1 / obsinfoDvec
        obsinfoC <- N * trigamma(alphaSum) - 
          sum(trigamma(m + alphaSum))
        # shrink c if necessary to make obs. info. pos def
        if (obsinfoC * sum(obsinfoDvecInv) >= 1) {
          obsinfoC <- 0.95 / sum(obsinfoDvecInv)
        }
        
        # compute Newton direction
        newtondir <- obsinfoDvecInv * score
        newtondir <- newtondir + 
       (sum(newtondir) / (1 / obsinfoC - sum(obsinfoDvecInv))) * obsinfoDvecInv
        
        # line search by step halving
        if (any(newtondir < 0)) {
          # make sure Newton iterate always lands within boundary
          stepsize <- min(- alphaHat[newtondir < 0] / newtondir[newtondir < 0])
          stepsize <- min(0.95 * stepsize, 1)
        } else {
          stepsize <- 1
        }
        for (btiter in 1:backtrackMaxiters) {
          alphaNew <- alphaHat + stepsize * newtondir
          loglNew <- sum(ddirmult(data, alphaNew, log = TRUE))
          # line search successful if improving log-L
          if (loglNew > loglIter) break
          else if (btiter == backtrackMaxiters) {
            warning("line search failed")
          } else {
            stepsize <- stepsize / 2
          }
        }
        alphaHat <- alphaNew
        alphaSum <- sum(alphaHat)
        loglOld <- loglIter
        loglIter <- loglNew
        
        # check convergence criterion
        if (abs(loglIter - loglOld) < tolfun * (abs(loglOld) + 1)) break
      } 
    }
    ##----------------------------------------##
    ## End of Newton loop
    
    # score, i.e., gradient
    alphaMat <- matrix(alphaHat, nrow = N, ncol = d, byrow = TRUE)
    score <- 
      colSums(digamma(data + alphaMat)) - 
      N * (digamma(alphaHat) - digamma(alphaSum)) -
      sum(digamma(alphaSum + m))
    
    # restore to original data size
    if (any(csum == 0)) {
      colidx <- (csum != 0)
      # parameter estimate
      tmp <- alphaHat
      alphaHat <- rep(0, ncol(X))
      alphaHat[colidx] <- tmp
    }
    
    # output
    return(alphaHat)
  }
  
  ##################################################
  
  MC_MLE <- function(I, d, S , N = 20) {
    #this is a function to calculate estimator for P using MLE method
    # Args: I: number of populations
    #       d: number of categories
    #       S: number of replicates
    #       N: batch size
    # Output Err MLE: estimated error for MLE estimator #
  set.seed(100)
  Err_MLE <- matrix(0, S) 
  temp <- matrix(0, I, 1)
    for ( s in 1:S) {    
      p <- matrix(runif(d * I), I, d)
      temp <-  rowSums(p)
      p <- p / temp
      data <- rmultinomial(I, N, p)  
      pr_MLE <- data / N
      Err_MLE[s] <- sum(abs(pr_MLE-p))/(2*I)
  }
    return(Err_MLE)
  }
  
  ##################################################
  
  MC_EB <- function(I, d, S, N = 20) {
    #
    # This is a function to calculate estimator for P using EB method
    # Args: I: number of populations
    #       d: number of categories
    #       S: number of replicates
    #       N: batch size
    # Output Err EB: estimated error for EB estimator 
    #
    set.seed(100)
    # initialize matrix
    data <- matrix(0, d ,I)
    temp1 <- matrix(0, I, 1)
    temp <- matrix(0, S, 1)
    # loop over different replicate
    for(s in 1:S){
      p <- matrix(runif(d * I), I, d)
      temp1 <-  rowSums(p)
      p <- p / temp1
      data <- rmultinomial(I, N, p)   
      # Newton method to estimate alpha
    alpha <- dirmultfit(data)
    pr_EB <- t((data + alpha) / (N + sum(t(alpha))))
    temp[s] <- sum(abs(pr_EB - t(p)))/(2 * I)
    }
    Err_EB <- temp
    return(Err_EB)
  }
  
  ##################################################
  
  MLE <- function(S, I, d) {
    # This is a function to calculate estimator error for different I and d
    # Args: I: number of populations (vector)
    #       d: number of categories (vector)
    #       S: number of replicates
    # Output Err MLE: estimated error for MLE estimator for different replicate 
    #
    tabledata_MLE <- matrix(0, S, length(I) * length(d))
  
    k = 1
    for (i in 1:length(I)){
      for (j in 1:length(d)) {
        tabledata_MLE[, k] <- MC_MLE(I[i], d[j], S)
        k <- k + 1
      }
    }
    return(tabledata_MLE)
  }
  
  ##################################################
  
  EB <- function(S, I, d){ 
    # This is a function to calculate estimator error for different I and d
    # Args: I: number of populations (vector)
    #       d: number of categories (vector)
    #       S: number of replicates
    # Output Err EB: estimated error for EB estimator for different replicate 
    #
    tabledata_EB <- matrix(0, S, length(I) * length(d))
    k <- 1
    for (i in 1:length(I)){
      for (j in 1:length(d)) {
        tabledata_EB[, k] <- MC_EB(I[i], d[j], S)
        k <- k + 1
      }
    }
    return(tabledata_EB)
  }
  
  # start simulation
  I <- c(10, 20, 50)
  d <- c(2, 3, 5, 10)
  S <- 100
  tabledata_MLE <- MLE(S, I, d)
  tabledata_EB <- EB(S, I, d)
  # Compute stardard error of estimated error
  SE_MLE <- matrix(apply(tabledata_MLE, 2, var), 4, 3)/S 
  SE_EB <- matrix(apply(tabledata_EB, 2, var), 4, 3)/S
  SE_MLE <- round(sqrt(SE_MLE), 3)
  SE_EB <- round(sqrt(SE_EB), 3)
  # Genearate table
  rownames(SE_EB) <- c("d = 2, SE", "d = 3, SE", "d = 5, SE", "d = 10, SE")
  rownames(SE_MLE) <- c("d = 2, SE", "d = 3, SE", "d = 5, SE", "d = 10, SE") 
  row.names(SE_EB) <- paste0(row.names(SE_EB), "(EB)")
  row.names(SE_MLE) <- paste0(row.names(SE_MLE), "(MLE)") 
  table_Mean_MLE <- colMeans(tabledata_MLE)
  table_Mean_MLE <- round(matrix(table_Mean_MLE, 4, 3), 2)
  ￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼￼
  colnames(table_Mean_MLE) <- c("I=10", "I=20", "I=50")
  rownames(table_Mean_MLE) <- c("d=2", "d=3", "d=5", "d=10")
  table_Mean_EB <- colMeans(tabledata_EB)
  table_Mean_EB <- round(matrix(table_Mean_EB, 4, 3), 2)
  colnames(table_Mean_EB) <- c("I=10", "I=20", "I=50")
  rownames(table_Mean_EB) <- c("d=2", "d=3", "d=5", "d=10") 
  row.names(table_Mean_EB) <- paste0(row.names(table_Mean_EB), "(EB)") 
  row.names(table_Mean_MLE) <- paste0(row.names(table_Mean_MLE), "(MLE)") 
  tableNew <- rbind(table_Mean_EB, SE_EB, table_Mean_MLE, SE_MLE) 
  tableNew1 <- tableNew[c(matrix(1:nrow(tableNew), nrow = 4, byrow = TRUE)), ]
                                                                                                                                                                                                                                                                                                 
library(knitr)
kable(tableNew1, format = 'pandoc')

library(glmnet)
library(MASS)
  
  
# Question 2
  
  # Generating data
genData <- function(n, beta, sigma) {
  # This is a function to generate data with size n, 
  #  beta and noise level of # sigma
  # Args: n : data size
  #     beta: true parameter
  # sigma: noise level
  # Outputs: list contains X and Y
  p <- length(beta)
  X <- matrix(rnorm(n*p), n, p)
  X[, 1] <- rep(1, n)
  y <- rnorm(n,X %*% beta, sigma ^ 2)
  return(list(X = X, y = y))
}

ridge <- function(data, beta) {
  # This is a function to find beta using ridge regression
  # Args: data: X and Y
  #       beta: true parameter
  # Outputs: list contains estimated coefficients
  # prediction error and estimation error
  #
  require(glmnet)
  cvfit <- cv.glmnet(data$X, data$y, family = "gaussian", alpha = 0)
  n <- length(data$y)
  coeff <- as.numeric(coef(cvfit, s=cvfit$lambda.min))[-2]
  residual <- predict(cvfit, newx = data$X)- data$y
  prederr <- sqrt(crossprod(residual) / n)
  return(list("coefficient" = coeff, "PredictionErr" = prederr, 
              "Err" =  sqrt(sum((coeff-beta)^2)/length(beta))))
}

lasso <- function(data, beta) {
  # This is a function to find beta using ridge regression
  # Args: data: X and Y
  #       beta: true parameter
  # Outputs: list contains estimated coefficients
  # prediction error and estimation error
  #
  require(glmnet)
  cvfit <- cv.glmnet(data$X, data$y, family = "gaussian", alpha = 1)
  n <- length(data$y)
  coeff <- as.numeric(coef(cvfit, s = cvfit$lambda.min))[-2]
  residual <- predict(cvfit, newx = data$X)- data$y
  prederr <- sqrt(crossprod(residual) / n)
  return(list("coefficient" = coeff, "PredictionErr" = prederr, 
              "Err" = sqrt(sum((coeff-beta)^2))/length(beta)))
}

OLS <- function(data, beta) {
  # This is a function to find beta using OLS regression
  # Args: data: X and Y
  #       beta: true parameter
  # Outputs: list contains estimated coefficients
  # prediction error and estimation error
  #
  fit <- lm(y ~ X, data)
  coeff <- coef(fit)[-2]
  prederr <- sqrt(crossprod(fit$residuals) / length(data$y))
  return(list("coefficient" = coeff, "PredictionErr" = prederr, 
              "Err" = sqrt(sum((coeff - beta)^2) /length(beta))))
}

set.seed(100)
beta <- runif(10, -5, 5)

simulation1 <- function(n, sigma, S, beta, i, j) {
  # This is the function do simulation across different combination of
  # sigma and n give specific value of S and beta
  # Args: n : data size
  #       sigma: noise level
  #       S : number of replicates
  #       beta : true value
  #       i,j: indextosetseed
  # Outputs : list contains the estimation
  # Err matrix for different method and replicates,
  # the prediction err matrix for different method and replicates,
  # true beta
p <- length(beta)
results <- array(NA,dim=c(S, p, 3),
                 dimnames=list(1:S, 1:p, c("Lasso","Ridge","OLS")))
predictionErr <- Err <- matrix(NA, S, 3)
set.seed(100)
for (s in 1:S) {
  set.seed((i - 1) * i + (j - 1) * j + s * (s - 1))
  # generate data
  Data <- genData(n, beta, sigma)
  lassfit <- lasso(Data, beta)
  ridgefit <- ridge(Data, beta)
  OLSfit <- OLS(Data, beta)
  
  # store results
  results[s, , 1] <- lassfit$coefficient
  results[s, , 2] <- ridgefit$coefficient
  results[s, , 3] <- OLSfit$coefficient
  predictionErr[s, 1] <- lassfit$PredictionErr
  predictionErr[s, 2] <- ridgefit$PredictionErr
  predictionErr[s, 3] <- OLSfit$PredictionErr
  Err[s, 1] <- lassfit$Err
  Err[s, 2] <- ridgefit$Err
  Err[s, 3] <- OLSfit$Err
}

return(list("Err" = Err, "PredictionErr" = predictionErr, "estimate" = results, 
            "true" = beta))  # return(Err_s vector)
}

S <- 700 # S = 700 can make error and prediction error significant  
# different settings in our study
  sigma <- c(sqrt(1), sqrt(0.8), sqrt(0.5), sqrt(0.1))
  N <- c(50, 100, 150, 200)
  k <- 1
  lasso_Err <- matrix(0, S, length(N) * length(sigma))
  ridge_Err <- lasso_Err
  OLS_Err <- lasso_Err
  lasso_pdErr <- matrix(0, S, length(N) * length(sigma))
  ridge_pdErr <- lasso_pdErr
  OLS_pdErr <- lasso_pdErr
  
  for ( sig in 1:length(sigma)) {
    for (Num in 1:length(N)) {
      sim <- simulation1(N[Num], sigma[sig], S, beta, sig, Num)
      lasso_Err[, k] <- sim$Err[, 1]
      ridge_Err[, k] <- sim$Err[, 2]
      OLS_Err[, k] <- sim$Err[, 3]
      lasso_pdErr[, k] <- sim$PredictionErr[, 1]
      ridge_pdErr[, k] <- sim$PredictionErr[, 2]
      OLS_pdErr[, k] <- sim$PredictionErr[, 3]
      k <- k +1
    }
  }
  
  table1_se_lassoErr <- matrix(apply(lasso_Err, 2, sd)/ sqrt(S), 4, 4)
  table2_se_ridgeErr <- matrix(apply(ridge_Err, 2, sd)/ sqrt(S), 4, 4)
  table3_se_OLSErr <- matrix(apply(OLS_Err, 2, sd)/ sqrt(S), 4, 4)
  table4_se_lassopderr <- matrix(apply(lasso_pdErr, 2, sd)/ sqrt(S), 4, 4)
  table5_se_ridgepderr <- matrix(apply(ridge_pdErr, 2, sd)/ sqrt(S), 4, 4)
  table6_se_OLSpderr <- matrix(apply(OLS_pdErr, 2, sd)/ sqrt(S), 4, 4)
  
  table1_lassoErr <- matrix(colMeans(lasso_Err), 4, 4)
  table2_ridgeErr <- matrix(colMeans(ridge_Err), 4, 4)
  table3_OLSErr <- matrix(colMeans(OLS_Err), 4, 4)
  table4_lassoErr <- matrix(colMeans(lasso_pdErr), 4, 4)
  table5_ridgeErr <- matrix(colMeans(ridge_pdErr), 4, 4)
  table6_OLSErr <- matrix(colMeans(OLS_pdErr), 4, 4)
  
# clear up
rm(list = ls())



