#' Acquisition functions for training set selection
#' 
#' @description
#' This function implements acquisition functions proposed in this Master's thesis,
#' which are designed to rank and select training sets for genomic selection
#' based on expected improvement-type and confidence upper bound-type criteria.
#' 
#' @param folds    Integer. Number of folds for cross-validation splitting.
#' @param candi    Vector. The list of candidates individuals. 
#' @param Nc       Integer. Total numbers of candidate.
#' @param kinshipA Matrix. Additive kinship matrix. 
#' @param kinshipD Matrix (Option). Dominance matrix.
#' @param p_true   Vector. Phenotypic value.
#' @param gamma    Numeric. penalty for EI calculate (default = 1).
#' @param kappa  Numeric. Kappa for upper bound calculate (default = 1.96).
#' 
#' @details
#' The function performs the following steps for each fold:
#' 1. Splits candidates into a training set and remaining set. 
#' 2. Fits a GBLUP model using the \code{BGLR} package to estimate genetic effects.
#' 3. Calculates the conditional posterior mean vector and variance-covariance matrix for the remaining sets.
#' 4. Compute the current maximum genetic value ($f_{Mg}$) in the training set.
#' 5. Derives EI and UCB scores based on the specified \code{gamma} and \code{kappa}.
#' 
#' @return A data frame with \code{Nc} rows and 2 columns:
#' \item{meanEI}{The average Expectated Improvment (EI) score across folds.}
#' \item{meanUCB}{The average Confidence Upper Bound (UCB) score across folds. }
#' 
#' @import BGLR
#' @importFrom stats pnorm dnorm
#' 
#' @export

AcqFun <- function(folds, candi, Nc, kinshipA, kinshipD = NULL, p_true, gamma = 1, kappa = 1.96){
  
  #---  Check dominance effect ---#
  dominance <- !is.null(kinshipD)
  
  #--- Store data ---#
  final <- data.frame(meanEI = rep(NA_real_, Nc), meanUCB = rep(NA_real_, Nc))
  EImat <- matrix(NA_real_, nrow = Nc, ncol = folds)
  UCBmat <- matrix(NA_real_, nrow = Nc, ncol = folds)
  rownames(EImat) <- rownames(UCBmat) <- rownames(kinshipA)
  
  #--- Cluster ---#
  part <- split(sample(candi), rep(seq_len(folds), length.out = Nc))
  
  #--- Main structure ---# 
  for (i in seq_len(folds)) {
    train <- part[[i]]
    remain <- setdiff(candi, train)
    p <- p_true[train]
    
    #--- sub-matrices for additive ---#
    ka11 <- as.matrix(kinshipA[train, train, drop = FALSE])
    ka12 <- as.matrix(kinshipA[train, remain, drop = FALSE])
    ka21 <- as.matrix(kinshipA[remain, train, drop = FALSE])
    ka22 <- as.matrix(kinshipA[remain, remain, drop = FALSE])
    
    #--- Construct ETA for BGLR ---#
    ETA_list <- list(A = list(K = ka11, model = "RKHS"))
    
    #--- sub-matrices for dominance ---#
    if(dominance){
      kd11 <- as.matrix(kinshipD[train, train, drop = FALSE])
      kd12 <- as.matrix(kinshipD[train, remain, drop = FALSE])
      kd21 <- as.matrix(kinshipD[remain, train, drop = FALSE])
      kd22 <- as.matrix(kinshipD[remain, remain, drop = FALSE])
      
      ETA_list$D <- list(K = kd11, model = "RKHS")
    }
    
    #--- BGLR GBLUP Model ---#
    fit_BGLR <- suppressWarnings(BGLR(y = p, 
                                      ETA = ETA_list, 
                                      verbose = FALSE, nIter = 10000, burnIn = 1000))
    varE <- fit_BGLR$varE
    mu <- fit_BGLR$mu
    
    #--- Additive Effect ---#
    varA <- fit_BGLR$ETA$A$varU
    g1A <- fit_BGLR$ETA$A$u; g1_total <- g1A

    inv_ka11 <- tryCatch({
      chol2inv(chol(ka11+diag(1e-6, nrow(ka11))))
    }, error = function(e){
      solve(ka11+diag(1e-6, nrow(ka11)))
    })
    mu_hat_g2 <- ka21%*%inv_ka11%*%g1A
    var_hat_g2 <- varA*diag(ka22-ka21%*%inv_ka11%*%ka12)
    var_total <- varA*diag(ka11)
    
    #--- Dominance Effect exist ---#
    if(dominance){
      varD <- fit_BGLR$ETA$D$varU
      g1D <- fit_BGLR$ETA$D$u
      
      g1_total <- g1_total+g1D
      
      inv_kd11 <- tryCatch({
        chol2inv(chol(kd11+diag(1e-6, nrow(kd11))))
      }, error = function(e){
        solve(kd11+diag(1e-6, nrow(kd11)))
      })
      mu_hat_g2_D <- kd21%*%inv_kd11%*%g1D
      var_hat_g2_D <- varD*diag(kd22-kd21%*%inv_kd11%*%kd12)
      
      mu_hat_g2 <- as.numeric(mu_hat_g2+mu_hat_g2_D)
      var_hat_g2 <- as.numeric(var_hat_g2+var_hat_g2_D)
      
      var_total <- var_total + varD*diag(kd11)
    }
    
    #--- Posterior for remaining set ---#
    mu_hat_g2 <- as.numeric(mu_hat_g2)
    var_hat_g2 <- as.numeric(var_hat_g2)
    
    #--- f*_Mg ---#
    uf <- g1_total - gamma*sqrt(var_total)
    fMg <- max(uf)
    
    #--- Expected Improvement ---#
    if(any(var_hat_g2 < 0)) var_hat_g2[var_hat_g2 < 0] <- 1e-8
    
    z <- (mu_hat_g2 - fMg) / sqrt(var_hat_g2)
    aug_ei <- ((mu_hat_g2 - fMg)*pnorm(z) + sqrt(var_hat_g2)*dnorm(z))* 
              (1 - sqrt(varE) / sqrt(var_hat_g2 + varE))
    
    #--- Confidence Upper Bound with Augmented ---#
    aug_ucb <- (mu_hat_g2 + kappa*sqrt(var_hat_g2))*
      (1 - sqrt(varE) / sqrt(var_hat_g2 + varE))
    
    #--- Store Result ---#
    EImat[remain, i] <- aug_ei
    UCBmat[remain, i] <- aug_ucb
    
  }
  
  #--- Average ---#
  final[, 1] <- rowMeans(EImat, na.rm = TRUE)
  final[, 2] <- rowMeans(UCBmat, na.rm = TRUE)
  
  #--- Return ---#
  return(final)
}
