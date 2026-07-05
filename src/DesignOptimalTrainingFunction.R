#' Generating the Optimal training set
#' 
#' @description
#' This function performs a Monte Carlo-based Bayesian Optimization to select 
#' the optimal training set. It simulates true breeding value based on kinship 
#' matrices to calculated acquisition function values (EI and UCB) for each 
#' candidate individuals.
#' 
#' @param kinshipA Matrix. Additive kinship matrix.
#' @param kinshipD Matrix (Option). Dominance kinship matrix.
#' @param nsim     Integer. Number of Monte Carlo simulations (default = 2500).
#' @param folds    Integer. Numbers of folds for cross-validation.
#' @param h        Numeric. The heritability (default = 0.5).
#' @param mu       Numeric. The overall population mean (default = 100).
#' @param varA     Numeric. The additive genetic variance component (default = 20).
#' @param rho      Numeric. The ratio of dominance variance to additive variance (varD/varA).
#'                          Used to calculate varD = rho * varA.
#' @param gamma    Numeric. A tuning parameter for the EI function (default = 1).
#' @param z.score  Numeric. The Z-score used for the Upper Confidence Bound (UCB) calculation 
#'                 (default = 1.96, corresponding to a 95% confidence interval).
#'                 
#' @return A list containing three elements:
#' \item{EI}{Data Frame. The ranking results based on the Expectated Improvement (EI) criterion. }
#' \item{UCB}{Data Frame. The ranking results based on the Confidence Upper Bound (UCB) criterion. }
#' \item{GV}{Data Frame. The A-optimality-like criterion. }
#' 
#' @import BGLR 
#' @import MASS
#' @import future
#' @import future.apply
#' @import progressr
#' @import cli
#' 
#' @export

OPTtrain <- function(kinshipA, kinshipD = NULL, nsim = 2500, folds = 5, h = 0.5, mu = 100, 
                     varA = 20, rho, gamma = 1, z.score = 1.96){
  
  suppressPackageStartupMessages({
    library(BGLR)
    library(MASS)
    requireNamespace("progressr", quietly = TRUE)
    requireNamespace("cli", quietly = TRUE)
    requireNamespace("future", quietly = TRUE)
    requireNamespace("future.apply", quietly = TRUE)
  })
  
  #--- Progress Bar Setting ---#
  progressr::handlers(progressr::handler_progress(
    format = " Processing [:bar] :percent | ETA: :eta | Time: :elapsed", 
    complete = "=", 
    incomplete = " ", 
    current = ">", 
    width = 100, 
    clear = FALSE
  ))
  
  #--- Basic Check ---#
  options(progressr.interval = 0.5, future.globals.maxSize = +Inf)
  cli::cli_h1("Starting Bayesian Optimization Approach")
  
  if(!is.matrix(kinshipA)) stop("{.var kinshipA} must be a matrix. ")
  
  if(nrow(kinshipA) != ncol(kinshipA)) stop("{.var kinshipA} must be square matrix. ")
  
  if(!is.null(kinshipD)){
    if(!is.matrix(kinshipD)) stop("{.var kinshipD} must be a matrix.")
    if(!identical(dim(kinshipA), dim(kinshipD))) stop("Dimensions of {.var kinshipA} and 
                                                      {.var kinshipD} must match.")
  }
  
  #--- Pre-process ---#
  Nc <- nrow(kinshipA)
  candi <- seq_len(Nc)
  
  #--- Variances & Phenotype simulation---#
  cli::cli_alert_info("Simulating phenotype data for {.val {nsim}} iterations...")
  
  if(!is.null(kinshipD)){
    varD <- rho * varA
    gD <- tryCatch(mvrnorm(nsim, mu = rep(0, Nc), Sigma = varD*kinshipD), 
                   error = function(e) stop("Error in generating gD: Matrix might not be positive definite."))
  } else {
    varD <- 0
    gD <- matrix(0, nrow = nsim, ncol = Nc)
  }
  
  varE <- (varA+varD)*(1-h)/h
  
  gA <- mvrnorm(nsim, mu = rep(0, Nc), Sigma = varA*kinshipA) # additive effect
  e <- mvrnorm(nsim, mu = rep(0, Nc), Sigma = varE*diag(Nc))
  p_new <- mu+gA+gD+e
  
  #--- Parallel Computing ---#
  ncore <- max(1, parallel::detectCores() - 2) # parallel core setting
  future::plan(future::multisession, workers = ncore)
  
  cli::cli_text(cli::col_yellow("Processing Cross-Validation Loop..."))
  
  #--- CV loop (nsim) ---#
  result <- progressr::with_progress({
    p <- progressr::progressor(steps = nsim)
    
    future.apply::future_lapply(seq_len(nsim), function(i){
      res <- AcqFun(folds = folds, candi = candi, Nc = Nc, kinshipA = kinshipA, kinshipD = kinshipD, 
                    p_true = p_new[i, ], gamma = gamma, z.score = z.score)
      
      p()
      return(res)
    }, future.seed = TRUE, 
    future.packages = c("BGLR", "MASS"), 
    future.globals = TRUE)
  })
  future::plan(future::sequential)
  
  cli::cli_alert_success("Cross-Validation Processing Finished.")
  
  #--- Summarize for Bayesian Optimization ---#
  
  if(is.null(rownames(kinshipA))){
    name_list <- as.character(seq_len(nrow(kinshipA)))
  } else {
    name_list <- rownames(kinshipA)  
  }
  
  mat_meanEI <- vapply(result, function(y) as.numeric(y$meanEI), numeric(Nc))
  mat_meanUCB <- vapply(result, function(y) as.numeric(y$meanUCB), numeric(Nc))
  rownames(mat_meanEI) <- rownames(mat_meanUCB) <- rownames(kinshipA)
  
  #--- Average 
  augEI <- rowMeans(mat_meanEI, na.rm = TRUE)
  ucb <- rowMeans(mat_meanUCB, na.rm = TRUE)
  
  #--- Ranking 
  ord_EI <- order(augEI, decreasing = TRUE)
  ord_ucb <- order(ucb, decreasing = TRUE)
  
  #--- Result
  res_EI <- data.frame(
    speciesName = rownames(kinshipA)[ord_EI], 
    aug.EI = augEI[ord_EI], 
    posIndex = ord_EI
  )
  
  res_UCB <- data.frame(
    speciesName = rownames(kinshipA)[ord_ucb], 
    UCB = ucb[ord_ucb], 
    posIndex = ord_ucb
  )
  
  #--- Summarize for A-optimality-like ---#
  res_GV <- GVave(kinshipA = kinshipA, kinshipD = kinshipD)
  
  #--- Return ---#
  cli::cli_alert_success("Optimization Done! ")
  cli::cli_rule()
  
  return(list(EI = res_EI, UCB = res_UCB, GV = res_GV))
}
