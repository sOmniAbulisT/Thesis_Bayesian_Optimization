#' Calculate Normalized Discounted Cumulative Gain (NDCG)
#'
#' @description
#' Evaluates the ranking quality of the predicted breeding values. NDCG measures the 
#' usefulness (gain) of a selected top-k set based on their positions in the ranked list.
#' The gain is accumulated from the top of the result list to the bottom, with the 
#' gain of each result discounted at lower ranks.
#'
#' @details
#' This function assumes that 'TBV' (True Breeding Value) acts as the relevance score.
#' It calculates the ratio of the Discounted Cumulative Gain (DCG) of the predicted 
#' order to the Ideal DCG (IDCG) of the perfect order.
#'
#' @param GEBV Numeric Vector. Genomic Estimated Breeding Values (Predictions).
#' @param TBV  Numeric Vector. True Breeding Values (True Labels/Relevance).
#' @param k    Integer. The cutoff rank (top-k) to evaluate.
#'
#' @return Numeric. The NDCG score ranging from 0 to 1 (assuming non-negative TBVs).
#'         A higher value indicates better ranking performance.
#' @export

NDCG <- function(GEBV, TBV, k){
  
  #--- Basic check ---#
  gain <- match.arg(gain)
  stopifnot(length(GEBV) == length(TBV))
  stopifnot(k <= length(GEBV) && k <= length(TBV))
  
  #--- Discounted function ---#
  d <- 1 / log2(seq_len(k) + 1)
  
  #--- DCG & IDCG ---#
  index_pred <- order(GEBV, decreasing = TRUE)[seq_len(k)]
  index_true <- order(TBV, decreasing = TRUE)[seq_len(k)]
  
  DCG <- sum(TBV[index_pred] * d)
  IDCG <- sum(TBV[index_true] * d)

  ndcg <- if(IDCG == 0) 0 else DCG / IDCG
  
  #--- Return ---#
  return(ndcg)
}

#' Calculate Rank Sum Ratio (RSR)
#'
#' @description
#' A metric to evaluate if the selected top-k individuals truly belong to the top tier.
#' It compares the sum of the ideal ranks (1 to k) with the sum of the true ranks 
#' of the individuals selected by the model.
#'
#' @details
#' RSR = (Sum of Ideal Ranks) / (Sum of True Ranks of Selected Individuals).
#' \itemize{
#'   \item An RSR of 1 indicates perfect selection (the selected top-k are exactly the true top-k).
#'   \item Lower values indicate that the selected individuals have poor true rankings.
#' }
#'
#' @param GEBV Numeric Vector. Genomic Estimated Breeding Values.
#' @param TBV  Numeric Vector. True Breeding Values.
#' @param k    Integer. The number of top individuals selected.
#'
#' @return Numeric. The Rank Sum Ratio.
#' @export

RSR <- function(GEBV, TBV, k) {
  #--- Basic check ---#
  stopifnot(length(GEBV) == length(TBV))
  stopifnot(k > 0, k <= length(TBV))
  
  #--- ideal rank sum ---#
  true_rs <- sum(seq_len(k))
  
  #--- pred. rank sum ---#
  pred_index <- order(GEBV, decreasing = TRUE)
  rs <- sum(rank(-TBV[pred_index])[seq_len(k)])
  
  #--- rank sum ratio ---#
  rsr <- true_rs / rs
  
  #--- return ---#
  return(rsr)
}

#' Calculate Spearman's Rank Correlation for Top-k (SRC)
#'
#' @description
#' Calculates the rank correlation coefficient specifically for the top-k individuals 
#' selected by the prediction model. This assesses how well the model preserves the 
#' relative order among the best-performing individuals.
#'
#' @details
#' Unlike standard Spearman correlation which considers the whole population, 
#' this function:
#' 1. Selects the top-k individuals based on GEBV.
#' 2. Ranks them locally based on GEBV.
#' 3. Ranks them locally based on TBV.
#' 4. Computes the correlation between these two sets of local ranks.
#'
#' @param GEBV Numeric Vector. Genomic Estimated Breeding Values.
#' @param TBV  Numeric Vector. True Breeding Values.
#' @param k    Integer. The number of top individuals to consider.
#'
#' @return Numeric. The correlation coefficient (between -1 and 1).
#' @export

SRC <- function(GEBV, TBV, k){
  
  #--- Basic check ---#
  stopifnot(length(GEBV) == length(TBV))
  stopifnot(k > 0, k <= length(TBV))
  
  index <- order(GEBV, decreasing = TRUE)[seq_len(k)]
  
  #--- Simulate breeding values and Predictive breeding values ---#
  bp_top_k <- GEBV[index]
  true_top_k <- TBV[index]
  
  #--- Ranking ---#
  Rx <- rank(-bp_top_k, na.last = NA)
  Ry <- rank(-true_top_k, na.last = NA)
  
  #--- Covariance ---#
  numerator <- sum((Rx - (k+1)/2) * (Ry - mean(Ry)))
  
  #--- Variance ---#
  denominator <- sqrt((sum((Rx - (k+1)/2)^2)) * (sum((Ry - mean(Ry))^2)))
  
  #--- Ranking correlation ---#
  if(denominator == 0){
    src <- 0
  }else{
    src <- numerator / denominator
  }
  
  #--- Return ---#
  return(src)
}

#' Generate true genotypic values and Simulated phenotypic values
#' 
#' @description
#' This function performs a data-driven simulation in two steps:
#' 1. Estimates variance components from real data using BGLR (RKHS model).
#' 2. Simulates new True Breeding Values (TBVs) and phenotypes based on the estimated parameters.
#' 
#' @param kinshipA Matrix. The additive kinship matrix (Nc * Nc).
#' @param kinshipD Matrix (option). The dominance kinship matrix (Nc * Nc).
#' @param yReal    Vector. The real phenotypic data used for parameter estimation.
#' @param dataSet  String. A prefix for the output filename (e.g. "Wheat").
#' @param nsim     Integer. Number of simulation replicates. Defaults to 2000.
#' @param output_path String. Suffix for the output filename. Defaults to "TrueValue.RData".
#' @param trait_name String (Optional). Name of the trait (e.g., "GrainYield"). Used for filename.
#' 
#' @return Invisible returns the file path. Saves an .RData file containing simulated matrices.
#' 
#' @import BGLR
#' @import MASS

GenerateTBV <- function(kinshipA, kinshipD = NULL, yReal, dataSet, nsim = 2000, 
                        trait = NULL){
  
  suppressPackageStartupMessages({
    library(BGLR)
    library(MASS)
    requireNamespace("cli", quietly = TRUE)
  })
  
  Nc <- nrow(kinshipA)
  
  if(any(is.na(yReal))){
    cat(sprintf("Warning: yReal contains %d NAs. BGLR will handle them.\n", sum(is.na(yReal))))
  }
  
  #--- Pre-process ---#
  cli::cli_h1("[0/3]: Data Preparation")
  trait_name <- ifelse(is.null(trait), "Unknown", trait)
  
  cli::cli_text("Current Task Info:")
  cli::cli_ul(c(
    "Dataset: {.val {dataSet}}",
    "Trait  : {.val {trait_name}}"
  ))
  
  #--- Estimated parameter ---#
  cli::cli_h1("[1/3]: Estimating Parameters...")
  
  ETA_list <- list(A = list(K = kinshipA, model = "RKHS")) # Additive
  
  if(!is.null(kinshipD)){
    ETA_list$D <- list(K = kinshipD, model = "RKHS") # Dominance
    model_type <- "Additive + Dominance model"
  }else{
    model_type <- "Additive model"
  }
  
  cli::cli_alert_info("Model detected: {.field {model_type}}")
  
  tmp_file <- file.path(tempdir(), paste0("BGLR_", Sys.getpid(), "_", sample(1e4, 1), "_"))
  fit_BGLR <- BGLR(y = yReal, 
                   ETA = ETA_list, 
                   nIter = 10000, 
                   burnIn = 1000, 
                   verbose = FALSE, 
                   saveAt = tmp_file)
  
  mu_hat <- fit_BGLR$mu
  gA_hat <- as.vector(fit_BGLR$ETA$A$u)
  gD_hat <- if(!is.null(kinshipD)) as.vector(fit_BGLR$ETA$D$u) else rep(0, Nc)
  
  #true_mean <- mu_hat+gA_hat+gD_hat
  
  VarA <- fit_BGLR$ETA$A$varU
  VarE <- fit_BGLR$varE
  VarD <- if(!is.null(kinshipD)) fit_BGLR$ETA$D$varU else 0
  
  unlink(list.files(tempdir(), pattern = basename(tmp_file), full.names = TRUE))
  
  cli::cli_text("Estimated mean: ")
  cli::cli_ul("Mean: {.val {round(mu_hat, 4)}}")
  
  cli::cli_text("Estimated Variance Components: ")
  cli::cli_ul(c(
    "Variance A: {.val {round(VarA, 4)}}",
    if(!is.null(kinshipD)) "Variance D: {.val {round(VarD, 4)}}",
    "Variance E: {.val {round(VarE, 4)}}"
  ))
  cli::cli_text("Heritability (h2): {.val {round((VarA+VarD)/(VarA+VarD+VarE), 4)}}")
  
  #--- Simulated true breeding value ---#
  cli::cli_h1("[2/3]: Simulation True Breeding Values")
  
  # Additive
  A <- cli::cli_progress_step("Simulating Additive Effects")
  gA_true <- mvrnorm(n = nsim, mu = rep(0, Nc), Sigma = VarA*kinshipA)
  cli::cli_progress_done(A)
  
  # Dominance (option)
  if(!is.null(kinshipD)){
    D <- cli::cli_progress_step("Simulating Dominance Effects")
    gD_true <- mvrnorm(n = nsim, mu = rep(0, Nc), Sigma = VarD*kinshipD)
    cli::cli_progress_done(D)
  } else {
    gD_true <- matrix(0, nrow = nsim, ncol = Nc)
  }
  
  # Residual 
  E <- cli::cli_progress_step("Simulating Residuals")
  e <- mvrnorm(n = nsim, mu = rep(0, Nc), Sigma = VarE*diag(Nc))
  cli::cli_progress_done(E)
  
  #--- Finalizing & Saving ---#
  cli::cli_h1("[3/3]: Saving Simulation Results")
  
  cli::cli_alert_success("Simulation Data Ready! ")
  
  result <- list(
    mu = mu_hat,
    gA_true = gA_true, 
    gD_true = gD_true, 
    residual = e, 
    params = list(VarA = VarA, VarD = VarD, VarE = VarE)
  )
  
  return(result)
}

#' Function of Evaluation Metrics
#' 
#' @param kinshipA Matrix. Additive kinship matrix (Nc * Nc).
#' @param kinshipD Matrix (option). Dominance kinship matrix (Nc * Nc). 
#' @param OptTRS   Dataframe. Optimization result.
#' @param TRSSize  Integer. Size of training set sizes.
#' @param p_true   Vector. Phenotypic values. (from simulation data)
#' @param mu_true  Numeric. Estimated fix effect.
#' @param g_true   Vector. Genotypic values. (from simulation data)
#' @param k        Integer. Top-k for ranking evaluation metrices.
#' 
#' @import BGLR
#' @import tibble
#' 
Evaluation <- function(kinshipA, kinshipD = NULL, OptTRS, TRSSize, p_true, mu_true, g_true, k){
  
  #--- Set selection strategy ---#
  TRS_index <- OptTRS$posIndex[seq_len(TRSSize)]
  
  #--- sub-matrices (training set) ---#
  kAtt <- as.matrix(kinshipA[TRS_index, TRS_index, drop = FALSE])
  y_train <- as.numeric(p_true[TRS_index])
  
  ETA_list <- list(A = list(K = kAtt, model = "RKHS"))
  if(!is.null(kinshipD)){
    kDtt <- as.matrix(kinshipD[TRS_index, TRS_index, drop = FALSE])
    ETA_list$D <- list(K = kDtt, model = "RKHS")
  }
  
  tmp_prefix <- tempfile(pattern = paste0("BGLR_PID", Sys.getpid(), "_"))
  
  on.exit({
    files <- list.files(dirname(tmp_prefix), pattern = basename(tmp_prefix), full.names = TRUE)
    if(length(files)>0) unlink(files)
  }, add = TRUE)
  
  #--- model building and predict ---#
  fit_BGLR <- tryCatch({
    BGLR(y = y_train, 
         ETA = ETA_list, 
         verbose = FALSE, nIter = 5000, burnIn = 1000, 
         saveAt = tmp_prefix)
  }, error = function(e){
    return(NULL)
  })
  if(is.null(fit_BGLR)) {
    return(tibble::tibble(metric = c("NDCG", "RSR", "SRC"), value = NA_real_))
  }
  
  #--- BLUPs on training set ---#
  #--- Additive
  g_hat_A <- as.numeric(fit_BGLR$ETA$A$u)
  kAct <- as.matrix(kinshipA[, TRS_index, drop = FALSE])
  
  inv_kAtt <- tryCatch({
    chol2inv(chol(kAtt + diag(1e-6, nrow(kAtt))))
  }, error = function(e) {
    solve(kAtt + diag(1e-6, nrow(kAtt)))
  })
  
  #--- additive effect of whole candidate individuals
  pred_A <- kAct%*%inv_kAtt%*%g_hat_A 
  
  #--- Dominance
  if(!is.null(kinshipD)){
    g_hat_D <- as.numeric(fit_BGLR$ETA$D$u)
    kDct <- as.matrix(kinshipD[, TRS_index, drop = FALSE]) 
    
    inv_kDtt <- tryCatch({
      chol2inv(chol(kDtt + diag(1e-6, nrow(kDtt))))
    }, error = function(e){
      solve(kDtt + diag(1e-6, nrow(kDtt)))
    })
    
    #--- dominance effect of whole candidate individuals
    pred_D <- kDct%*%inv_kDtt%*%g_hat_D
  } else {
    pred_D <- 0
  }
  
  #--- Calculate GEBV of whole candidate population ---#
  g_pred <- pred_A + pred_D
  est_mu <- fit_BGLR$mu
  gebv <- est_mu + g_pred
  
  #--- TBV 
  tbv <- mu_true + g_true
  
  #--- Evaluation metrics ---#
  #--- NDCG
  ndcg <- NDCG(GEBV = gebv, TBV = tbv, k = k, gain = "linear")
  
  #--- Rank sum ratio
  rsr <- RSR(GEBV = gebv, TBV = tbv, k = k)
  
  #--- Spearman's ranking correlation
  src <- SRC(GEBV = gebv, TBV = tbv, k = k)
  
  #--- Summary metrics ---#
  result <- tibble::tibble(
    metric = c("NDCG", "RSR", "SRC"),
    value  = c(ndcg, rsr, src)
  )
  
  #--- Return results ---#
  return(result)
}

#' Simulate different heritability and training set size scenarios
#' 
#' @description
#' This function performs a comprehensive simulation study to validate the performance
#' of optimized training sets (OptTRS) against random selection. It uses a Data-Driven
#' approach where simulation parameters (variance componenets) are estimated  from 
#' real phenotypic data. 
#' 
#' @details
#' The simulation process involves:
#' 1. Estimating variance component (varA, varD and varE) from 'yReal'.
#' 2. Generating 'nsim' replicates of True Breeding Values and Phenotypic values. 
#' 3. Running parallel evaluations to compare EI, UCB, GV and Random strategies.
#' 4. Calculating ranking based metrics: NDCG, RS_ratio and RSR
#' 
#' @param kinshipA   Matrix. The additive kinship matrix. 
#' @param kinshipD   Matrix (Option). The dominance kinship matrix. 
#' @param yReal      Vector. The phenotypic values from real data uses to estimated parameters.
#' @param OptTRS     List. A list containing data frames for different strategies (EI, UCB, GV). 
#'                   Each data frame must contain a column 'posIndex'.
#' @param nsim       Integer. Number of simulation replicates (default = 2000). 
#' @param TRSSize    Integer Vector. A sequence of training set sizes to evaluate (e.g. c(50, 100)). 
#' @param k          Integer. The 'k' parameter for Top-k metrics. 
#' @param dataSet    String. Name of data set. 
#' @param trait      String. Name of trait. 
#' @param output_dir String. Directory path to save intermediate and final results.
#' 
#' @return A list containing:
#' \item {Optimal}{Data Frame. Summary statistics (Mean and SD) of metrics aggregated by strategies and size.}
#' 
#' @import BGLR
#' @import dplyr
#' @import tibble
#' @import purrr
#' @import MASS
#' @import progressr
#' @import cli
#' @import future
#' @import future.apply
#' 
#' @export

Simulation <- function(kinshipA, kinshipD = NULL, yReal, OptTRS, nsim = 2000, TRSSize, k, 
                       dataSet = "Unknown_DataSet", trait = "Unknown_Trait", 
                       output_dir = "Simulation_Results"){
  
  suppressPackageStartupMessages({
    library(BGLR)
    library(dplyr)
    library(tibble)
    library(purrr)
    library(MASS)
    requireNamespace("progressr", quietly = TRUE)
    requireNamespace("cli", quietly = TRUE)
    requireNamespace("future", quietly = TRUE)
    requireNamespace("future.apply", quietly = TRUE)
  })
  
  #--- progress bar setting ---#
  pb <- progressr::handler_progress(
    format = " Simulation Process [:bar] :percent | ETA: :eta", 
    complete = "=", 
    incomplete = " ", 
    current = ">", 
    width = 100, 
    clear = FALSE
  )
  
  options(future.globals.maxSize = +Inf)
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  Nc <- nrow(kinshipA)
  
  #--- Generate Data (Data-Driven) ---#
  cli::cli_h1(sprintf("Simulation Start: [%s] - %s", dataSet, trait))
  cli::cli_alert_info("Generating True Breeding Values...")
  
  Sim_TBV <- GenerateTBV(kinshipA = kinshipA, 
                         kinshipD = kinshipD, 
                         yReal = yReal, 
                         dataSet = dataSet, 
                         trait = trait, 
                         nsim = nsim)
  
  # parameters
  g_true <- Sim_TBV$gA_true+Sim_TBV$gD_true
  e_true <- Sim_TBV$residual
  mu_true <- Sim_TBV$mu
  
  p_true <- sweep(g_true, MARGIN = 2, STATS = mu_true, FUN = "+") + e_true # true breeding value
  
  # estimated heritability
  VarA <- Sim_TBV$params$VarA; VarE <- Sim_TBV$params$VarE
  VarD <- Sim_TBV$params$VarD
  h2 <- (VarA+VarD)/(VarA+VarD+VarE)
  
  cli::cli_alert_success(sprintf("Estimated Heritability: %.3f", h2))
  
  #--- Parallel Plan Setting ---#
  future::plan(future::multisession)
  
  #--- Simulation Loop ---#
  final <- list()
  
  progressr::handlers(pb)
  
  cli::cli_h1("Running Evaluation Loop")
  
  for(size in TRSSize){
    
    cli::cli_h2(sprintf("Processing Training Set Size: %d", size))
    
    res_size <- progressr::with_progress({
      p <- progressr::progressor(steps = nsim)
      
      res_list <- future.apply::future_lapply(seq_len(nsim), function(s){
        
        tryCatch({
          ps <- p_true[s, ]; gs <- g_true[s, ]
          
          ## Random Strategy
          random_idx <- sample(Nc, size)
          RandomTRS <- data.frame(posIndex = random_idx)
          
          ## Evaluations
          res_eva <- dplyr::bind_rows(
            # EI
            Evaluation(kinshipA = kinshipA, kinshipD = kinshipD, OptTRS = OptTRS$EI, TRSSize = size, 
                       p_true = ps, g_true = gs, mu_true = mu_true, k = k) |> 
              dplyr::mutate(strategy = "EI"), 
            
            # UCB
            Evaluation(kinshipA = kinshipA, kinshipD = kinshipD, OptTRS = OptTRS$UCB, TRSSize = size, 
                       p_true = ps, g_true = gs, mu_true = mu_true, k = k) |> 
              dplyr::mutate(strategy = "UCB"), 
            
            # GVaverage
            Evaluation(kinshipA = kinshipA, kinshipD = kinshipD, OptTRS = OptTRS$GV, TRSSize = size, 
                       p_true = ps, g_true = gs, mu_true = mu_true, k = k) |> 
              dplyr::mutate(strategy = "GVaverage"), 
            
            # Random baseline
            Evaluation(kinshipA = kinshipA, kinshipD = kinshipD, OptTRS = RandomTRS, TRSSize = size, 
                       p_true = ps, g_true = gs, mu_true = mu_true, k = k) |> 
              dplyr::mutate(strategy = "Random"), 
          )
          
          res_eva$TRSSize <- size
          p()
          return(res_eva)
          
        }, error = function(e){
          message(sprintf("\n[ERROR] Sim %d failed: %s", s, e$message))
          
          return(tibble::tibble(
            metric = c("NDCG", "RSR", "SRC"), 
            value = NA_real_, 
            strategy = "Error", 
            TRSSize = size
          ))
        })
        
      }, 
      future.seed = TRUE, 
      future.packages = c("BGLR", "MASS", "dplyr", "tibble"), 
      future.globals = TRUE)
      
      dplyr::bind_rows(res_list)
    })
    
    res_size <- res_size |> 
      dplyr::mutate(dataset = dataSet, trait = trait, heritability = h2, k = k, .before = 1)
    
    saveName <- file.path(output_dir, sprintf("Res_%s_%s_Size%d.RData", dataSet, trait, size))
    save(res_size, file = saveName)
    cli::cli_alert_success(sprintf("Saved results for size %d to %s", size, basename(saveName)))
    
    final[[as.character(size)]] <- res_size
  }
  future::plan(future::sequential)
  
  #--- Combine all heritability results ---#
  cli::cli_h1("Summarizing...")
  
  all_table <- dplyr::bind_rows(final)
  
  Optimal <- all_table |>  
    dplyr::group_by(dataset, trait, strategy, heritability, TRSSize, metric, k) |> 
    dplyr::summarize(
      Mean = mean(value, na.rm = TRUE), 
      SD = sd(value, na.rm = TRUE), 
      .groups = "drop"
    )
  
  cli::cli_alert_success(sprintf("Simulation Done! All results saved in '%s'", output_dir))
  cli::cli_rule()
  
  #--- Return ---#
  return(list(Optimal = Optimal, true_mean = mu_true))
}


#--- Simulation2 --- #
Simulation2 <- function(kinshipA, kinshipD = NULL, OptTRS, nsim = 2000, h = c(0.2, 0.5, 0.8), varA = 20, mu = 100,
                        gamma = c(0.5, 1, 2, 4), TRSSize, k){
  
  suppressPackageStartupMessages({
    library(BGLR)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(purrr)
    library(MASS)
    requireNamespace("progressr", quietly = TRUE)
    requireNamespace("cli", quietly = TRUE)
    requireNamespace("future", quietly = TRUE)
    requireNamespace("future.apply", quietly = TRUE)
  })
  
  total_steps <- length(gamma) * length(h) * nsim
  
  #--- progress bar setting ---#
  pb <- progressr::handler_progress(
    format = " Simulation Process [:bar] :percent | ETA: :eta", 
    complete = "=", 
    incomplete = " ", 
    current = ">", 
    width = 100, 
    clear = FALSE
  )
  
  options(future.globals.maxSize = +Inf)
  
  Nc <- nrow(kinshipA)
  
  gA_true <- mvrnorm(n = nsim, mu = rep(0, Nc), Sigma = varA*kinshipA)
  progressr::handlers(pb)
  future::plan(future::multisession, workers = parallelly::availableCores() - 1)
  
  final_result <- progressr::with_progress({
    p <- progressr::progressor(steps = total_steps)
    
    lapply(gamma, function(g){
      varD <- g*varA
      
      if(!is.null(kinshipD)){
        gD_true <- mvrnorm(n = nsim, mu = rep(0, Nc), Sigma = varD*kinshipD)
      } else {
        gD_true <- matrix(0, nrow = nsim, ncol = Nc)
      }
      
      lapply(h, function(h2){
        varE <- ((varA+varD)/h2)-(varA+varD)
        
        res_nsim <- future.apply::future_lapply(seq_len(nsim), function(s){
          gs <- gA_true[s, ]+gD_true[s, ]
          ps <- mu+gs+rnorm(Nc, mean = 0, sd = sqrt(varE))
          
          size_res <- lapply(TRSSize, function(sz){
            
            # random
            random_idx <- sample(Nc, sz)
            Random_TRS <- data.frame(posIndex = random_idx)
            
            dplyr::bind_rows(
              Evaluation(kinshipA=kinshipA, kinshipD=kinshipD, OptTRS=OptTRS$EI, 
                         TRSSize=sz, p_true=ps, g_true=gs, mu_true=mu, k=k) |> 
                mutate(strategy="EI"),
              Evaluation(kinshipA=kinshipA, kinshipD=kinshipD, OptTRS=OptTRS$UCB, 
                         TRSSize=sz, p_true=ps, g_true=gs, mu_true=mu, k=k) |> 
                mutate(strategy="UCB"),
              Evaluation(kinshipA=kinshipA, kinshipD=kinshipD, OptTRS=OptTRS$GV, 
                         TRSSize=sz, p_true=ps, g_true=gs, mu_true=mu, k=k) |> 
                mutate(strategy="GVaverage"),
              Evaluation(kinshipA=kinshipA, kinshipD=kinshipD, OptTRS=Random_TRS, 
                         TRSSize=sz, p_true=ps, g_true=gs, mu_true=mu, k=k) |> 
                mutate(strategy="Random"),
            ) |> mutate(TRSSize = sz)
          })
          p()
          return(dplyr::bind_rows(size_res))
        }, future.packages = c("dplyr"))
        
        return(dplyr::bind_rows(res_nsim) |> mutate(heritability=h2, gamma=g))
      }) |> dplyr::bind_rows()
    }) |> dplyr::bind_rows()
  })
  future::plan(future::sequential)
  
  Result <- final_result |> 
    dplyr::group_by(gamma, heritability, TRSSize, strategy, metric) |> 
    dplyr::summarize(
      Mean = mean(value, na.rm = TRUE),
      SD = sd(value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(k = k)
  
  return(list(Result = Result))
}

#--- Simulation3 (discussion) ---#
Simulation3 <- function(kinshipA, kinshipD = NULL, yReal, OptTRS, TRSSize, k, mu_true, g_true, e_true, p_true, 
                        dataSet = "Unknown_DataSet", trait = "Unknown_Trait", output_dir = "Simulation_Results"){
  
  suppressPackageStartupMessages({
    library(BGLR)
    library(dplyr)
    library(tibble)
    library(purrr)
    library(MASS)
    requireNamespace("progressr", quietly = TRUE)
    requireNamespace("cli", quietly = TRUE)
    requireNamespace("future", quietly = TRUE)
    requireNamespace("future.apply", quietly = TRUE)
  })
  
  #--- progress bar setting ---#
  pb <- progressr::handler_progress(
    format = " Simulation Process [:bar] :percent | ETA: :eta", 
    complete = "=", 
    incomplete = " ", 
    current = ">", 
    width = 100, 
    clear = FALSE
  )
  
  options(future.globals.maxSize = +Inf)
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  Nc <- nrow(kinshipA)
  
  #--- Parallel Plan Setting ---#
  future::plan(future::multisession)
  
  #--- Simulation Loop ---#
  final <- list()
  
  progressr::handlers(pb)
  
  cli::cli_h1("Running Evaluation Loop")
  
  for(size in TRSSize){
    
    cli::cli_h2(sprintf("Processing Training Set Size: %d", size))
    
    res_size <- progressr::with_progress({
      p <- progressr::progressor(steps = 2000)
      
      res_list <- future.apply::future_lapply(seq_len(2000), function(s){
        
        tryCatch({
          ps <- p_true[s, ]; gs <- g_true[s, ]
          
          ## Evaluations
          res_eva <- dplyr::bind_rows(
            # EI
            Evaluation(kinshipA = kinshipA, kinshipD = kinshipD, OptTRS = OptTRS$EI, TRSSize = size, 
                       p_true = ps, g_true = gs, mu_true = mu_true, k = k) |> 
              dplyr::mutate(strategy = "EI"), 
            
            # UCB
            Evaluation(kinshipA = kinshipA, kinshipD = kinshipD, OptTRS = OptTRS$UCB, TRSSize = size, 
                       p_true = ps, g_true = gs, mu_true = mu_true, k = k) |> 
              dplyr::mutate(strategy = "UCB")
          )
          
          res_eva$TRSSize <- size
          p()
          return(res_eva)
          
        }, error = function(e){
          message(sprintf("\n[ERROR] Sim %d failed: %s", s, e$message))
          p()
          return(tibble::tibble(
            metric = c("NDCG", "RSR", "SRC"), 
            value = NA_real_, 
            strategy = "Error", 
            TRSSize = size
          ))
        })
        
      }, 
      future.seed = TRUE, 
      future.packages = c("BGLR", "MASS", "dplyr", "tibble"), 
      future.globals = TRUE)
      
      dplyr::bind_rows(res_list)
    })
    
    res_size <- res_size |> 
      dplyr::mutate(dataset = dataSet, trait = trait, heritability = NA_real_, k = k, .before = 1)
    
    saveName <- file.path(output_dir, sprintf("Res_%s_%s_Size%d.RData", dataSet, trait, size))
    save(res_size, file = saveName)
    cli::cli_alert_success(sprintf("Saved results for size %d to %s", size, basename(saveName)))
    
    final[[as.character(size)]] <- res_size
  }
  future::plan(future::sequential)
  
  #--- Combine all heritability results ---#
  cli::cli_h1("Summarizing...")
  
  all_table <- dplyr::bind_rows(final)
  
  Optimal <- all_table |>  
    dplyr::group_by(dataset, trait, strategy, heritability, TRSSize, metric, k) |> 
    dplyr::summarize(
      Mean = mean(value, na.rm = TRUE), 
      SD = sd(value, na.rm = TRUE), 
      .groups = "drop"
    )
  
  cli::cli_alert_success(sprintf("Simulation Done! All results saved in '%s'", output_dir))
  cli::cli_rule()
  
  #--- Return ---#
  return(list(Optimal = Optimal, true_mean = mu_true))
}
