CalGEBV <- function(kinshipA, kinshipD = NULL, OptTRS, TRSSize, p_true){
  Nc <- nrow(kinshipA)
  TRS_index <- OptTRS$posIndex[seq_len(TRSSize)]
  
  kAtt <- kinshipA[TRS_index, TRS_index] |> as.matrix()
  y_train <- p_true[TRS_index]
  
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
  fit_BGLR <- BGLR(y = y_train, 
                   ETA = ETA_list, 
                   verbose = FALSE, nIter = 5000, burnIn = 1000, 
                   saveAt = tmp_prefix)
  
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
  
  ord_gebv <- order(gebv, decreasing = TRUE)
  res_gebv <- data.frame(
    speciesName = rownames(kinshipA)[ord_gebv], 
    GEBV = gebv[ord_gebv], 
    posIndex = ord_gebv, 
    stringsAsFactors = FALSE
  )
  
  return(res_gebv)
}

Validation <- function(kinshipA, kinshipD = NULL, yReal, OptTRS, nsim = 2000, TRSSize, k, 
                       dataSet = "Unknown_DataSet", trait = "Unknown_Trait", 
                       output_dir = "Validation_Results"){
  
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
  OptTRS <- data.frame(posIndex = OptTRS)
  
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
          
          ## Evaluations
          res_eva <- Evaluation(kinshipA = kinshipA, kinshipD = kinshipD, OptTRS = OptTRS, TRSSize = size, 
                                p_true = ps, g_true = gs, mu_true = mu_true, k = k) |> 
                      dplyr::mutate(strategy = "EI")
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
  
  cli::cli_alert_success(sprintf("Validation Done! All results saved in '%s'", output_dir))
  cli::cli_rule()
  
  #--- Return ---#
  return(list(Optimal = Optimal, true_mean = mu_true))
}
