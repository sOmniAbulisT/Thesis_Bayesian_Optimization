#' Calculate A-optimality-like Criterion
#' 
#' @description
#' This function calculates a ranking criterion based on the diagonal elements of the 
#' combined genetic covariance matrix (Additive + optional Dominance).
#' By selecting individuals with higher diagonal values (higher genetic variance), 
#' this criterion approximates an A-optimality approach to maximize the information context
#' of the training set. 
#' 
#' @param kinshipA Matrix. Additive kinship matrix.
#' @param kinshipD Matrix (Option). Dominance kinship matrix.
#' @param sigmaA   Numeric. The additive variance component (weight for additive matrix, default = 1).
#' @param sigmaD   Numeric. The dominance variance component (weight for dominance matrix, default = 1).
#' 
#' @return A data frame sorted by the GV score in decreasing order, containing:
#' \item{speciesName}{Character. The names of the individuals.}
#' \item{averageGV}{Numeric. The calculated criterion score (diagonal element). }
#' \item{posIndex}{Integer. The original index of the individual in the input matrix.}
#' 
#' @export

GVave <- function(kinshipA, kinshipD = NULL, sigmaA = 1, sigmaD = 1){
  
  if(!is.null(kinshipD)){
    K <- sigmaA * kinshipA + sigmaD * kinshipD
  } else {
    K <- sigmaA * kinshipA
  }
  ave <- diag(K)
  
  ord_GV <- order(ave, decreasing = TRUE)
  res_GV <- data.frame(
    speciesName = rownames(kinshipA)[ord_GV], 
    averageGV = ave[ord_GV], 
    posIndex = ord_GV
  )
  
  return(res_GV)
}
