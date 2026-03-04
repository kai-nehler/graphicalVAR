Rothmana <-
function(X, Y, lambda_beta, lambda_kappa, penalty, 
         regression_gamma_nonconvex = NULL, contemp_gamma_nonconvex = NULL,
         regularize_mat_beta, regularize_mat_kappa, convergence = 1e-4, gamma = 0.5, maxit.in = 100, maxit.out = 100,
         penalize.diagonal, # if FALSE, penalizes the first diagonal (assumed to be auto regressions), even when ncol(X) != ncol(Y) !
         interceptColumn = 1, # Set to NULL or NA to omit
         mimic = "current",
         start_beta = c("empty", "ridge"),
         ridge_correction = c("none", "sample_size"),
         likelihood = c("unpenalized","penalized")
         ){
  # Algorithm 2 of Rothmana, Levinaa & Ji Zhua
  
  likelihood <- match.arg(likelihood)
  
  # Defaults for nonconvex shape parameter (used for BOTH atan and scad)
  if (penalty %in% c("atan","scad")) {
    if (is.null(regression_gamma_nonconvex)) {
      regression_gamma_nonconvex <- if (penalty == "atan") 0.5 else 3.7
    }
    if (is.null(contemp_gamma_nonconvex)) {
      contemp_gamma_nonconvex <- if (penalty == "atan") 0.5 else 3.7
    }
  }
  
  nY <- ncol(Y)
  nX <- ncol(X)
 
  if (missing(penalize.diagonal)){
    if (mimic == "0.1.2"){
      penalize.diagonal <- nY != nX
    } else {
      penalize.diagonal <- (nY != nX-1) & (nY != nX ) 
    }
  }
  
  # Add regularization matrix for regression coefficients:
  # Initial lambda_mat is redundant for nonconvex penalties; 
  # it is recomputed from beta and masked by regularize_mat_beta at the start of 
  # each iteration.
  # Without regularization matrix - everything is subject to regularization
  if (missing(regularize_mat_beta)){
    lambda_mat <- matrix(lambda_beta,nX, nY)
    if (!penalize.diagonal){
      if (nY == nX){
        add <- 0
      } else if (nY == nX - 1){
        add <- 1
      } else {
        stop("Beta is not P x P or P x P+1, cannot detect diagonal.")
      }
      for (i in 1:min(c(nY,nX))){
        lambda_mat[i+add,i] <- 0
      }
    }
    
    # by default, regularize everything (later: diagonal & intercept handled below)
    regularize_mat_beta <- matrix(TRUE, nY, nX) 
    
  } else {
    lambda_mat <- lambda_beta * t(regularize_mat_beta)
    if (nrow(lambda_mat) == nX-1){
      lambda_mat <- rbind(FALSE,lambda_mat)
    }
    if (nrow(lambda_mat) != nX){
      browser()
      stop("Number of rows in 'regularize_mat_beta' is incorrect.")
    }
    
    if (ncol(lambda_mat) != nY){
      stop("Number of columns in 'regularize_mat_beta' is incorrect.")
    }
  }
  
  
  if (!is.null(interceptColumn) && !is.na(interceptColumn)){
    lambda_mat[interceptColumn,] <- 0
  }
 
  n <- nrow(X)
  # ridge correction for starting values:
  if (ridge_correction == "sample_size") {
    beta_ridge <- beta_ridge_C(X, Y, lambda_beta * n)
  } else if (ridge_correction == "none") {
    beta_ridge <- beta_ridge_C(X, Y, lambda_beta)
  }
  
  # Starting values:
  if (start_beta == "ridge") {
    beta <- beta_ridge
  } else if (start_beta == "empty") {
    beta <- matrix(0, nX, nY)  
  }
  
  # Algorithm:
  it <- 0

  repeat{
    it <- it + 1
    kappa <- Kappa(beta, X, Y, penalty, contemp_gamma_nonconvex, lambda_kappa, regularize_mat_kappa)
    
    # Update lambda_mat only for nonconvex penalties
    if (penalty %in% c("atan","scad")) {
      
      if (penalty == "atan") {
        lambda_mat <- (regression_gamma_nonconvex * (regression_gamma_nonconvex + 2/pi)) /
          (regression_gamma_nonconvex^2 + beta^2)
        
      } else if (penalty == "scad") {
        t <- pmax(abs(beta), 1e-12)
        lambda_mat <- ifelse(
          t <= lambda_beta,
          lambda_beta,
          ifelse(t <= regression_gamma_nonconvex * lambda_beta,
                 (regression_gamma_nonconvex * lambda_beta - t) /
                   ((regression_gamma_nonconvex - 1)),
                 0)
        )
      }
      
      # apply either default or user mask
      lambda_mat <- lambda_mat * t(regularize_mat_beta)
      
      
      # allow mask without intercept row
      if (nrow(lambda_mat) == nX - 1) {
        lambda_mat <- rbind(FALSE, lambda_mat)
      }
      if (nrow(lambda_mat) != nX || ncol(lambda_mat) != nY) {
        stop("lambda_mat has wrong dimensions after applying regularize_mat_beta.")
      }
      
      # diagonal handling
      if (!penalize.diagonal) {
        if (!exists("add", inherits = FALSE)) {
          # only needed if missing(regularize_mat_beta) AND !penalize.diagonal
          if (nY == nX) {
            add <- 0
          } else if (nY == nX - 1) {
            add <- 1
          } else {
            stop("Beta is not P x P or P x P+1, cannot detect diagonal.")
          }
        }
        for (i in 1:min(nY, nX)) lambda_mat[i + add, i] <- 0
      }
      
      # intercept row never penalized
      if (!is.null(interceptColumn) && !is.na(interceptColumn)) {
        lambda_mat[interceptColumn, ] <- 0
      }
    }
    
    
    beta_old <- beta
    beta <- Beta_C(kappa, beta, X, Y, lambda_beta, lambda_mat, convergence, maxit.in) 
    
    if (sum(abs(beta - beta_old)) < (convergence * sum(abs(beta_ridge)))){
      break
    }
    
    if (it > maxit.out){
      warning("Model did NOT converge in outer loop")
      break
    }
  }
  
  ## Compute unconstrained kappa (codes from SparseTSCGM):
  ZeroIndex <- which(kappa==0, arr.ind=TRUE) ## Select the path of zeros
  WS <-  (t(Y)%*%Y - t(Y) %*% X  %*% beta - t(beta) %*% t(X)%*%Y + t(beta) %*% t(X)%*%X %*% beta)/(nrow(X))
  
  if (any(eigen(WS,only.values = TRUE)$values < -sqrt(.Machine$double.eps))){
    stop("Residual covariance matrix is not non-negative definite")
  }
  
  if (likelihood == "unpenalized"){
    if (nrow(ZeroIndex)==0){
      out4 <- suppressWarnings(glasso(WS, rho = 0, trace = FALSE))
    } else {
      out4 <- suppressWarnings(glasso(WS, rho = 0, zero = ZeroIndex,
                                      trace = FALSE))
    }
    lik1  <- determinant( out4$wi)$modulus[1]
    lik2 <- sum(diag( out4$wi%*%WS))
  } else {
    lik1  <- determinant( kappa )$modulus[1]
    lik2 <- sum(diag( kappa%*%WS))
  }

  pdO = sum(sum(kappa[upper.tri(kappa,diag=FALSE)] !=0))
  if (mimic == "0.1.2"){
    pdB = sum(sum(beta !=0))
  } else {
    pdB = sum(sum(beta[lambda_mat!=0] !=0)) 
  }
  
  LLk <-  (n/2)*(lik1-lik2) 
  LLk0 <-  (n/2)*(-lik2)
  
  EBIC <-  -2*LLk + (log(n))*(pdO +pdB) + (pdO  + pdB)*4*gamma*log(2*nY)

  
  ### TRANSPOSE BETA!!!
  return(list(beta=t(beta), kappa=kappa, EBIC = EBIC, it = it))
}