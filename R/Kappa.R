Kappa <-
function(beta, X, Y, penalty,contemp_gamma_nonconvex, lambda_kappa,regularize_mat_kappa){
  if (missing(regularize_mat_kappa)){
    regularize_mat_kappa <- matrix(TRUE, ncol(Y), ncol(Y))
    diag(regularize_mat_kappa) <- FALSE
  }
  n <- nrow(Y)  
  SigmaR <- 1/n * t(Y - X %*% beta) %*% (Y - X %*% beta)
  if (any(eigen(SigmaR,only.values = TRUE)$values < -sqrt(.Machine$double.eps))){
    stop("Residual covariance matrix is not non-negative definite")
  } 
  if (penalty == "lasso") {
    
    lambda_mat <- regularize_mat_kappa * lambda_kappa
  } else {
    
    kappa_prep <- abs(solve(SigmaR))
    
    if (penalty == "atan") {
      lambda_mat <- lambda_kappa * (contemp_gamma_nonconvex * (contemp_gamma_nonconvex + 2/pi)) /
        (contemp_gamma_nonconvex^2 + kappa_prep^2)
      
    } else if (penalty == "scad") {
      t <- pmax(kappa_prep, 1e-12)
      lambda_mat <- ifelse(
        t <= lambda_kappa,
        lambda_kappa,
        ifelse(t <= contemp_gamma_nonconvex * lambda_kappa,
               (contemp_gamma_nonconvex * lambda_kappa - t) /
                 ((contemp_gamma_nonconvex - 1)),
               0)
      )
    }
    
    # apply regularization matrix
    lambda_mat <- lambda_mat * regularize_mat_kappa
  }
  
  res <- glasso(SigmaR, lambda_mat)
  return(as.matrix(forceSymmetric(res$wi)))
}

