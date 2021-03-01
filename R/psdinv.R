psdinv <- 
  function(x){
    # Positive Semi-Definite Matrix Inverse
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: February 15, 2021
    
    x <- as.matrix(x)
    tol <- .Machine$double.eps
    xeig <- eigen(x, symmetric = TRUE)
    xrnk <- sum(xeig$values > xeig$values[1] * tol * nrow(x))
    xinv <- tcrossprod(xeig$vectors[, 1:xrnk, drop = FALSE] %*% 
                         diag(1/sqrt(xeig$values[1:xrnk]), nrow = xrnk, ncol = xrnk))
    xinv
    
  } # end psdinv