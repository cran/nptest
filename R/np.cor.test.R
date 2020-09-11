np.cor.test <- 
  function(x, y, z = NULL,
           alternative = c("two.sided", "less", "greater"),
           rho = 0, independent = FALSE, partial = TRUE,
           R = 9999, parallel = FALSE, cl = NULL,
           perm.dist = TRUE){
    # Nonparametric Tests of Correlation Coefficients
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: July 8, 2020
    
    
    ### check x and y
    x <- as.numeric(x)
    y <- as.numeric(y)
    n <- length(x)
    if(length(y) != n) stop("Inputs 'x' and 'y' must have same length.")
    
    ### check z
    if(!is.null(z)){
      z <- cbind(1, z)
      if(nrow(z) != n) stop("Input 'z' is incompatible with inputs 'x' and 'y'.")
      zsvd <- svd(z, nv = 0)
      x <- x - zsvd$u %*% crossprod(zsvd$u, x)
      if(partial) y <- y - zsvd$u %*% crossprod(zsvd$u, y)
    }
    
    ### permutation test
    pt <- rand.test.cor(x = x, y = y, 
                        alternative = alternative,
                        rho = rho, independent = independent,
                        R = R, parallel = parallel, cl = cl,
                        perm.dist = perm.dist)
    
    ### correct se?
    if(independent && !is.null(z)){
      correct <- sqrt((n - ncol(z) - 1) / (n - 2))
      pt$statistic <- pt$statistic * correct
      if(perm.dist) pt$perm.dist <- pt$perm.dist * correct
    }
    
    ### return results
    class(pt) <- "np.cor.test"
    return(pt)
    
  } # end np.cor.test.R
