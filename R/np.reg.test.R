np.reg.test <-
  function(x, y, z = NULL, method = NULL,
           beta = NULL, homosced = FALSE, lambda = 0,
           R = 9999, parallel = FALSE, cl = NULL,
           perm.dist = TRUE){
    # Nonparametric Tests of Regression Coefficients
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: February 17, 2021
    
    
    ### check z
    znull <- is.null(z)
    
    ### check method
    if(is.null(method)) method <- ifelse(znull, "perm", "HJ")
    
    ### permutation test
    if(znull){
      pt <- rand.test.lm1(x = x, y = y,
                          method = method,
                          beta = beta, homosced = homosced,
                          lambda = lambda, R = R, 
                          parallel = parallel, cl = cl,
                          perm.dist = perm.dist)
    } else {
      pt <- rand.test.lm2(x = x, y = y, z = z,
                          method = method,
                          beta = beta, homosced = homosced,
                          lambda = lambda, R = R, 
                          parallel = parallel, cl = cl,
                          perm.dist = perm.dist)
      if(homosced && method %in% c("KC", "SW") && R > 0){
        x <- as.matrix(x)
        z <- as.matrix(z)
        n <- nrow(x)
        p <- ncol(x)
        correct <- (n - p - ncol(z) - 1) / (n - p - 1)
        pt$statistic <- pt$statistic * correct
        if(perm.dist) pt$perm.dist <- pt$perm.dist * correct
      }
    }
    
    ### return results
    class(pt) <- "np.reg.test"
    return(pt)
    
  } # end np.reg.test.R