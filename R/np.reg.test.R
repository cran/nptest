np.reg.test <-
  function(x, y, z = NULL, method = NULL,
           beta = NULL, homosced = FALSE,
           R = 9999, parallel = FALSE, cl = NULL,
           perm.dist = TRUE){
    # Nonparametric Tests of Regression Coefficients
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: November 4, 2018
    
    
    ### check z
    znull <- is.null(z)
    
    ### check method
    if(is.null(method)) method <- ifelse(znull, "perm", "HJ")
    
    ### permutation test
    if(znull){
      pt <- rand.test.lm1(x = x, y = y,
                          method = method,
                          beta = beta, homosced = homosced,
                          R = R, parallel = parallel, cl = cl,
                          perm.dist = perm.dist)
    } else {
      pt <- rand.test.lm2(x = x, y = y, z = z,
                          method = method,
                          beta = beta, homosced = homosced,
                          R = R, parallel = parallel, cl = cl,
                          perm.dist = perm.dist)
    }
    
    ### return results
    class(pt) <- "np.reg.test"
    return(pt)
    
  } # end np.reg.test.R