rand.test.one <-
  function(x, 
           alternative = c("two.sided", "less", "greater"),
           mu = 0, symmetric = TRUE, median.test = FALSE,
           R = 9999, parallel = FALSE, cl = NULL,
           perm.dist = TRUE){
    # One-Sample Randomization Tests for Location (Mean/Median)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: September 9, 2020
    
    
    #########   INITIAL CHECKS   #########
    
    ### check x
    x <- as.matrix(x)
    n <- nrow(x)
    nvar <- ncol(x)
    
    ### check alternative
    alternative <- as.character(alternative[1])
    alternative <- pmatch(alternative, c("two.sided", "less", "greater"))
    if(is.na(alternative)) stop("Invalid 'alternative' input.")
    alternative <- c("two.sided", "less", "greater")[alternative]
    
    ### check mu
    mu <- as.numeric(mu)
    if(length(mu) != nvar) mu <- rep(mu, length.out = nvar)
    
    ### check symmetric
    symmetric <- as.logical(symmetric[1])
    
    ### check median.test
    median.test <- as.logical(median.test[1])
    
    ### check R
    R <- as.integer(R)
    if(R < 1) stop("Input 'R' must be a positive integer.")
    
    ### check parallel
    parallel <- as.logical(parallel[1])
    
    ### check 'cl'
    make.cl <- FALSE
    if(parallel){
      if(is.null(cl)){
        make.cl <- TRUE
        cl <- makeCluster(detectCores())
      } else {
        if(!any(class(cl) == "cluster")) stop("Input 'cl' must be an object of class 'cluster'.")
      }
    }
    
    ### exact or approximate?
    suppressWarnings( nperm <- 2^n )
    exact <- ifelse(nperm <= R + 1L, TRUE, FALSE)
    if(exact){
      ix <- flipn(n)
      nperm <- ncol(ix)
    }
    
    #########   LOCATION TEST   #########
    
    ### univariate or multivariate
    if(nvar == 1L){
      
      ## UNIVARIATE TEST
      
      x <- as.numeric(x)
      
      ## non-zero null hypothesis?
      if(abs(mu) > 0) {
        x <- x - mu
      }
      
      ## observed test statistic
      Tstat <- Tstat.one(x = x, symmetric = symmetric,
                         median.test = median.test)
      
      ## permutation distribution
      if(exact){
        
        # parallel or sequential computation?
        if(parallel){
          permdist <- parCapply(cl = cl, x = ix, 
                                FUN = Tperm.one, xvec = x, 
                                symmetric = symmetric,
                                median.test = median.test,
                                exact = exact)
        } else {
          permdist <- apply(X = ix, MARGIN = 2, 
                            FUN = Tperm.one, xvec = x, 
                            symmetric = symmetric,
                            median.test = median.test,
                            exact = exact)
        } # end if(parallel)
        
      } else {
        
        # approximate permutation test (given input R)
        nperm <- R + 1L
        permdist <- rep(0, nperm)
        permdist[1] <- Tstat
        
        # parallel or sequential computation?
        if(parallel){
          permdist[2:nperm] <- parSapply(cl = cl, X = integer(R), 
                                         FUN = Tperm.one, xvec = x, 
                                         symmetric = symmetric,
                                         median.test = median.test,
                                         exact = exact)
        } else {
          permdist[2:nperm] <- sapply(X = integer(R),
                                      FUN = Tperm.one, xvec = x, 
                                      symmetric = symmetric,
                                      median.test = median.test,
                                      exact = exact)
        } # end if(parallel)
        
      } # end if(exact)
      
      ## permutation p-value
      if(alternative == "less"){
        p.value <- mean(permdist <= Tstat)
      } else if(alternative == "greater"){
        p.value <- mean(permdist >= Tstat)
      } else {
        p.value <- mean(abs(permdist) >= abs(Tstat))
      }
      
      ## estimate
      if(median.test){
        if(symmetric){
          ltmat <- lower.tri(diag(n), diag = TRUE)
          walsh.avg <- outer(x, x, FUN = "+")[ltmat] / 2
          est <- median(walsh.avg) + mu
        } else {
          est <- median(x) + mu
        }
      } else {
        est <- mean(x) + mu
      }
      
      
    } else {
      
      ## MULTIVARIATE TEST
      
      ## non-zero null hypothesis?
      if(max(abs(mu)) > 0) {
        for(v in 1:nvar) x[,v] <- x[,v] - mu[v]
      }
      
      ## observed test statistic
      Tuni <- Tstat.one.mv(x = x, symmetric = symmetric,
                           median.test = median.test, combine = FALSE)
      Tstat <- ifelse(alternative == "two.sided", max(abs(Tuni)),
                      ifelse(alternative == "greater", max(Tuni), min(Tuni)))
      
      ## permutation distribution
      if(exact){
        
        # parallel or sequential computation?
        if(parallel){
          permdist <- parCapply(cl = cl, x = ix, 
                                FUN = Tperm.one.mv, xmat = x, 
                                symmetric = symmetric,
                                median.test = median.test,
                                alternative = alternative,
                                exact = exact)
        } else {
          permdist <- apply(X = ix, MARGIN = 2, 
                            FUN = Tperm.one.mv, xmat = x, 
                            symmetric = symmetric,
                            median.test = median.test,
                            alternative = alternative,
                            exact = exact)
        } # end if(parallel)
        
      } else {
        
        # approximate permutation test (given input R)
        nperm <- R + 1L
        permdist <- rep(0, nperm)
        permdist[1] <- Tstat
        
        # parallel or sequential computation?
        if(parallel){
          permdist[2:nperm] <- parSapply(cl = cl, X = integer(R), 
                                         FUN = Tperm.one.mv, xmat = x, 
                                         symmetric = symmetric,
                                         median.test = median.test,
                                         alternative = alternative,
                                         exact = exact)
        } else {
          permdist[2:nperm] <- sapply(X = integer(R),
                                      FUN = Tperm.one.mv, xmat = x, 
                                      symmetric = symmetric,
                                      median.test = median.test,
                                      alternative = alternative,
                                      exact = exact)
        } # end if(parallel)
        
      } # end if(exact)
      
      ## permutation p-value
      uni.p.value <- rep(NA, nvar)
      if(alternative == "less"){
        p.value <- mean(permdist <= Tstat)
        for(v in 1:nvar) uni.p.value[v] <- mean(permdist <= Tuni[v])
      } else if(alternative == "greater"){
        p.value <- mean(permdist >= Tstat)
        for(v in 1:nvar) uni.p.value[v] <- mean(permdist >= Tuni[v])
      } else {
        absperm <- abs(permdist)
        p.value <- mean(absperm >= abs(Tstat))
        for(v in 1:nvar) uni.p.value[v] <- mean(absperm >= abs(Tuni[v]))
      }
      
      ## estimate
      est <- rep(NA, nvar)
      if(median.test){
        if(symmetric){
          ltmat <- lower.tri(diag(n), diag = TRUE)
          for(v in 1:nvar){
            walsh.avg <- outer(x[,v], x[,v], FUN = "+")[ltmat] / 2
            est[v] <- median(walsh.avg) + mu[v]
          }
        } else {
          for(v in 1:nvar) est[v] <- median(x[,v]) + mu[v]
        }
      } else {
        for(v in 1:nvar) est[v] <- mean(x[,v]) + mu[v]
      }
      
    } # end if(nvar == 1L)
    
    ### return results
    if(make.cl) stopCluster(cl)
    if(!perm.dist) permdist <- NULL
    res <- list(statistic = Tstat, p.value = p.value,
                perm.dist = permdist, alternative = alternative, 
                null.value = mu, symmetric = symmetric, 
                median.test = median.test, 
                R = nperm - ifelse(exact, 0, 1), exact = exact,
                estimate = est)
    if(nvar > 1L) {
      res$univariate <- Tuni
      res$adj.p.values <- uni.p.value
    }
    class(res) <- "rand.test.one"
    return(res)
    
  } # end rand.test.one.R


### permutation replication (univariate)
Tperm.one <-
  function(i, xvec, symmetric = FALSE, 
           median.test = FALSE, exact = FALSE){
    if(!exact) i <- sample(c(-1, 1), size = length(xvec), replace = TRUE)
    Tstat <- Tstat.one(x = xvec * i, 
                       symmetric = symmetric, 
                       median.test = median.test)
    return(Tstat)
  } # end Tperm.one.R


### permutation replication (multivariate)
Tperm.one.mv <-
  function(i, xmat, symmetric = FALSE, median.test = FALSE, 
           alternative = "two.sided", combine = TRUE, exact = FALSE){
    if(!exact) i <- sample(c(-1, 1), size = nrow(xmat), replace = TRUE)
    Tstat <- Tstat.one.mv(x = xmat * i, 
                          symmetric = symmetric, 
                          median.test = median.test,
                          alternative = alternative,
                          combine = combine)
    return(Tstat)
  } # end Tperm.one.mv.R


### test statistic (univariate)
Tstat.one <- 
  function(x, symmetric = FALSE, median.test = FALSE){
    n <- length(x)
    if(median.test){
      if(symmetric){
        rx <- rank(abs(x))
        Tplus <- sum(rx[x > 0])
        E.Tplus <- n * (n + 1) / 4
        V.Tplus <- n * (n + 1) * (2*n + 1) / 24
        Tstat <- (Tplus - E.Tplus) / sqrt(V.Tplus)
      } else {
        Tplus <- sum(x > 0)
        E.Tplus <- n / 2
        V.Tplus <- n / 4
        Tstat <- (Tplus - E.Tplus) / sqrt(V.Tplus)
      }
    } else {
      xbar <- mean(x)
      xvar <- sum((x - xbar)^2) / (n - 1)
      Tstat <- xbar / sqrt(xvar / n)
      if(!symmetric){
        skew <- sum((x - xbar)^3) / (n - 1)
        Tshift <- ((xbar / xvar)^2 / 3) + (1 / (6 * xvar * n))
        Tstat <- Tstat + skew * Tshift / sqrt(xvar / n)
      } 
    } # end if(median.test)
    Tstat
  } # end Tstat.one.R


### test statistic  (multivariate)
Tstat.one.mv <- 
  function(x, symmetric = FALSE, median.test = FALSE,
           alternative = "two.sided", combine = TRUE){
    n <- nrow(x)
    nvar <- ncol(x)
    if(median.test){
      if(symmetric){
        Tstat <- rep(0, nvar)
        E.Tplus <- n * (n + 1) / 4
        S.Tplus <- sqrt(n * (n + 1) * (2*n + 1) / 24)
        for(v in 1:nvar){
          rx <- rank(abs(x[,v]))
          Tstat[v] <- (sum(rx * (x[,v] > 0)) - E.Tplus) / S.Tplus
        }
      } else {
        Tplus <- colSums(x > 0)
        E.Tplus <- n / 2
        V.Tplus <- n / 4
        Tstat <- (Tplus - E.Tplus) / sqrt(V.Tplus)
      }
    } else {
      xbars <- colMeans(x)
      Tstat <- rep(0, nvar)
      if(symmetric){
        xvars <- colSums((x - matrix(xbars, nrow = n, ncol = nvar, byrow = TRUE))^2) / (n - 1)
        Tstat <- xbars / sqrt(xvars / n)
      } else {
        x <- x - matrix(xbars, nrow = n, ncol = nvar, byrow = TRUE)
        xsq <- x * x
        xvars <- colSums(xsq) / (n - 1)
        skews <- colSums(xsq * x) / (n - 1)
        shift <- ((xbars / xvars)^2 / 3) + (1 / (6 * xvars * n))
        Tstat <- (xbars + skews * shift) / sqrt(xvars / n)
      }
    } # end if(median.test)
    if(combine) {
      Tstat <- ifelse(alternative == "two.sided", max(abs(Tstat)),
                      ifelse(alternative == "greater", max(Tstat), min(Tstat)))
    }
    as.numeric(Tstat)
  } # end Tstat.one.mv.R

