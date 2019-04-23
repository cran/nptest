rand.test.lm1 <-
  function(x, y, 
           method = c("perm", "flip", "both"),
           beta = NULL, homosced = FALSE,
           R = 9999, parallel = FALSE, cl = NULL,
           perm.dist = TRUE){
    # Randomization Test for Regression (w/o covariates)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: November 4, 2018
    
    
    #########   INITIAL CHECKS   #########
    
    ### check x and y
    x <- as.matrix(x)
    y <- as.matrix(y)
    nvar <- ncol(y)
    n <- nrow(x)
    p <- ncol(x)
    if(nrow(y) != n) stop("Inputs 'x' and 'y' are incompatible.")
    
    ### check method
    method <- as.character(method[1])
    method.options <- c("perm", "flip", "both")
    method <- method.options[pmatch(method, method.options)]
    if(is.na(method)) stop("Invalid 'method' input.")
    
    ### check beta
    if(!is.null(beta)){
      beta <- as.matrix(beta)
      if(nrow(beta) != p) stop("Inputs 'beta' and 'x' are incompatible.")
      if(ncol(beta) != nvar) stop("Inputs 'beta' and 'y' are incompatible.")
    }
    
    ### check homosced
    homosced <- as.logical(homosced[1])
    
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
    if(method == "perm"){
      suppressWarnings( nperm <- factorial(n) )
    } else if(method == "sign"){
      suppressWarnings( nperm <- 2^n )
    } else {
      suppressWarnings( np <- factorial(n) )
      suppressWarnings( ns <- 2^n )
      nperm <- np * ns
    }
    exact <- ifelse(nperm <= R + 1L, TRUE, FALSE)
    if(exact){
      if(method == "perm"){
        ix <- permn(n)
        nperm <- ncol(ix)
      } else if(method == "flip"){
        ix <- flipn(n)
        nperm <- ncol(ix)
      } else {
        ixp <- permn(n)
        ixs <- flipn(n)
        igrid <- expand.grid(perm = 1:ncol(ixp), sign = 1:ncol(ixs))
        ix <- apply(igrid, 1, function(u) c(ixp[,u[1]], ixs[,u[2]]))
      }
    }
    
    
    #########   REGRESSION TEST   #########
    
    ### non-zero null hypothesis?
    if(!is.null(beta)) y <- y - x %*% beta
    
    ### center data
    xbar <- colMeans(x)
    ybar <- colMeans(y)
    x <- scale(x, center = xbar, scale = FALSE)
    y <- scale(y, center = ybar, scale = FALSE)
    
    ### calculate crossproducts and inverses
    xtx <- crossprod(x)
    xtxi <- solve(xtx)
    xinv <- tcrossprod(xtxi, x)
    coefs <- xinv %*% y
    ssy <- NULL
    if(homosced) ssy <- colSums(y^2)
    
    
    ### univariate or multivariate
    if(nvar == 1L){
      
      ## UNIVARIATE TEST
      
      y <- as.numeric(y)
      
      ## observed test statistic
      Tstat <- Tstat.lm1(x = x, y = y, homosced = homosced,
                         xtx = xtx, xtxi = xtxi, xinv = xinv, 
                         ssy = ssy, beta = coefs)
      
      ## permutation distribution
      if(exact){
        
        # parallel or sequential computation?
        if(parallel){
          permdist <- parCapply(cl = cl, x = ix, 
                                FUN = Tperm.lm1, 
                                xmat = x, yvec = y, 
                                method = method, 
                                homosced = homosced,
                                exact = exact, xtx = xtx, 
                                xtxi = xtxi, xinv = xinv, 
                                ssy = ssy)
        } else {
          permdist <- apply(X = ix, MARGIN = 2, 
                            FUN = Tperm.lm1, 
                            xmat = x, yvec = y, 
                            method = method, 
                            homosced = homosced,
                            exact = exact, xtx = xtx, 
                            xtxi = xtxi, xinv = xinv, 
                            ssy = ssy)
        } # end if(parallel)
        
      } else {
        
        # approximate permutation test (given input R)
        nperm <- R + 1L
        permdist <- rep(0, nperm)
        permdist[1] <- Tstat
        
        # parallel or sequential computation?
        if(parallel){
          permdist[2:nperm] <- parSapply(cl = cl, X = integer(R), 
                                         FUN = Tperm.lm1, 
                                         xmat = x, yvec = y, 
                                         method = method, 
                                         homosced = homosced,
                                         exact = exact, xtx = xtx, 
                                         xtxi = xtxi, xinv = xinv, 
                                         ssy = ssy)
        } else {
          permdist[2:nperm] <- sapply(X = integer(R),
                                      FUN = Tperm.lm1,
                                      xmat = x, yvec = y, 
                                      method = method, 
                                      homosced = homosced,
                                      exact = exact, xtx = xtx, 
                                      xtxi = xtxi, xinv = xinv, 
                                      ssy = ssy)
        } # end if(parallel)
        
      } # end if(exact)
      
      ## permutation p-value
      p.value <- mean(permdist >= Tstat)
      
      ## correct coefficients
      if(!is.null(beta)) {
        ybar <- ybar + sum(xbar * beta)
        coefs <- coefs + beta
      }
      
      ## intercept
      alpha <- ybar - sum(xbar * coefs)
      coefs <- c(alpha, coefs)
      
      ## name coefficients
      xnames <- colnames(x)
      if(is.null(xnames)) xnames <- paste0("x", 1:p)
      names(coefs) <- c("(Intercept)", xnames)
      
    } else {
      
      ## MULTIVARIATE TEST
      
      ## observed test statistic
      Tuni <- Tstat.lm1.mv(x = x, y = y, homosced = homosced,
                           xtx = xtx, xtxi = xtxi, xinv = xinv, 
                           ssy = ssy, beta = coefs, combine = FALSE)
      Tstat <- max(Tuni)
      
      ## permutation distribution
      if(exact){
        
        # parallel or sequential computation?
        if(parallel){
          permdist <- parCapply(cl = cl, x = ix, 
                                FUN = Tperm.lm1.mv, 
                                xmat = x, ymat = y, 
                                method = method, 
                                homosced = homosced,
                                exact = exact, xtx = xtx, 
                                xtxi = xtxi, xinv = xinv, ssy = ssy)
        } else {
          permdist <- apply(X = ix, MARGIN = 2, 
                            FUN = Tperm.lm1.mv, 
                            xmat = x, ymat = y, 
                            method = method, 
                            homosced = homosced,
                            exact = exact, xtx = xtx, 
                            xtxi = xtxi, xinv = xinv, ssy = ssy)
        } # end if(parallel)
        
      } else {
        
        # approximate permutation test (given input R)
        nperm <- R + 1L
        permdist <- rep(0, nperm)
        permdist[1] <- Tstat
        
        ## parallel or sequential computation?
        if(parallel){
          permdist[2:nperm] <- parSapply(cl = cl, X = integer(R), 
                                         FUN = Tperm.lm1.mv, 
                                         xmat = x, ymat = y, 
                                         method = method, 
                                         homosced = homosced,
                                         exact = exact, xtx = xtx, 
                                         xtxi = xtxi, xinv = xinv, ssy = ssy)
        } else {
          permdist[2:nperm] <- sapply(X = integer(R),
                                      FUN = Tperm.lm1.mv,
                                      xmat = x, ymat = y, 
                                      method = method, 
                                      homosced = homosced,
                                      exact = exact, xtx = xtx, 
                                      xtxi = xtxi, xinv = xinv, ssy = ssy)
        } # end if(parallel)
        
      } # end if(exact)
      
      ## permutation p-value
      p.value <- mean(permdist >= Tstat)
      uni.p.value <- rep(NA, nvar)
      for(v in 1:nvar) uni.p.value[v] <- mean(permdist >= Tuni[v])
      
      ## correct coefficients
      if(!is.null(beta)) {
        ybar <- ybar + t(xbar) %*% beta
        coefs <- coefs + beta
      }
      
      ## intercept
      alpha <- ybar - t(xbar) %*% coefs
      coefs <- rbind(alpha, coefs)
      
      ## name coefficients
      xnames <- colnames(x)
      if(is.null(xnames)) xnames <- paste0("x", 1:p)
      rownames(coefs) <- c("(Intercept)", xnames)
      ynames <- colnames(y)
      if(is.null(ynames)) ynames <- paste0("y", 1:nvar)
      colnames(coefs) <- ynames
      
    } # end if(nvar == 1L)
    
    ### return results
    if(make.cl) stopCluster(cl)
    if(!perm.dist) permdist <- NULL
    res <- list(statistic = Tstat, p.value = p.value,
                perm.dist = permdist, method = method, 
                null.value = beta, homosced = homosced, 
                R = nperm - ifelse(exact, 0, 1), exact = exact,
                coefficients = coefs)
    if(nvar > 1L) {
      res$univariate <- Tuni
      res$adj.p.values <- uni.p.value
    }
    class(res) <- "rand.test.lm1"
    return(res)
    
  } # end rand.test.lm1.R


### permutation replication (univariate)
Tperm.lm1 <-
  function(i, xmat, yvec, method = "perm", 
           homosced = FALSE, exact = FALSE,
           xtx = NULL, xtxi = NULL, 
           xinv = NULL, ssy = NULL){
    n <- nrow(xmat)
    if(method == "perm"){
      if(!exact) i <- sample.int(n)
      yvec <- yvec[i]
    } else if(method == "flip"){
      if(!exact) i <- sample(c(-1, 1), size = n, replace = TRUE)
      yvec <- yvec * i
    } else {
      if(!exact) i <- c(sample.int(n), sample(c(-1, 1), size = n, replace = TRUE))
      yvec <- yvec * i[(n+1):(2*n)]   # flip sign
      yvec <- yvec[i[1:n]]            # permute
    }
    Tstat <- Tstat.lm1(x = xmat, y = yvec, homosced = homosced,
                       xtx = xtx, xtxi = xtxi, xinv = xinv, ssy = ssy)
    return(Tstat)
  } # end Tperm.lm1.R

### permutation replication (multivariate)
Tperm.lm1.mv <-
  function(i, xmat, ymat, method = "perm", 
           homosced = FALSE, exact = FALSE,
           xtx = NULL, xtxi = NULL, 
           xinv = NULL, ssy = NULL){
    n <- nrow(xmat)
    if(method == "perm"){
      if(!exact) i <- sample.int(n)
      ymat <- ymat[i,]
    } else if(method == "flip"){
      if(!exact) i <- sample(c(-1, 1), size = n, replace = TRUE)
      ymat <- ymat * i
    } else {
      if(!exact) i <- c(sample.int(n), sample(c(-1, 1), size = n, replace = TRUE))
      ymat <- ymat * i[(n+1):(2*n)]    # flip sign
      ymat <- ymat[i[1:n],]            # permute
    }
    Tstat <- Tstat.lm1.mv(x = xmat, y = ymat, homosced = homosced,
                          xtx = xtx, xtxi = xtxi, xinv = xinv, ssy = ssy)
    return(Tstat)
  } # end Tperm.lm1.mv.R


### test statistic (univariate) 
Tstat.lm1 <- 
  function(x, y, homosced = FALSE,
           xtx = NULL, xtxi = NULL, 
           xinv = NULL, ssy = NULL,
           beta = NULL){
    # assumes x and y are centered
    
    # check inputs
    p <- ncol(x)
    n <- length(y)
    
    # coefficients
    if(is.null(xtx)) xtx <- crossprod(x)
    if(is.null(xtxi)) xtxi <- solve(xtx)
    if(is.null(xinv)) xinv <- tcrossprod(xtxi, x)
    if(is.null(beta)) beta <- xinv %*% y
    
    # test statistic
    if(homosced){
      if(is.null(ssy)) ssy <- sum(y^2)
      rss <- ssy - sum(crossprod(x, y) * beta)
      sig <- rss / (n - p - 1)
      Tstat <- crossprod(beta, xtx %*% beta) / (p * sig)
    } else {
      omega <- crossprod(abs(y) * x)
      Tstat <- crossprod(beta, solve(xtxi %*% omega %*% xtxi) %*% beta)
    }
    
    as.numeric(Tstat)
    
  } # end Tstat.lm1.R

### test statistic (multivariate)
Tstat.lm1.mv <- 
  function(x, y, homosced = FALSE, xtx = NULL, 
           xtxi = NULL, xinv = NULL, ssy = NULL,
           beta = NULL, combine = TRUE){
    # assumes x and y are centered
    
    # check inputs
    p <- ncol(x)
    n <- nrow(y)
    nvar <- ncol(y)
    
    # coefficients
    if(is.null(xtx)) xtx <- crossprod(x)
    if(is.null(xtxi)) xtxi <- solve(xtx)
    if(is.null(xinv)) xinv <- tcrossprod(xtxi, x)
    if(is.null(beta)) beta <- xinv %*% y
    
    # test statistic
    Tstat <- rep(NA, nvar)
    if(homosced){
      if(is.null(ssy)) ssy <- colSums(y^2)
      rss <- ssy - colSums(crossprod(x, y) * beta)
      sig <- rss / (n - p - 1)
      for(v in 1:nvar) Tstat[v] <- (t(beta[,v]) %*% (xtx %*% beta[,v])) / (p * sig[v])
    } else {
      for(v in 1:nvar){
        omega <- crossprod(abs(y[,v]) * x)
        Tstat[v] <- t(beta[,v]) %*% (solve(xtxi %*% omega %*% xtxi) %*% beta[,v])
      }
    }
    
    if(combine) Tstat <- max(Tstat)
    as.numeric(Tstat)
    
  } # end Tstat.lm1.mv.R

