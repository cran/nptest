rand.test.lm2 <-
  function(x, y, z, 
           method = c("HJ", "KC", "SW", "TB",
                      "FL", "MA", "OS", "DS"),
           beta = NULL, homosced = FALSE, lambda = 0,
           R = 9999, parallel = FALSE, cl = NULL,
           perm.dist = TRUE){
    # Randomization Test for Regression (w/ covariates)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: 2023-04-14
    
    # Methods with Z:
    # ---  permute X
    # DS:  Y = P %*% X %*% beta + Z %*% gamma + e
    # OS:  Y = P %*% Rz %*% X %*% beta + Z %*% gamma + e
    # ---  permute Y
    # MA:  P %*% Y = X %*% beta + Z %*% gamma + e
    # FL:  (P %*% Rz + Hz) %*% Y = X %*% beta + Z %*% gamma + e
    # TB:  (P %*% Rm + Hm) %*% Y = X %*% beta + Z %*% gamma + e
    
    # Methods without Z (all permute Y):
    # SW:  P %*% Rz %*% Y = X %*% beta + e
    # KC:  P %*% Rz %*% Y = Rz %*% X %*% beta + e
    # HJ:  P %*% t(Q) %*% Rz %*% Y = t(Q) %*% Rz %*% X %*% beta + e
    
    
    #########   INITIAL CHECKS   #########
    
    ### check x, y, and z
    x <- as.matrix(x)
    z <- as.matrix(z)
    y <- as.matrix(y)
    nvar <- ncol(y)
    n <- nrow(x)
    p <- ncol(x)
    q <- ncol(z)
    r <- p + q
    if(nrow(y) != n) stop("Inputs 'x' and 'y' are incompatible.")
    if(nrow(z) != n) stop("Inputs 'x' and 'z' are incompatible.")
    
    ### check method
    method <- as.character(method[1])
    method.options <- c("DS", "OS", "MA", "FL",
                        "TB", "SW", "KC", "HJ")
    method <- method.options[pmatch(method, method.options)]
    if(is.na(method)) stop("Invalid 'method' input.")
    use.z <- ifelse(any(method == c("DS", "OS", "MA", "FL", "TB")), TRUE, FALSE)
    
    ### check beta
    if(!is.null(beta)){
      beta <- as.matrix(beta)
      if(nrow(beta) != p) stop("Inputs 'beta' and 'x' are incompatible.")
      if(ncol(beta) != nvar) stop("Inputs 'beta' and 'y' are incompatible.")
    }
    
    ### check homosced
    homosced <- as.logical(homosced[1])
    
    ### check lambda
    nlambda <- length(lambda)
    if(nlambda == 1L){
      lambda <- rep(lambda, r)
    } else if(nlambda != r){
      warning("length(lambda) != ncol(cbind(x,z))\nUsing lambda <- rep(lambda[1], ncol(cbind(x,z)))")
      lambda <- rep(lambda[1], r)
    }
    if(any(lambda < 0)) stop("Input lambda must contain non-negative penalty weights.")
    lambda <- n * lambda
    lambda.x <- lambda[1:p]
    lambda.z <- lambda[(p+1):r]
    
    ### check R
    R <- as.integer(R)
    if(R < 0) {
      stop("Input 'R' must be a non-negative integer.")
    } else if(R == 0 && nvar > 1L){
      stop("Input 'R' must be a positive integer for multivariate tests.")
    }
    
    ### check parallel
    parallel <- as.logical(parallel[1])
    
    ### check 'cl'
    make.cl <- FALSE
    if(parallel){
      if(is.null(cl)){
        make.cl <- TRUE
        cl <- parallel::makeCluster(2L)
      } else {
        if(!any(class(cl) == "cluster")) stop("Input 'cl' must be an object of class 'cluster'.")
      }
    }
    
    ### exact or approximate?
    suppressWarnings( nperm <- factorial(n) )
    exact <- ifelse(nperm <= R + 1L, TRUE, FALSE)
    if(exact){
      ix <- permn(n)
      nperm <- ncol(ix)
    }
    
    
    #########   REGRESSION TEST   #########
    
    ### estimate coefficients
    m <- cbind(1, x, z)
    coefs <- psdinv(crossprod(m) + diag(c(0, lambda), nrow = r + 1, ncol = r + 1)) %*% crossprod(m, y)
    xnames <- colnames(x)
    if(is.null(xnames)) xnames <- paste0("x", 1:p)
    znames <- colnames(z)
    if(is.null(znames)) znames <- paste0("z", 1:q)
    if(nvar == 1L){
      coefs <- as.numeric(coefs)
      names(coefs) <- c("(Intercept)", xnames, znames)
    } else {
      ynames <- colnames(y)
      if(is.null(ynames)) ynames <- paste0("y", 1:nvar)
      rownames(coefs) <- c("(Intercept)", xnames, znames)
      colnames(coefs) <- ynames
    }
    
    ### non-zero null hypothesis?
    if(!is.null(beta)) y <- y - x %*% beta
    
    ### center data
    xbar <- colMeans(x)
    ybar <- colMeans(y)
    zbar <- colMeans(z)
    x <- scale(x, center = xbar, scale = FALSE)
    y <- scale(y, center = ybar, scale = FALSE)
    z <- scale(z, center = zbar, scale = FALSE)
    
    ### process data (if needed)
    if(method == "OS"){
      x <- x - z %*% psdinv(crossprod(z) + diag(lambda.z, nrow = q, ncol = q)) %*% crossprod(z, x)
    } else if(method == "FL"){
      fit <- z %*% psdinv(crossprod(z) + diag(lambda.z, nrow = q, ncol = q)) %*% crossprod(z, y)
      res <- y - fit
    } else if(method == "TB"){
      m <- cbind(x, z)
      mtm <- crossprod(m) + diag(lambda, nrow = r, ncol = r)
      mtmi <- psdinv(mtm)
      minv <- tcrossprod(mtmi, m)
      coefs0 <- minv %*% y
      fit <- m %*% coefs0
      res <- y - fit
    } else if(method == "SW"){
      y <- y - z %*% psdinv(crossprod(z) + diag(lambda.z, nrow = q, ncol = q)) %*% crossprod(z, y)
    } else if(method == "KC"){
      zinv <- tcrossprod(psdinv(crossprod(z) + diag(lambda.z, nrow = q, ncol = q)), z)
      x <- x - z %*% (zinv %*% x)
      y <- y - z %*% (zinv %*% y)
    } else if(method == "HJ"){
      if(max(lambda.z) == 0){
        zsvd <- svd(z, nu = n, nv = 0)
      } else {
        ztzeig <- eigen(crossprod(z) + diag(lambda.z, nrow = q, ncol = q))
        ztzisqrt <- ztzeig$vectors %*% diag(1/sqrt(ztzeig$values), nrow = q, ncol = q)
        zsvd <- svd(z %*% ztzisqrt, nu = n, nv = 0)
      }
      x <- crossprod(zsvd$u[, (q+1):n, drop = FALSE], x)
      y <- crossprod(zsvd$u[, (q+1):n, drop = FALSE], y)
      n <- n - q
    }
    
    ### set z to NULL (if no longer needed)
    if(!use.z) {
      z <- NULL
      lambda <- lambda[1:p]
      r <- p
    }
    
    ### make m and crossprod matrices
    if(method != "TB"){
      m <- cbind(x, z)
      mtm <- crossprod(m) + diag(lambda, nrow = r, ncol = r)
      mtmi <- psdinv(mtm)
      minv <- tcrossprod(mtmi, m)
    }
    
    ### univariate or multivariate?
    if(nvar == 1L){
      
      ## UNIVARIATE TEST
      
      y <- as.numeric(y)
      
      ## observed test statistic
      if(use.z){
        Tstat <- Tstat.lm2(m = m, y = y, homosced = homosced,
                           mtm = mtm, mtmi = mtmi, minv = minv, p = p)
      } else {
        Tstat <- Tstat.lm1(x = m, y = y, homosced = homosced,
                           xtx = mtm, xtxi = mtmi, xinv = minv)
      }
      
      ## further processing (if needed)
      if(any(method == c("FL", "TB"))) {
        if(method == "TB") {
          xbeta <- x %*% coefs0[1:p]
          fit <- fit - xbeta
          res <- res + xbeta
        }
        y <- cbind(fit, res)
      }
      
      ## permutation distribution
      if(R == 0){
        # parametric test
        if(homosced){
          errdf <- ifelse(method == "HJ", n - p - 1, n - p - q - 1)
          if(method %in% c("KC", "SW")){
            Tstat <- Tstat * errdf / (n - p - 1)
          } 
          p.value <- 1 - pf(Tstat, df1 = p, df2 = errdf)
        } else {
          p.value <- 1 - pchisq(Tstat, df = p)
        }
        method <- "parametric"
        perm.dist <- FALSE
      } else if(exact){
        
        # parallel or sequential computation?
        if(parallel){
          permdist <- parallel::parCapply(cl = cl, x = ix, 
                                          FUN = Tperm.lm2, 
                                          m = m, y = y, p = p,
                                          method = method, 
                                          homosced = homosced, exact = exact, 
                                          lambda = lambda, use.z = use.z,
                                          mtm = mtm, mtmi = mtmi, minv = minv)
        } else {
          permdist <- apply(X = ix, MARGIN = 2, 
                            FUN = Tperm.lm2, 
                            m = m, y = y, p = p,
                            method = method, 
                            homosced = homosced, exact = exact, 
                            lambda = lambda, use.z = use.z,
                            mtm = mtm, mtmi = mtmi, minv = minv)
        } # end if(parallel)
        
      } else {
        
        # approximate permutation test (given input R)
        nperm <- R + 1L
        permdist <- rep(0, nperm)
        permdist[1] <- Tstat
        
        # parallel or sequential computation?
        if(parallel){
          permdist[2:nperm] <- parallel::parSapply(cl = cl, X = integer(R), 
                                                   FUN = Tperm.lm2, 
                                                   m = m, y = y, p = p,
                                                   method = method, 
                                                   homosced = homosced, exact = exact, 
                                                   lambda = lambda, use.z = use.z,
                                                   mtm = mtm, mtmi = mtmi, minv = minv)
        } else {
          permdist[2:nperm] <- sapply(X = integer(R),
                                      FUN = Tperm.lm2, 
                                      m = m, y = y, p = p,
                                      method = method, 
                                      homosced = homosced, exact = exact, 
                                      lambda = lambda, use.z = use.z,
                                      mtm = mtm, mtmi = mtmi, minv = minv)
        } # end if(parallel)
        
      } # end if(exact)
      
      ## permutation p-value
      if(R > 0) p.value <- mean(permdist >= Tstat)
      
    } else {
      
      ## observed test statistic
      if(use.z){
        Tuni <- Tstat.lm2.mv(m = m, y = y, homosced = homosced,
                             mtm = mtm, mtmi = mtmi, minv = minv, 
                             p = p, combine = FALSE)
      } else {
        Tuni <- Tstat.lm1.mv(x = m, y = y, homosced = homosced,
                             xtx = mtm, xtxi = mtmi, xinv = minv,
                             combine = FALSE)
      }
      Tstat <- max(Tuni)
      
      ## further processing (if needed)
      if(any(method == c("FL", "TB"))) {
        if(method == "TB") {
          xbeta <- x %*% coefs0[1:p,,drop=FALSE]
          fit <- fit - xbeta
          res <- res + xbeta
        }
        y <- array(NA, dim = c(n, nvar, 2))
        y[,,1] <- fit
        y[,,2] <- res
      }
      
      ## permutation distribution
      if(exact){
        
        # parallel or sequential computation?
        if(parallel){
          permdist <- parallel::parCapply(cl = cl, x = ix, 
                                          FUN = Tperm.lm2.mv, 
                                          m = m, y = y, method = method, 
                                          homosced = homosced, exact = exact,
                                          lambda = lambda, mtm = mtm, mtmi = mtmi, 
                                          minv = minv, use.z = use.z, p = p)
        } else {
          permdist <- apply(X = ix, MARGIN = 2, 
                            FUN = Tperm.lm2.mv, 
                            m = m, y = y, method = method, 
                            homosced = homosced, exact = exact, 
                            lambda = lambda, mtm = mtm, mtmi = mtmi, 
                            minv = minv, use.z = use.z, p = p)
        } # end if(parallel)
        
      } else {
        
        # approximate permutation test (given input R)
        nperm <- R + 1L
        permdist <- rep(0, nperm)
        permdist[1] <- Tstat
        
        # parallel or sequential computation?
        if(parallel){
          permdist[2:nperm] <- parallel::parSapply(cl = cl, X = integer(R), 
                                                   FUN = Tperm.lm2.mv, 
                                                   m = m, y = y, method = method, 
                                                   homosced = homosced, exact = exact, 
                                                   lambda = lambda, mtm = mtm, mtmi = mtmi, 
                                                   minv = minv, use.z = use.z, p = p)
        } else {
          permdist[2:nperm] <- sapply(X = integer(R),
                                      FUN = Tperm.lm2.mv, 
                                      m = m, y = y, method = method, 
                                      homosced = homosced, exact = exact, 
                                      lambda = lambda, mtm = mtm, mtmi = mtmi, 
                                      minv = minv, use.z = use.z, p = p)
        } # end if(parallel)
        
      } # end if(exact)
      
      ## permutation p-value
      p.value <- mean(permdist >= Tstat)
      uni.p.value <- rep(NA, nvar)
      for(v in 1:nvar) uni.p.value[v] <- mean(permdist >= Tuni[v])
      
      
    } # end if(nvar == 1L)
    
    ### return results
    if(make.cl) parallel::stopCluster(cl)
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
    class(res) <- "rand.test.lm2"
    return(res)
    
  } # end rand.test.lm2.R


### permutation replication (univariate)
Tperm.lm2 <-
  function(i, m, y, p = 1, method = "HJ", 
           homosced = FALSE, exact = FALSE, 
           lambda = 0, use.z = NULL, mtm = NULL, 
           mtmi = NULL, minv = NULL){
    if(!exact) i <- sample.int(nrow(m))
    if(any(method == c("DS", "OS"))){
      m[,1:p] <- m[i, 1:p, drop = FALSE]
      mtm <- crossprod(m) + diag(lambda, nrow = ncol(m))
      mtmi <- psdinv(mtm)
      minv <- tcrossprod(mtmi, m)
    } else if(any(method == c("FL", "TB"))) {
      y <- y[,1] + y[i,2]
    } else {
      y <- y[i]
    }
    if(is.null(use.z)) use.z <- ifelse(any(method == c("DS", "OS", "MA", "FL", "TB")), TRUE, FALSE)
    if(use.z){
      Tstat <- Tstat.lm2(m = m, y = y, homosced = homosced,
                         mtm = mtm, mtmi = mtmi, minv = minv, p = p)
    } else {
      Tstat <- Tstat.lm1(x = m, y = y, homosced = homosced,
                         xtx = mtm, xtxi = mtmi, xinv = minv)
    }
    return(Tstat)
  } # end Tperm.lm2.R

### permutation replication (multivariate)
Tperm.lm2.mv <-
  function(i, m, y, 
           method = "HJ", homosced = FALSE, exact = FALSE, 
           lambda = 0, mtm = NULL, mtmi = NULL, minv = NULL, 
           use.z = NULL, p = 1){
    if(!exact) i <- sample.int(nrow(m))
    if(any(method == c("DS", "OS"))){
      m[,1:p] <- m[i, 1:p, drop = FALSE]
      mtm <- crossprod(m) + diag(lambda, nrow = ncol(m))
      mtmi <- psdinv(mtm)
      minv <- tcrossprod(mtmi, m)
    } else if(any(method == c("FL", "TB"))) {
      y <- as.matrix(y[,,1] + y[i,,2])
    } else {
      y <- y[i,,drop=FALSE]
    }
    if(is.null(use.z)) use.z <- ifelse(any(method == c("DS", "OS", "MA", "FL", "TB")), TRUE, FALSE)
    if(use.z){
      Tstat <- Tstat.lm2.mv(m = m, y = y, homosced = homosced,
                            mtm = mtm, mtmi = mtmi, minv = minv, p = p)
    } else {
      Tstat <- Tstat.lm1.mv(x = m, y = y, homosced = homosced,
                            xtx = mtm, xtxi = mtmi, xinv = minv)
    }
    return(Tstat)
  } # end Tperm.lm2.mv.R

## test statistic (univariate)
Tstat.lm2 <- function(m, y, p = 1, homosced = FALSE,
                      mtm = NULL, mtmi = NULL, minv = NULL){
  # assumes m and y are centered
  
  # check inputs
  r <- ncol(m)
  n <- length(y)
  
  # coefficients
  if(is.null(mtm)) mtm <- crossprod(m)
  if(is.null(mtmi)) mtmi <- psdinv(mtm)
  if(is.null(minv)) minv <- tcrossprod(mtmi, m)
  theta <- minv %*% y
  
  # fitted values and residuals
  fit <- as.numeric(m %*% theta)
  res <- y - fit
  
  # test statistic
  if(homosced){
    top <- t(theta[1:p]) %*% psdinv(mtmi[1:p, 1:p, drop = FALSE]) %*% theta[1:p]
    sigsq <- sum(res^2) / (n - r - 1)
    Tstat <- top / (p * sigsq)
  } else {
    omega <- crossprod(abs(res) * m)
    siosi <- mtmi %*% omega %*% mtmi
    Tstat <- t(theta[1:p]) %*% psdinv(siosi[1:p, 1:p, drop = FALSE]) %*% theta[1:p]
  } # end if(independent)
  
  as.numeric(Tstat)
  
} # end Tstat.lm2.R

## test statistic (multivariate)
Tstat.lm2.mv <- 
  function(m, y, homosced = FALSE, mtm = NULL, 
           mtmi = NULL, minv = NULL, 
           p = 1, combine = TRUE){
    # assumes m and y are centered
    
    # check inputs
    r <- ncol(m)
    n <- nrow(y)
    nvar <- ncol(y)
    
    # coefficients
    if(is.null(mtm)) mtm <- crossprod(m)
    if(is.null(mtmi)) mtmi <- psdinv(mtm)
    if(is.null(minv)) minv <- tcrossprod(mtmi, m)
    theta <- minv %*% y
    
    # fitted values and residuals
    fit <- m %*% theta
    res <- y - fit
    
    # test statistic
    Tstat <- rep(NA, nvar)
    if(homosced){
      mtmi.xi <- psdinv(mtmi[1:p, 1:p, drop = FALSE])
      sig <- colSums(res^2) / (n - r - 1)
      for(v in 1:nvar){
        Tstat[v] <- (t(theta[1:p,v]) %*% mtmi.xi %*% theta[1:p,v]) / (p * sig[v])
      }
    } else {
      for(v in 1:nvar){
        omega <- crossprod(abs(res[,v]) * m)
        siosi <- mtmi %*% omega %*% mtmi
        Tstat[v] <- t(theta[1:p,v]) %*% psdinv(siosi[1:p, 1:p, drop = FALSE]) %*% theta[1:p,v]
      }
    } # end if(independent)
    
    if(combine) Tstat <- max(Tstat)
    as.numeric(Tstat)
    
  } # end Tstat.lm2.mv.R

