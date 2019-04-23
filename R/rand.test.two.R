rand.test.two <- 
  function(x, y,
           alternative = c("two.sided", "less", "greater"),
           mu = 0, var.equal = FALSE, median.test = FALSE,
           R = 9999, parallel = FALSE, cl = NULL,
           perm.dist = TRUE){
    # Two-Sample Randomization Test for Location (Mean/Median)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: November 4, 2018
    
    
    #########   INITIAL CHECKS   #########
    
    ### check x and y
    x <- as.matrix(x)
    y <- as.matrix(y)
    m <- nrow(x)
    n <- nrow(y)
    nvar <- ncol(x)
    if(ncol(y) != nvar) stop("Inputs 'x' and 'y' must have the same number of columns.")
    N <- m + n
    
    ### check alternative
    alternative <- as.character(alternative[1])
    alternative <- pmatch(alternative, c("two.sided", "less", "greater"))
    if(is.na(alternative)) stop("Invalid 'alternative' input.")
    alternative <- c("two.sided", "less", "greater")[alternative]
    
    ### check mu
    mu <- as.numeric(mu)
    if(length(mu) != nvar) mu <- rep(mu, length.out = nvar)
    
    ### check var.equal
    var.equal <- as.logical(var.equal[1])
    
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
    suppressWarnings( nperm <- choose(N, n) )
    exact <- ifelse(nperm <= R + 1L, TRUE, FALSE)
    if(exact){
      ix <- combn(N, n)
      nperm <- ncol(ix)
    }
    
    
    #########   LOCATION TEST   #########
    
    ### univariate or multivariate
    if(nvar == 1L){
      
      ## UNIVARIATE TEST
      
      ## non-zero null hypothesis?
      if(abs(mu) > 0){
        x <- x - mu
      }
      
      ## combine x and y
      z <- c(x, y)
      
      ## observed test statistic
      Tstat <- Tstat.two(x = x, y = y, var.equal = var.equal,
                         median.test = median.test)
      
      ## permutation distribution
      if(exact){
        
        # parallel or sequential computation?
        if(parallel){
          permdist <- parCapply(cl = cl, x = ix, 
                                FUN = Tperm.two, 
                                z = z, ny = n,
                                var.equal = var.equal,
                                median.test = median.test,
                                exact = exact)
        } else {
          permdist <- apply(X = ix, MARGIN = 2, 
                            FUN = Tperm.two, 
                            z = z, ny = n,
                            var.equal = var.equal,
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
                                         FUN = Tperm.two, 
                                         z = z, ny = n,
                                         var.equal = var.equal,
                                         median.test = median.test,
                                         exact = exact)
        } else {
          permdist[2:nperm] <- sapply(X = integer(R),
                                      FUN = Tperm.two, 
                                      z = z, ny = n,
                                      var.equal = var.equal,
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
        est <- median(outer(x, y, FUN = "-")) + mu
      } else {
        est <- mean(x) - mean(y) + mu
      }
      
      
    } else {
      
      ## MULTIVARIATE TEST
      
      ## non-zero null hypothesis?
      if(max(abs(mu)) > 0){
        for(v in 1:nvar) x[,v] <- x[,v] - mu[v]
      }
      
      ## combine x and y
      z <- rbind(x, y)
      
      ## observed test statistic
      Tuni <- Tstat.two.mv(x = x, y = y, var.equal = var.equal,
                           median.test = median.test, combine = FALSE)
      Tstat <- ifelse(alternative == "two.sided", max(abs(Tuni)),
                      ifelse(alternative == "greater", max(Tuni), min(Tuni)))
      
      ## permutation distribution
      if(exact){
        
        # parallel or sequential computation?
        if(parallel){
          permdist <- parCapply(cl = cl, x = ix, 
                                FUN = Tperm.two.mv, 
                                z = z, ny = n,
                                var.equal = var.equal,
                                median.test = median.test,
                                alternative = alternative, 
                                exact = exact)
        } else {
          permdist <- apply(X = ix, MARGIN = 2, 
                            FUN = Tperm.two.mv, 
                            z = z, ny = n,
                            var.equal = var.equal,
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
                                         FUN = Tperm.two.mv, 
                                         z = z, ny = n,
                                         var.equal = var.equal,
                                         median.test = median.test,
                                         alternative = alternative, 
                                         exact = exact)
        } else {
          permdist[2:nperm] <- sapply(X = integer(R),
                                      FUN = Tperm.two.mv, 
                                      z = z, ny = n,
                                      var.equal = var.equal,
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
        for(v in 1:nvar) est[v] <- median(outer(x[,v], y[,v], FUN = "-")) + mu[v]
      } else {
        for(v in 1:nvar) est[v] <- mean(x[,v]) - mean(y[,v]) + mu[v]
      }
      
    } # end if(nvar == 1L)
    
    ### return results
    if(make.cl) stopCluster(cl)
    if(!perm.dist) permdist <- NULL
    res <- list(statistic = Tstat, p.value = p.value,
                perm.dist = permdist, alternative = alternative, 
                null.value = mu, var.equal = var.equal, 
                median.test = median.test, 
                R = nperm - ifelse(exact, 0, 1), exact = exact,
                estimate = est)
    if(nvar > 1L) {
      res$univariate <- Tuni
      res$adj.p.values <- uni.p.value
    }
    class(res) <- "rand.test.two"
    return(res)
    
  } # end rand.test.two.R


### permutation replication (univariate)
Tperm.two <-
  function(i, z, ny, var.equal = FALSE, 
           median.test = FALSE, exact = FALSE){
    if(!exact) i <- sample.int(n = length(z), size = ny)
    Tstat <- Tstat.two(x = z[-i], y = z[i],
                       var.equal = var.equal, 
                       median.test = median.test)
    return(Tstat)
  } # end Tperm.two.R

### permutation replication (multivariate)
Tperm.two.mv <-
  function(i, z, ny, var.equal = FALSE, median.test = FALSE, 
           alternative = "two.sided", combine = TRUE, exact = FALSE){
    if(!exact) i <- sample.int(n = nrow(z), size = ny)
    Tstat <- Tstat.two.mv(x = z[-i,], y = z[i,],
                         var.equal = var.equal, 
                         median.test = median.test,
                         alternative = alternative,
                         combine = combine)
    return(Tstat)
  } # end Tperm.two.mv.R

### test statistic (univariate)
Tstat.two <- 
  function(x, y, var.equal = FALSE, median.test = FALSE){
    m <- length(x)
    n <- length(y)
    N <- m + n
    if(median.test){
      mn <- m * n
      rz <- rank(c(x, y))
      U <- (sum(rz[1:m]) - m * (m + 1) / 2) / mn
      if(var.equal){
        se <- sqrt( (N + 1) / (12 * mn) )
      } else {
        rx <- (rz[1:m] - rank(x)) / n
        ry <- (rz[(m + 1):N] - rank(y)) / m
        xi.x <- sum((rx - mean(rx))^2) / (m - 1)
        xi.y <- sum((ry - mean(ry))^2) / (n - 1)
        se <- sqrt( (xi.x / m) + (xi.y / n) )
      } # end if(var.equal)
      Tstat <- (U - 1/2) / se
    } else {
      mux <- mean(x)
      muy <- mean(y)
      if(var.equal){
        ssx <- sum((x - mux)^2)
        ssy <- sum((y - muy)^2)
        sigsq <- (ssx + ssy) / (N - 2)
        se <- sqrt( sigsq * ((1 / m) + (1 / n)) )
      } else {
        vx <- sum((x - mux)^2) / (m - 1)
        vy <- sum((y - muy)^2) / (n - 1)
        se <- sqrt( (vx / m) + (vy / n) )
      } # end if(var.equal)
      Tstat <- (mux - muy) / se
    } # end if(median.test)
    Tstat
  } # end Tstat.two.R

### test statistic (multivariate)
Tstat.two.mv <- 
  function(x, y, var.equal = FALSE, median.test = FALSE,
           alternative = "two.sided", combine = TRUE){
    m <- nrow(x)
    n <- nrow(y)
    nvar <- ncol(x)
    N <- m + n
    if(median.test){
      mn <- m * n
      U <- se <- rep(0, nvar)
      for(v in 1:nvar){
        rz <- rank(c(x[,v], y[,v]))
        U[v] <- (sum(rz[1:m]) - m * (m + 1) / 2) / mn
        if(var.equal){
          se[v] <- sqrt( (N + 1) / (12 * mn) )
        } else {
          rx <- (rz[1:m] - rank(x[,v])) / n
          ry <- (rz[(m + 1):N] - rank(y[,v])) / m
          xi.x <- sum((rx - mean(rx))^2) / (m - 1)
          xi.y <- sum((ry - mean(ry))^2) / (n - 1)
          se[v] <- sqrt( (xi.x / m) + (xi.y / n) )
        } # end if(var.equal)
      }
      Tstat <- (U - 1/2) / se
    } else {
      mux <- colMeans(x)
      muy <- colMeans(y)
      if(var.equal){
        ssx <- ssy <- rep(0, nvar)
        for(v in 1:nvar){
          ssx[v] <- sum((x[,v] - mux[v])^2)
          ssy[v] <- sum((y[,v] - muy[v])^2)
        }
        sigsq <- sum(ssx + ssy) / (N - 2)
        se <- sqrt( sigsq * ((1 / m) + (1 / n)) )
      } else {
        vx <- vy <- rep(0, nvar)
        for(v in 1:nvar){
          vx[v] <- sum((x[,v] - mux[v])^2) / (m - 1)
          vy[v] <- sum((y[,v] - muy[v])^2) / (n - 1)
        }
        se <- sqrt( (vx / m) + (vy / n) )
      } # end if(var.equal)
      Tstat <- (mux - muy) / se
    } # end if(median.test)
    if(combine) {
      Tstat <- ifelse(alternative == "two.sided", max(abs(Tstat)),
                      ifelse(alternative == "greater", max(Tstat), min(Tstat)))
    }
    as.numeric(Tstat)
  } # end Tstat.two.mv.R
