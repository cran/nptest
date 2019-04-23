rand.test.cor <- 
  function(x, y, 
           alternative = c("two.sided", "less", "greater"),
           rho = 0, independent = FALSE,
           R = 9999, parallel = FALSE, cl = NULL,
           perm.dist = TRUE){
    # Randomization Tests of Correlation Coefficients
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: November 4, 2018
    
    
    #########   INITIAL CHECKS   #########
    
    ### check x and y
    x <- as.numeric(x)
    y <- as.numeric(y)
    n <- length(x)
    if(length(y) != n) stop("Inputs 'x' and 'y' must have same length.")
    
    ### check alternative
    alternative <- as.character(alternative[1])
    alternative <- pmatch(alternative, c("two.sided", "less", "greater"))
    if(is.na(alternative)) stop("Invalid 'alternative' input.")
    alternative <- c("two.sided", "less", "greater")[alternative]
    
    ### check rho
    rho <- as.numeric(rho[1])
    
    ### check independent
    independent <- as.logical(independent[1])
    
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
    
    
    #########   CORRELATION TEST   #########
    
    ### non-zero null hypothesis?
    xycor <- cor(x,y)
    if(!is.null(rho)) {
      beta <- rho * sd(x) / sd(y)
      y <- y - as.numeric(x * beta)
    }
    
    ### observed test statistic
    Tstat <- Tstat.cor(x, y, independent = independent)
    
    ### exact or approximate?
    suppressWarnings( nperm <- factorial(n) )
    exact <- ifelse(nperm <= R + 1L, TRUE, FALSE)
    
    ### permutation distribution
    if(exact){
      
      ## exact permutation test (R = factorial(n))
      ix <- permn(n)
      nperm <- ncol(ix)
      
      ## parallel or sequential computation?
      if(parallel){
        permdist <- parCapply(cl = cl, x = ix, 
                              FUN = Tperm.cor, 
                              xvec = x, yvec = y, 
                              independent = independent, 
                              exact = exact)
      } else {
        permdist <- apply(X = ix, MARGIN = 2, 
                          FUN = Tperm.cor,
                          xvec = x, yvec = y, 
                          independent = independent, 
                          exact = exact)
        # permdist <- rep(0, nperm)
        # for(j in 1:nperm) {
        #   permdist[j] <- Tstat.cor(x = x, y = y[ix[,j]], 
        #                            independent = independent, rho = rho)
        # }
      } # end if(parallel)
      
    } else {
      
      ## approximate permutation test (given input R)
      nperm <- R + 1L
      permdist <- rep(0, nperm)
      permdist[1] <- Tstat
      
      ## parallel or sequential computation?
      if(parallel){
        permdist[2:nperm] <- parSapply(cl = cl, X = integer(R), 
                                       FUN = Tperm.cor, 
                                       xvec = x, yvec = y, 
                                       independent = independent, 
                                       exact = exact)
      } else {
        permdist[2:nperm] <- sapply(X = integer(R),
                                    FUN = Tperm.cor,
                                    xvec = x, yvec = y,
                                    independent = independent,
                                    exact = exact)
         # for(j in 2:nperm) {
         #   permdist[j] <- Tstat.cor(x = x, y = sample(y),
         #                            independent = independent, rho = rho)
         # }
      } # end if(parallel)
      
    } # end if(exact)
    
    ### permutation p-value
    if(alternative == "less"){
      p.value <- mean(permdist <= Tstat)
    } else if(alternative == "greater"){
      p.value <- mean(permdist >= Tstat)
    } else {
      p.value <- mean(abs(permdist) >= abs(Tstat))
    }
    
    ### return results
    if(make.cl) stopCluster(cl)
    if(!perm.dist) permdist <- NULL
    res <- list(statistic = Tstat, p.value = p.value,
                perm.dist = permdist, alternative = alternative, 
                null.value = rho, independent = independent, 
                R = nperm - ifelse(exact, 0, 1), exact = exact,
                estimate = xycor)
    class(res) <- "rand.test.cor"
    return(res)
    
  } # end rand.test.cor.R


### test statistic 
Tstat.cor <- 
  function(x, y, independent = FALSE){
    x <- x - mean(x)
    y <- y - mean(y)
    ssx <- sum(x^2)
    ssy <- sum(y^2)
    ssxy <- sum(x * y)
    r <- ssxy / sqrt(ssx * ssy)
    if(independent){
      se <- sqrt( (1 - r^2) / (length(x) - 2) )
    } else {
      ssxy2 <- sum((x^2) * (y^2))
      se <- sqrt(ssxy2 / (ssx * ssy))
    }
    Tstat <- r / se
    return(Tstat)
  } # end Tstat.cor.R


### permutation replication
Tperm.cor <-
  function(i, xvec, yvec, independent = FALSE, exact = FALSE){
    if(exact){
      Tstat <- Tstat.cor(x = xvec, y = yvec[i],
                         independent = independent)
    } else {
      Tstat <- Tstat.cor(x = xvec, y = sample(yvec),
                         independent = independent)
    }
    return(Tstat)
  } # end Tperm.cor.R
