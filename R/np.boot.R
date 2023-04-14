np.boot <- 
  function(x, statistic, ..., R = 9999, level = c(0.9, 0.95, 0.99),
           method = c("norm", "basic", "perc", "stud", "bca")[-4], 
           sdfun = NULL, sdrep = 99, jackknife = NULL, 
           parallel = FALSE, cl = NULL, boot.dist = TRUE){
    # Nonparametric Bootstrap (SE, Bias, and CIs)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: 2023-04-14
    
    
    ######***######   CHECKS   ######***######
    
    ### check x
    x <- as.vector(x)
    nobs <- length(x)
    
    ### check statistic
    statistic <- as.function(statistic)
    t0 <- statistic(x, ...)
    if(!is.vector(t0))
      stop("Output of the 'statistic' function must be a vector")
    nstat <- length(t0)
    varnames <- names(t0)
    
    ### check R
    R <- as.integer(R[1])
    if(R < 1L) stop("Input 'R' must be a positive integer")
    nboot <- R + 1L
    
    ### check level
    level <- as.numeric(level)
    if(any(level <= 0) | any(level >= 1)) 
      stop("Input 'level' must be between 0 and 1 for each element")
    nlevel <- length(level)
    alphas <- 1 - level
    
    ### check method
    if(!is.null(method)){
      method <- as.character(method)
      method <- pmatch(method, c("norm", "perc", "basic", "stud", "bca"))
      if(any(is.na(method))) stop("Input 'method' contains an invalid method.")
      method <- c("norm", "perc", "basic", "stud", "bca")[method]
    }
    
    ### check sdfun and sdrep
    make.sdfun <- FALSE
    if(any(method == "stud")){
      if(!is.null(sdfun)){
        sdfun <- as.function(sdfun)
        sd0 <- sdfun(x, ...)
        if(!is.vector(sd0))
          stop("Output of the 'sdfun' function must be a vector")
        if(length(sd0) != nstat) stop("Outputs of 'statistic' and 'sdfun' functions must have the same length")
      } else {
        make.sdfun <- TRUE
        sdrep <- as.integer(sdrep[1])
        if(sdrep < 1L) stop("Input 'sdrep' must be a positive integer")
        nbootsd <- sdrep + 1L
        sdfun <- function(x, ...){
          #bootx <- matrix(sample(x, size = nobs * sdrep, replace = TRUE), nrow = sdrep, ncol = nobs)
          bootx <- matrix(x[sdfun.sample.index], nrow = sdrep, ncol = nobs)
          if(nstat == 1L){
            return(sd(apply(X = rbind(x, bootx), MARGIN = 1, FUN = statistic, ...)))
          } else {
            return(apply(apply(X = rbind(x, bootx), MARGIN = 1, FUN = statistic, ...), MARGIN = 1, FUN = sd))
          }
        } # end sdfun
      }
    } # end if(any(method == "stud"))
    
    ### check jackknife
    make.jack <- FALSE
    if(any(method == "bca")){
      if(is.null(jackknife)){
        make.jack <- TRUE
        jackknife <- statistic
      } else {
        jackknife <- as.function(jackknife)
        j0 <- jackknife(x, ...)
        if(length(j0) != nstat) stop("Outputs of 'statistic' and 'jackknife' functions must have the same length")
      }
    }
    
    ### check parallel and cl inputs
    parallel <- as.logical(parallel[1])
    make.cl <- FALSE
    if(parallel){
      if(is.null(cl)){
        make.cl <- TRUE
        cl <- parallel::makeCluster(2L)
      } else {
        if(!any(class(cl) == "cluster")) stop("Input 'cl' must be an object of class 'cluster'.")
      }
    }
    
    
    ######***######   BOOTSTRAP   ######***######
    
    ### bootstrap x
    bootx <- matrix(sample(x, size = nobs * R, replace = TRUE), nrow = R, ncol = nobs)
    
    ### bootstrap distribution
    bootdist <- matrix(0, nboot, nstat)
    colnames(bootdist) <- varnames
    bootdist[1,] <- t0
    if(parallel){
      bootdist[2:nboot,] <- t(parallel::parApply(cl = cl, X = bootx,
                                                 MARGIN = 1, FUN = statistic, ...))
    } else {
      bootdist[2:nboot,] <- t(apply(X = bootx, MARGIN = 1, FUN = statistic, ...))
    }
    if(!any(method == "stud")) rm(bootx)
    
    ### bootstrap standard error
    bootse <- apply(X = bootdist, MARGIN = 2, FUN = sd)
    
    ### bootstrap bias
    bootbias <- apply(X = bootdist, MARGIN = 2, FUN = mean) - t0
    
    ### bootstrap covariance?
    if(nstat > 1L) {
      bootcov <- cov(bootdist)
    } else {
      bootcov <- NULL
    }
    
    
    ######***######   INTERVALS   ######***######
    
    ### define dimensions
    dims <- c(nlevel, 2, nstat)
    dnames <- list(paste0(as.character(level * 100), "%"), 
                   c("lower", "upper"), varnames)
    
    ### normal interval
    normal <- NULL
    if(any(method == "norm")){
      normal <- array(NA, dim = dims, dimnames = dnames)
      for(j in 1:nstat){
        normal[,1,j] <- t0[j] - bootbias[j] - qnorm(1 - alphas/2) * bootse[j]
        normal[,2,j] <- t0[j] - bootbias[j] - qnorm(alphas/2) * bootse[j]
      }
      normal <- normal[,,,drop=TRUE]
    } # end if(any(method == "norm"))
    
    ### percent interval
    percent <- NULL
    if(any(method == "perc")){
      percent <- array(NA, dim = dims, dimnames = dnames)
      probs <- c(alphas/2, 1 - alphas/2)
      for(j in 1:nstat){
        quant <- quantile(bootdist[,j], probs = probs, na.rm = TRUE)
        percent[,1,j] <- quant[1:nlevel]
        percent[,2,j] <- quant[(nlevel+1):length(probs)]
      }
      percent <- percent[,,,drop=TRUE]
    } # end if(any(method == "perc"))
    
    ### basic interval
    basic <- NULL
    if(any(method == "basic")){
      basic <- array(NA, dim = dims, dimnames = dnames)
      probs <- c(alphas/2, 1 - alphas/2)
      for(j in 1:nstat){
        quant <- quantile(bootdist[,j], probs = probs, na.rm = TRUE)
        basic[,1,j] <- 2 * t0[j] - quant[(nlevel+1):length(probs)]
        basic[,2,j] <- 2 * t0[j] - quant[1:nlevel]
      }
      basic <- basic[,,,drop=TRUE]
    } # end if(any(method == "basic"))
    
    ### student interval
    student <- NULL
    if(any(method == "stud")){
      student <- array(NA, dim = dims, dimnames = dnames)
      bootdistse <- matrix(0, nboot, nstat)
      sdfun.sample.index <- sample.int(nobs, size = nobs * sdrep, replace = TRUE)
      if(parallel){
        bootdistse[1:nboot,] <- t(parallel::parApply(cl = cl, X = rbind(x, bootx),
                                                     MARGIN = 1, FUN = sdfun, ...))
      } else {
        bootdistse[1:nboot,] <- t(apply(X = rbind(x, bootx), MARGIN = 1, FUN = sdfun, ...))
      }
      for(j in 1:nstat){
        tstat <- (bootdist[,j] - t0[j]) / bootdistse[,j]
        student[,1,j] <- t0[j] - quantile(tstat, probs = 1 - alphas/2, na.rm = TRUE) * bootse[j]
        student[,2,j] <- t0[j] - quantile(tstat, probs = alphas/2, na.rm = TRUE) * bootse[j]
      }
      student <- student[,,,drop=TRUE]
    } # end if(any(method == "stud"))
    
    ### bca interval
    bca <- z0 <- acc <- NULL
    if(any(method == "bca")){
      bca <- array(NA, dim = dims, dimnames = dnames)
      z1 <- qnorm(alphas/2)
      z2 <- qnorm(1 - alphas/2)
      jackstat <- matrix(NA, nobs, nstat)
      for(i in 1:nobs) jackstat[i,] <- jackknife(x[-i], ...)
      meanjack <- colMeans(jackstat)
      z0 <- acc <- rep(0, nstat)
      for(j in 1:nstat){
        z0[j] <- qnorm(mean(bootdist[,j] < t0[j]))
        acc[j] <- sum((meanjack[j] - jackstat[,j])^3) / (6 * sum((meanjack[j] - jackstat[,j])^2)^(3/2))
        alphas1 <- pnorm(z0[j] + (z0[j] + z1) / (1 - acc[j]*(z0[j] + z1)))
        alphas2 <- pnorm(z0[j] + (z0[j] + z2) / (1 - acc[j]*(z0[j] + z2)))
        probs <- c(alphas1, alphas2)
        quant <- quantile(bootdist[,j], probs = probs, na.rm = TRUE)
        bca[,1,j] <- quant[1:nlevel]
        bca[,2,j] <- quant[(nlevel+1):length(probs)]
      }
      bca <- bca[,,,drop=TRUE]
    } # end if(any(method == "bca"))
    
    
    ######***######   RESULTS   ######***######
    
    ### clean up
    if(make.cl) parallel::stopCluster(cl)
    if(make.sdfun) sdfun <- NULL
    if(make.jack) jackknife <- NULL
    if(!boot.dist) {
      bootdist <- NULL
    } else {
      if(nstat == 1L) bootdist <- c(bootdist)
    }
    
    ### return result
    res <- list(t0 = t0, se = bootse, bias = bootbias, cov = bootcov,
                normal = normal, basic = basic, percent = percent,
                student = student, bca = bca, z0 = z0, acc = acc,
                boot.dist = bootdist, statistic = statistic, 
                R = R, level = level, sdfun = sdfun, sdrep = sdrep,
                jackknife = jackknife)
    class(res) <- "np.boot"
    return(res)
    
} # end np.boot