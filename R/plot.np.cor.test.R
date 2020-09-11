plot.np.cor.test <-
  function(x, alpha = 0.05, col = "grey", col.rr = "red",
           col.stat = "black", lty.stat = 2, lwd.stat = 2,
           xlab = "Test Statistic", main = "Permutation Distribution", 
           breaks = "scott", border = NA, box = TRUE, ...){
    
    if(!is.null(x$perm.dist)){
      
      # check "alpha"
      alpha <- as.numeric(alpha[1])
      if(alpha <= 0 | alpha >= 1) stop("Input 'alpha' must be between 0 and 1")
      
      # get breaks and midpoints
      h <- suppressWarnings(hist(x$perm.dist, breaks = breaks, plot = FALSE, ...))
      
      # define colors
      if(x$alternative == "two.sided"){
        cval <- quantile(x$perm.dist, probs = c(alpha / 2, 1 - alpha / 2))
        hcol <- ifelse(h$mids <= cval[1] | h$mids >= cval[2], col.rr, col)
      } else if(x$alternative == "less"){
        cval <- quantile(x$perm.dist, probs = alpha)
        hcol <- ifelse(h$mids <= cval, col.rr, col)
      } else if(x$alternative == "greater"){
        cval <- quantile(x$perm.dist, probs = 1 - alpha)
        hcol <- ifelse(h$mids >= cval, col.rr, col)
      }
      
      # plot perm.test
      hist(x$perm.dist, breaks = breaks, col = hcol, 
           border = border, main = main, xlab = xlab, ...)
      if(box) box()
      
      # plot statistic
      abline(v = x$statistic, lty = lty.stat, lwd = lwd.stat, col = col.stat)
      
    } else {
      
      warning("Plot unavailable because perm.dist is NULL\nSet perm.dist = TRUE when calling np.cor.test")
      
    }
    
  } # end plot.np.cor.test