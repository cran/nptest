mcse <- 
  function(R, delta, conf.level = 0.95, sig.level = 0.05,
           alternative = c("two.sided", "one.sided")){
    # Monte Carlo Standard Errors for Nonparametric Tests
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: November 3, 2018
    
    # check R
    if(!missing(R)){
      R <- as.integer(R)
      if(R < 1L) stop("Input 'R' must be a positive integer.")
    }
    
    # check delta
    if(!missing(delta)){
      delta <- as.numeric(delta)
      if(delta <= 0 | delta >= 1) stop("Input 'delta' must be a positive scalar between 0 and 1.")
    }
    
    # missing both?
    if(missing(R) & missing(delta)) stop("Either 'R' or 'delta' must be provided.")
    
    # check conf.level
    conf.level <- as.numeric(conf.level[1])
    if(conf.level <= 0 | conf.level >= 1) stop("Input 'conf.level' must be a positive scalar between 0 and 1.")
    crit.val <- qnorm(1 - (1 - conf.level)/2)
    
    # check sig.level
    sig.level <- as.numeric(sig.level)
    if(sig.level <= 0 | sig.level >= 1) stop("Input 'sig.level' must be a positive scalar between 0 and 1.")
    
    # check alternative
    alternative <- as.character(alternative[1])
    alt.opts <- c("two.sided", "one.sided")
    alternative <- pmatch(alternative, alt.opts)
    if(is.na(alternative)) stop("The 'alternative' must be either 'two.sided' or 'one.sided'")
    alternative <- alt.opts[alternative]
    if(alternative == "two.sided") sig.level <- sig.level / 2
    
    # calculate mcse
    if(missing(R)){
      mcse <- (1 / crit.val) * sig.level * delta
      R <- ceiling(sig.level * (1 - sig.level) / (mcse^2))
    } else {
      mcse <- sqrt(sig.level * (1 - sig.level) / R)
      delta <- sqrt((1 - sig.level) / sig.level) * crit.val / sqrt(R)
    }
    
    # return results
    if(alternative == "two.sided") sig.level <- sig.level * 2
    res <- list(mcse = mcse, R = R, delta = delta,
                conf.level = conf.level, 
                sig.level = sig.level,
                alternative = alternative)
    class(res) <- "mcse"
    return(res)
    
} # mcse.R