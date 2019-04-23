print.mcse <- 
  function(x, digits = 4, ...){
    
    cat("\nMonte Carlo Standard Errors for Nonparametric Tests\n", sep = "")
    cat("alternative hypothesis:", x$alternative, "\n")
    cat("Monte Carlo Std Err:", round(x$mcse, digits), "\n")
    cat("number of resamples:", x$R, "\n")
    cat("accuracy of approx:", round(x$delta, digits), "\n")
    cat("  confidence level:", round(x$conf.level, digits), "\n")
    cat("significance level:", round(x$sig.level, digits), "\n\n")
    
  } 