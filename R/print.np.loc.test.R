print.np.loc.test <-
  function(x, digits = 4, n = 6, ...){
    
    # get info
    statistic <- paste("t =", round(x$statistic, digits))
    p.value <- round(x$p.value, digits)
    twosamp <- c("Two Sample t-test",
                 "Welch Two Sample t-test",
                 "Wilcoxon Rank Sum Test",
                 "Studentized Wilcoxon Test")
    paired <- substr(x$method, start = 1, stop = 6) == "Paired"
    univar <- is.null(x$univariate)
    
    parameter <- ifelse(any(x$method == twosamp) | paired,
                        ifelse(x$median.test, "median of differences", "difference of means"),
                        ifelse(x$median.test, "median", "mean"))
    if(univar | (max(x$null.value) - min(x$null.value)) == 0){
      H1 <- ifelse(x$alternative == "two.sided", 
                   "not equal to ",
                   ifelse(x$alternative == "greater", 
                          "greater than ", "less than "))
      H1 <- paste0("alternative hypothesis: true ", parameter," is ",
                   H1, x$null.value[1])
    } else {
      H1 <- paste0("alternative hypothesis: ", x$alternative, ";   parameter: ", parameter)
    }
    
    estimate <- x$estimate
    if(univar){
      type <- ifelse(any(x$method == twosamp),
                     ifelse(x$median.test, 
                            "median of the differences", 
                            "difference of the means"),
                     ifelse(paired,
                            ifelse(x$median.test, 
                                   ifelse(x$symmetric, "(pseudo)median of the differences", 
                                          "median of the differences"), 
                                   "mean of the differences"),
                            ifelse(x$median.test, 
                                   ifelse(x$symmetric, "(pseudo)median", "median"), 
                                   "mean")))
      names(estimate) <- type
    } else {
      names(estimate) <- paste0("mu", 1:length(estimate))
    }
    
    
    # print info
    if(univar){
      cat("\nNonparametric Location Test (", x$method, ")\n", sep = "")
    } else {
      cat("\nMultivariate Nonparametric Location Test (", x$method, ")\n", sep = "")
    }
    cat(H1, "\n")
    cat(statistic, ",  p-value = ", p.value, "\n", sep = "")
    if(length(estimate) < n){
      cat(ifelse(univar, "sample estimate:", 
                 "sample estimates:"), "\n", sep = "")
      print(estimate)
    }
    cat("\n")
    
  }