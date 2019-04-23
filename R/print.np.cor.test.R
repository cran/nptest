print.np.cor.test <-
  function(x, digits = 4, ...){
    
    # get info
    independent <- ifelse(x$independent, "Independence", "Correlation")
    H1 <- ifelse(x$alternative == "two.sided", 
                 "not equal to ",
                 ifelse(x$alternative == "greater", 
                        "greater than ", "less than "))
    H1 <- paste0("alternative hypothesis: true correlation is ",
                 H1, x$null.value)
    statistic <- paste("t =", round(x$statistic, digits))
    p.value <- round(x$p.value, digits)
    estimate <- x$estimate
    names(estimate) <- "cor"
    
    # print info
    cat("\nNonparametric ", independent ," Test\n", sep = "")
    cat(H1, "\n")
    cat(statistic, ",  p-value = ", p.value, "\n", sep = "")
    cat("sample estimate:", "\n", sep = "")
    print(estimate)
    cat("\n")
    
    
  }