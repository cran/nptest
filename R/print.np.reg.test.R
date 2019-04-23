print.np.reg.test <-
  function(x, digits = 4, n = 10, ...){
    
    # get info
    homoscedastic <- ifelse(x$homosced, "Homoscedastic", "Heteroscedastic")
    method <- ifelse(x$method == "perm", "Permute",
                     ifelse(x$method == "flip", "Flip Sign", 
                            ifelse(x$method == "both", 
                                   "Flip Sign & Permute", NA)))
    univar <- is.null(x$univariate)
    if(is.na(method)){
      method.full <- c("Draper-Stoneman", "O'Gorman-Smith",
                       "Manly", "Freedman-Lane",
                       "ter Braak", "Still-White",
                       "Kennedy-Cade", "Huh-Jhun")
      method.shrt <- c("DS", "OS", "MA", "FL", 
                       "TB", "SW", "KC", "HJ")
      method <- method.full[match(x$method, method.shrt)] 
    }
    statistic <- paste(ifelse(x$homosced, "F =", "W ="), round(x$statistic, digits))
    p.value <- round(x$p.value, digits)
    
    # print info
    if(univar){
      cat("\nNonparametric Regression Test (", homoscedastic,")\n", sep = "")
    } else {
      cat("\nMultivariate Nonparametric Regression Test (", homoscedastic,")\n", sep = "")
    }
    cat("randomization method:", method, "\n")
    cat(statistic, ",  p-value = ", p.value, "\n", sep = "")
    if(length(x$coef) <= n){
      cat("sample estimates:", "\n", sep = "")
      print(round(x$coef, digits))
    }
    cat("\n")
    
  }