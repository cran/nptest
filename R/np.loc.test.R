np.loc.test <-
  function(x, y = NULL,
           alternative = c("two.sided", "less", "greater"),
           mu = 0, paired = FALSE, var.equal = FALSE, 
           median.test = FALSE, symmetric = TRUE,
           R = 9999, parallel = FALSE, cl = NULL,
           perm.dist = TRUE){
    # Nonparameteric Tests of Location Parameters
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: November 4, 2018
    
    
    ### check x
    x <- as.matrix(x)
    nvar <- ncol(x)
    
    ### check y
    if(!is.null(y)){
      y <- as.matrix(y)
      if(ncol(y) != nvar) stop("Inputs 'x' and 'y' must have the same number of columns.")
    }
    
    ### check type
    type <- ifelse(is.null(y), "One Sample",
                   ifelse(paired, "Paired", "Two Sample"))
    
    ### paired test?
    if(type == "Paired"){
      if(nrow(x) != nrow(y)) stop("Arguments 'x' and 'y' must have the same ", 
                                  ifelse(nvar == 1L, "length","number of rows"), 
                                  " when 'paired = TRUE'.")
      x <- x - y
    }
    
    ### permutation test
    if(type == "Two Sample"){
      pt <- rand.test.two(x = x, y = y, 
                          alternative = alternative, 
                          mu = mu, var.equal = var.equal, 
                          median.test = median.test, R = R,
                          parallel = parallel, cl = cl,
                          perm.dist = perm.dist)
    } else {
      pt <- rand.test.one(x = x,
                          alternative = alternative, 
                          mu = mu, symmetric = symmetric, 
                          median.test = median.test, R = R,
                          parallel = parallel, cl = cl,
                          perm.dist = perm.dist)
    } # end if(method == "Two Sample")
    
    ### get method
    if(median.test){
      symmetric <- ifelse(symmetric,
                          "Wilcoxon Signed Rank Test",
                          "Fisher Sign Test")
      if(type == "Paired") symmetric <- paste("Paired", symmetric)
      var.equal <- ifelse(var.equal, 
                          "Wilcoxon Rank Sum Test",
                          "Studentized Wilcoxon Test")
      method <- ifelse(any(type == c("One Sample", "Paired")),
                       symmetric,
                       var.equal)
    } else {
      var.equal <- ifelse(var.equal, 
                          "Two Sample t-test",
                          "Welch Two Sample t-test")
      method <- ifelse(type == "One Sample",
                       ifelse(symmetric,
                              "One Sample t-test",
                              "Johnson t-test"),
                       ifelse(type == "Paired",
                              paste("Paired", 
                                    ifelse(symmetric,
                                           "t-test",
                                           "Johnson t-test")),
                              var.equal))
    }
    
    ### return results
    class(pt) <- "np.loc.test"
    pt$method <- method
    return(pt)
    
  } # end np.loc.test.R
