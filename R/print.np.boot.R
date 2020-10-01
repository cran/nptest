print.np.boot <-
  function(x, digits = 4, n = 6, ...){
    
    # get info
    t0 <- round(x$t0, digits)
    se <- round(x$se, digits)
    bias <- round(x$bias, digits)
    nstat <- length(t0)
    
    # print basic info
    if(nstat == 1L){
      cat("\nNonparametric Bootstrap of Univariate Statistic\n", sep = "")
      cat("using R = ", x$R, " bootstrap replicates\n", sep = "")
    } else {
      cat("\nNonparametric Bootstrap of Multivariate Statistic\n", sep = "")
      cat("using R = ", x$R, " bootstrap replicates\n", sep = "")
    }
    
    # se and bias
    if(nstat == 1L){
      cat("\n  t0:", t0)
      cat("\n  SE:", se)
      cat("\nBias:", bias, "\n") 
    } else if (nstat < n){
      mat <- rbind(t0, se, bias)
      rownames(mat) <- c("  t0:", "  SE:", "Bias:")
      cat("\n")
      print(mat)
    }
    
    # confint
    if(nstat == 1L){
      if(!is.null(x$bca)){
        if(!is.matrix(x$bca)){
          cat("\n", x$level * 100, "% BCa Confidence Interval:\n", sep = "")
        } else {
          cat("\nBCa Confidence Intervals:\n")
        }
        print(round(x$bca, digits = digits))
        cat("\n")
      } else if(!is.null(x$student)){
        if(!is.matrix(x$student)){
          cat("\n", x$level * 100, "% Studentized Confidence Interval:\n", sep = "")
        } else {
          cat("\nStudentized Confidence Intervals:\n")
        }
        print(round(x$student, digits = digits))
        cat("\n")
      } else if(!is.null(x$percent)){
        if(!is.matrix(x$percent)){
          cat("\n", x$level * 100, "% Percentile Confidence Interval:\n", sep = "")
        } else {
          cat("\nPercentile Confidence Intervals:\n")
        }
        print(round(x$percent, digits = digits))
        cat("\n")
      } else if(!is.null(x$basic)){
        if(!is.matrix(x$basic)){
          cat("\n", x$level * 100, "% Basic Confidence Interval:\n", sep = "")
        } else {
          cat("\nBasic Confidence Intervals:\n")
        }
        print(round(x$basic, digits = digits))
        cat("\n")
      } else if(!is.null(x$normal)){
        if(!is.matrix(x$normal)){
          cat("\n", x$level * 100, "% Normal Confidence Interval:\n", sep = "")
        } else {
          cat("\nNormal Confidence Intervals:\n")
        }
        print(round(x$normal, digits = digits))
        cat("\n")
      }
      
    } else if(nstat < n){
      if(!is.null(x$bca)){
        if(is.matrix(x$bca)){
          cat("\n", x$level * 100, "% BCa Confidence Intervals:\n", sep = "")
          print(round(x$bca, digits = digits))
        } else {
          cidim <- dim(x$bca)
          mid <- ceiling(median(1:cidim[1]))
          cat("\n", x$level[mid] * 100, "% BCa Confidence Intervals:\n", sep = "")
          print(round(x$bca[mid,,], digits = digits))
        }
        cat("\n")
      } else if(!is.null(x$student)){
        if(is.matrix(x$student)){
          cat("\n", x$level * 100, "% Studentized Confidence Intervals:\n", sep = "")
          print(round(x$student, digits = digits))
        } else {
          cidim <- dim(x$student)
          mid <- ceiling(median(1:cidim[1]))
          cat("\n", x$level[mid] * 100, "% Studentized Confidence Intervals:\n", sep = "")
          print(round(x$student[mid,,], digits = digits))
        }
        cat("\n")
      } else if(!is.null(x$percent)){
        if(is.matrix(x$percent)){
          cat("\n", x$level * 100, "% Percentile Confidence Intervals:\n", sep = "")
          print(round(x$percent, digits = digits))
        } else {
          cidim <- dim(x$percent)
          mid <- ceiling(median(1:cidim[1]))
          cat("\n", x$level[mid] * 100, "% Percentile Confidence Intervals:\n", sep = "")
          print(round(x$percent[mid,,], digits = digits))
        }
        cat("\n")
      } else if(!is.null(x$basic)){
        if(is.matrix(x$basic)){
          cat("\n", x$level * 100, "% Basic Confidence Intervals:\n", sep = "")
          print(round(x$basic, digits = digits))
        } else {
          cidim <- dim(x$basic)
          mid <- ceiling(median(1:cidim[1]))
          cat("\n", x$level[mid] * 100, "% Basic Confidence Intervals:\n", sep = "")
          print(round(x$basic[mid,,], digits = digits))
        }
        cat("\n")
      } else if(!is.null(x$normal)){
        if(is.matrix(x$normal)){
          cat("\n", x$level * 100, "% Normal Confidence Intervals:\n", sep = "")
          print(round(x$normal, digits = digits))
        } else {
          cidim <- dim(x$normal)
          mid <- ceiling(median(1:cidim[1]))
          cat("\n", x$level[mid] * 100, "% Normal Confidence Intervals:\n", sep = "")
          print(round(x$normal[mid,,], digits = digits))
        }
        cat("\n")
      }
    } # end if(nstat == 1L)
    
  } # end print.np.boot