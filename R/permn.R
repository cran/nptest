permn <- function(n){
  # Generate all n! Permutation Vectors
  # Nathaniel E. Helwig (helwig@umn.edu)
  # last updated: February 25, 2018
  # adapted from "permutations" in e1071 R package
  
  if(n == 1){
    return(matrix(1))
  } else {
    P <- 1
    for(i in 2:n){
      x <- rbind(P, i)
      iseq <- 1:i
      indx <- c(iseq, iseq[-i])
      ncx <- ncol(x)
      P <- matrix(0, nrow = nrow(x), ncol = i * ncx)
      P[, 1:ncol(x)] <- x
      for(j in iseq[-i]){
        P[, j * ncx + 1:ncx] <- x[indx[1:i + j] ,]
      }
    }
    P
  }
}
