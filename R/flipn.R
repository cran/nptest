flipn <- function(n){
  # Generate all 2^n Sign-Flipping Vectors
  # Nathaniel E. Helwig (helwig@umn.edu)
  # last updated: October 17, 2018
  # adapted from "bincombinations" in e1071 R package
  
  ns <- 2^n
  sc <- matrix(0, n, ns)
  for(i in 1:n) sc[i,] <- rep(rep(c(-1, 1), each = ns/(2^i)), length = ns)
  sc
  
}