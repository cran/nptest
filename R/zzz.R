# Some code (for startup message) is re-purposed from the R package
# "mclust" (Scrucca et al., 2016) https://cran.r-project.org/package=mclust

nptestStartupMessage <- 
  function(){
    msg <- c(paste0("             _            _   
 _ __  _ __ | |_ ___  ___| |_ 
| '_ \\| '_ \\| __/ _ \\/ __| __|
| | | | |_) | ||  __/\\__ \\ |_ 
|_| |_| .__/ \\__\\___||___/\\__|
      |_|       version ", 
                    packageVersion("nptest"), "\n"),
             "\nType 'citation(\"nptest\")' to cite this package.\n")
    return(msg)
  }

.onAttach <- 
  function(lib, pkg){
    msg <- nptestStartupMessage()
    if(!interactive()) msg[1] <- paste("Package 'nptest' version", packageVersion("nptest"))
    packageStartupMessage(msg)
    invisible()
  }