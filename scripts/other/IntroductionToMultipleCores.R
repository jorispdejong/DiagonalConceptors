# clear global environment
rm(list=ls(all.names = TRUE))

# libraries
suppressPackageStartupMessages(library(foreach)) # parallel 
suppressPackageStartupMessages(library(doParallel)) # parallel

###############
### EXAMPLE ###
###############
# a list of 30 1000x1000 matrices
matrixList <- lapply(1:30, FUN = function(x) matrix(runif(1e6),1e3,1e3))

### SINGLE CORE ###
system.time({
  lapply(matrixList, FUN = function(x) solve(x))
})

### MULTIPLE CORES ###
cl <- makeCluster(detectCores()[1]-1,outfile = NULL, verbose = TRUE)
registerDoParallel(cl)
system.time({
  foreach(i=1:length(matrixList)) %dopar% {solve(matrixList[[i]])}
})
stopCluster(cl)
