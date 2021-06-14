### CLEAR GLOBAL ENVIRONMENT ###
rm(list=ls(all.names = TRUE))

### REPRODUCIBILITY ###
set.seed(6)

### SOURCES ###
source('source/Libraries.R')
source('source/GeneralFunctions.R')
source('source/DiagonalConceptorsFunctions.R')
source('source/PreprocessingFunctions.R')

#################
### LOAD DATA ###
#################
load('output/reservoirs/reservoir_dc_ca.RData')

#####################
### SET VARIABLES ###
#####################
patterns <- reservoir$patterns
cs <- reservoir$cs
W <- reservoir$W
b <- reservoir$b
W_out <- reservoir$W_out
n_washout <- reservoir$n_washout
leaking_rate <- reservoir$leaking_rate

# system dimensions
N <- nrow(W)
M <- nrow(W_out)

# number of patterns
n_pattern <- length(patterns)
n_points <- sapply(patterns, function(x) nrow(x))

###################################
### LET RESERVOIR SELF GENERATE ###
###################################
n_run <- n_points
sg_outputs <- list()
for(p in 1:n_pattern){
  pattern <- patterns[[p]]
  print(paste('pattern',p))
  
  output_collector <- matrix(0,M,n_run[p])
  
  # initial state
  r <- matrix(runif(N),N,1)
  z <- cs[[p]] * r
  
  # create progress bar
  pb <- txtProgressBar(1,n_run[p]+n_washout[p], style = 3)
  for(n in 1:(n_run[p]+n_washout[p])){
    r <- tanh(W %*% z + b)
    z <- (1-leaking_rate) * z + leaking_rate * cs[[p]]*r
    
    if(n > n_washout[p]){
      output_collector[,n-n_washout[p]] <- W_out %*% r
    }
    
    # update progress bar
    setTxtProgressBar(pb, n)
  }
  
  sg_outputs[[p]] <- t(output_collector)
}

# plot the patterns
par(mfrow=c(2,4), mar=c(2.1, 3, 2.1, 2.1))
for(p in 1:4){
  plot(patterns[[p]][,1], patterns[[p]][,2], type = ifelse(p==4, 'p', 'l'), 
                   xlim = c(0,1), ylim = c(0,1), 
                   xlab = '', ylab = '')
  plot(sg_outputs[[p]][,1], sg_outputs[[p]][,2], type = ifelse(p==4, 'p', 'l'),
       col=2,
       xlim = c(0,1), ylim = c(0,1), 
       xlab = '', ylab = '')
}
