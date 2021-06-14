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
load('output/reservoirs/reservoir_dc_pp.RData')

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

###################################
### LET RESERVOIR SELF GENERATE ###
###################################
# let reservoir self-generate with diagonal conceptors
n_run <- sapply(patterns, nrow)
sg_outputs <- list()
for(p in 1:n_pattern){
  pattern <- patterns[[p]]
  print(paste('pattern',p))
  
  output_collector <- matrix(0,M,n_run[p])
  
  # initial state
  r <- matrix(runif(N),N,1)
  z <- cs[[p]] * r
  
  # create progress bar
  pb <- txtProgressBar(1,n_run[p]+2*n_washout[p], style = 3)
  for(n in 1:(n_run[p]+2*n_washout[p])){
    r <- tanh(W %*% z + b) + rnorm(N,0,0.01)
    z <- (1-leaking_rate) * z + leaking_rate * cs[[p]]*r
    
    if(n > 2*n_washout[p]){
      output_collector[,n-2*n_washout[p]] <- W_out %*% r
    }
    
    # update progress bar
    setTxtProgressBar(pb, n)
  }
  
  sg_outputs[[p]] <- t(output_collector)
}

# shift the output of the sg_outputs for comparison
shifted_patterns <- lapply(1:n_pattern, function(p) shiftPattern(sg_outputs[[p]], 
                                                                 patterns[[p]],
                                                                 max_shift = 100,
                                                                 shift_along = 'rows'))

####################
### PLOT RESULTS ###
####################
# compare target outputs and self-generated outputs (plotting and Normalized Root Mean Square Error (NRMSE))
par(mfrow=c(2,2), mar=rep(3,4))
n_plot <- 50
nrmses <- rep(NA, n_pattern)
st <- 1
for(p in 1:n_pattern){
  nrmses[p] <- shifted_patterns[[p]]$nrmse
  print(paste0('pattern ',p,' nrmse : ', round(nrmses[p],5)))
  
  # plot outputs
  plot(shifted_patterns[[p]]$target_pattern[st:(st+n_plot)], type = 'l', 
       xlab = '', ylab = '', main = paste0('Pattern ',p),
       ylim = c(-1,1))
  lines(shifted_patterns[[p]]$pattern[st:(st+n_plot)], col=2, lwd=3, lty=2)
  legend('bottomleft', legend=c(paste0('error=',round(nrmses[p],5))))
}



