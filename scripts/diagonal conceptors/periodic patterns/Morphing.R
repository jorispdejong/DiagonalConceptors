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
leaking_rate <- reservoir$leaking_rate

# system dimensions
N <- nrow(W)
M <- nrow(W_out)

#########################################
### MORPH FROM PATTERN i TO PATTERN j ###
#########################################
pattern_i <- 3
pattern_j <- 1

# period lengths
n_washout <- 100
n_before_morph <- 50
n_morph <- 30
n_after_morph <- 50
n_run <- n_before_morph + n_morph + n_after_morph

# collectors
output_collector <- matrix(0,n_run,M)
mu_collector <- rep(NA,n_run)

# initial state
r <- matrix(runif(N),N,1)
z <- cs[[pattern_i]] * r

# create progress bar
pb <- txtProgressBar(1,n_run+n_washout, style = 3)
for(n in 1:(n_run+n_washout)){
  # set mu
  if(n < n_washout + n_before_morph){
    mu <- 0
  }else if(n >= n_washout + n_before_morph && n <= n_washout + n_before_morph + n_morph){
    mu <- (n - n_washout - n_before_morph)/n_morph
    if(n == n_washout + n_before_morph) print(paste0('begin mu=',mu))
    if(n == n_washout + n_before_morph + n_morph) print(paste0('end mu=',mu))
  }else{
    mu <- 1
  }
  
  # update equation
  r <- tanh(W %*% z + b)
  z <- (1-leaking_rate) * z + leaking_rate * ((1-mu)*cs[[pattern_i]] + mu*cs[[pattern_j]]) * r
  
  # collect output
  if(n > n_washout){
    output_collector[n-n_washout,] <- W_out %*% r
    mu_collector[n-n_washout] <- mu
  }
  
  # update progress bar
  setTxtProgressBar(pb, n)
}

# plot the patterns
par(mfrow=c(2,1), mar=c(5.1, 4.1, 4.1, 2.1))
plot(output_collector, type = 'l', ylim = c(-1,1), xlab = '', ylab = '')
plot(mu_collector, col=2, lty=2, lwd=2)
