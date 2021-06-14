### CLEAR GLOBAL ENVIRONMENT ###
rm(list=ls(all.names = TRUE))

### REPRODUCIBILITY ###
set.seed(6)

### SOURCES ###
source('source/Libraries.R')
source('source/GeneralFunctions.R')

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
leaking_rate <- reservoir$leaking_rate

# system dimensions
N <- nrow(W)
M <- nrow(W_out)

#########################################
### MORPH FROM PATTERN i TO PATTERN j ###
#########################################
pattern_i <- 1
pattern_j <- 2

# period lengths
n_washout <- 100
n_before_morph <- 1000
n_morph <- 50
n_after_morph <- 1000
n_run <- n_before_morph + n_morph + n_after_morph

# output collector
output_collector <- matrix(0,n_run,M)

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
  }
  
  # update progress bar
  setTxtProgressBar(pb, n)
}

# plot a few outputs to see how the figure changes over time
n_outputs <- 9
n_plot <- 250
max_plot <- nrow(output_collector)-n_plot
interval <- max_plot/(n_outputs-1)
begin_end_points <- t(sapply(1:n_outputs, function(i) round(c(1+(i-1)*interval,(i-1)*interval+n_plot))))

# plot the patterns
par(mfrow=c(ceiling(sqrt(n_outputs)),round(sqrt(n_outputs))), mar=rep(1,4))
for(i in 1:n_outputs){
  plot(output_collector[begin_end_points[i,1]:begin_end_points[i,2],1], 
       output_collector[begin_end_points[i,1]:begin_end_points[i,2],2], 
       type = 'l', 
       xlab = '', ylab = '', xaxt='n', yaxt='n')
}
