### CLEAR GLOBAL ENVIRONMENT ###
rm(list=ls(all.names = TRUE))

### REPRODUCIBILITY ###
set.seed(6)

### SOURCES ###
source('source/Libraries.R')
source('source/GeneralFunctions.R')

# load reservoir
load('output/reservoirs/reservoir_dc_hm_H.RData')

# set variables
patterns <- reservoir$patterns
cs <- reservoir$cs
W <- reservoir$W
b <- reservoir$b
W_out <- reservoir$W_out
start_z <- reservoir$start_z
n_learn <- reservoir$n_learn
leaking_rate <- reservoir$leaking_rate

n_pattern <- length(patterns) # number of patterns
n_points <- sapply(patterns, FUN = function(x) nrow(x)) # length of each pattern

# system dimensions
N <- nrow(W)
M <- nrow(W_out)

# patterns
pattern_i <- 1
pattern_j <- 4

# mu
mu_min <- 0
mu_max <- 1

# period lengths
n_before_morph <- n_learn[pattern_i]
n_morph <- 150
n_after_morph <- n_learn[pattern_j]
n_run <- n_before_morph + n_morph + n_after_morph

# output collector
output_collector <- matrix(0,n_run,M)
mu_collector <- rep(NA,n_run)

# initial state
z <- start_z[[pattern_i]]

# create progress bar
pb <- txtProgressBar(1,n_run, style = 3)
for(n in 1:(n_run)){
  # set mu
  if(n < n_before_morph){
    mu <- mu_min
  }else if(n >= n_before_morph && n <= n_before_morph + n_morph){
    mu <- mu_min + ((n - n_before_morph)/n_morph) * (mu_max - mu_min)
    # if(n == n_before_morph) print(paste0('begin mu=',mu))
    # if(n == n_before_morph + n_morph) print(paste0('end mu=',mu))
  }else{
    mu <- mu_max
  }
  
  # nudge reservoir to correct state after morphing
  if(n >= n_before_morph && n <= n_before_morph + n_morph){
    z <- (1-mu)*z + mu*start_z[[pattern_j]]
  }
  
  # update equation
  r <- tanh(W %*% z + b)
  c_morph <- (1-mu)*cs[[pattern_i]] + mu*cs[[pattern_j]]
  z <- (1-leaking_rate) * z + leaking_rate * c_morph * r
  
  # collect output
  output_collector[n,] <- W_out %*% r
  mu_collector[n] <- mu
  
  # update progress bar
  setTxtProgressBar(pb, n)
}

# plot a few outputs to see how the figure changes over time
which_dims <- c(3,4,5)
main_title <- ''

# align 
target <- patterns[[pattern_j]][(n_points[pattern_j]-n_learn[pattern_j]+1):n_points[pattern_j],]
observed <- output_collector[n_learn[pattern_i]:nrow(output_collector),]
pa <- phaseAlignPatterns(observed, target)

####################
### PLOT RESULTS ###
####################
colfunc<-colorRampPalette(c("royalblue","springgreen"))
par(mfrow=c(length(which_dims),3), mar=c(1,1,3,1), oma=c(0,0,3,0))
for(i in 1:length(which_dims)){
  # plot the target pattern i in the left column
  plot(patterns[[pattern_i]][,which_dims[i]], type = 'l', lwd=2, 
       xaxt='n', yaxt = 'n', col=colfunc(3)[1],
       main = if(i==1) 'Boxing', cex.main=1.7)
  
  # plot morphed pattern in the middle column
  plot(output_collector[,which_dims[i]], type = 'l', 
       xlab = '', ylab = '', xaxt='n', yaxt = 'n', 
       col=colfunc(3)[2], lwd=3,
       main = if(i==1) 'Morph between boxing and cart wheel', cex.main=1.7)
  lines(tail(patterns[[pattern_i]][,which_dims[i]],n_learn[pattern_i]), 
        lwd=3, lty=2, col=colfunc(3)[1])
  lines((n_learn[pattern_i]+pa$shift+1):(n_learn[pattern_i]+pa$shift+n_learn[pattern_j]),
        tail(patterns[[pattern_j]][,which_dims[i]],n_learn[pattern_j]), 
        lwd=3, lty=2, col=colfunc(3)[3])
  
  # plot the target pattern j in the right column
  plot(patterns[[pattern_j]][,which_dims[i]], type = 'l', lwd=3, 
       col=colfunc(3)[3], 
       xaxt='n', yaxt = 'n',
       main = if(i==1) 'Cart wheel', cex.main=1.7)
}
title(main = 'Diagonal conceptors', outer = T, cex.main=3)
