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
load('output/reservoirs/reservoir_dc_hm_H.RData')

#####################
### SET VARIABLES ###
#####################
patterns <- reservoir$patterns
cs <- reservoir$cs
W <- reservoir$W
b <- reservoir$b
W_out <- reservoir$W_out
n_learn <- reservoir$n_learn
start_z <- reservoir$start_z
leaking_rate <- reservoir$leaking_rate

# system dimensions
N <- nrow(W)
M <- nrow(W_out)

# number of patterns
n_pattern <- length(patterns)
n_points <- sapply(patterns, nrow)

###################################
### LET RESERVOIR SELF GENERATE ###
###################################
multiplier <- 4
n_run <- multiplier * n_learn
sg_outputs <- list()
for(p in 1:n_pattern){
  pattern <- patterns[[p]]
  print(paste('pattern',p))
  
  # output collector
  output_collector <- matrix(0,M,n_run[p])
  
  # initial state
  z <- start_z[[p]]
  
  # re starting points
  restarting_points <- seq(1, multiplier*n_learn[p], n_learn[p])
  
  # create progress bar
  pb <- txtProgressBar(1,n_run[p], style = 3)
  for(n in 1:n_run[p]){
    if(n %in% restarting_points) z <- start_z[[p]]
    
    # update state
    r <- tanh(W %*% z + b)
    z <- (1-leaking_rate) * z + leaking_rate * cs[[p]]*r
    
    # collect output
    output_collector[,n] <- W_out %*% r
    
    # update progress bar
    setTxtProgressBar(pb, n)
  }
  
  sg_outputs[[p]] <- t(output_collector)
}

####################
### PLOT RESULTS ###
####################
# plot (part of the) patterns
k_patterns <- c(1,13,4,3)
names_p <- c('boxing 1', 'standing up', 'cart wheel', 'boxing 2')
n_plots <- 1
par(mfrow=c(n_plots,length(k_patterns)), mar=c(3.1,2.6,2.5,1.1), oma=c(0,0,0,0))
for(i in 1:n_plots){
  for(p in 1:length(k_patterns)){
    plot(sg_outputs[[k_patterns[p]]][,i], type = 'l', ylim = if(p==1) c(-1,1),
         xlab = '', ylab = '', yaxt = 'n', col='#FF000088', lwd=3, 
         cex.axis=1.3)
    lines(tail(patterns[[k_patterns[p]]][,i],n_learn[k_patterns[p]]), lwd=2)
    axis(side = 2, at = c(-1,0,1), labels = c(-1,0,1), cex.axis=1.3)
    title(main = if(i==1) names_p[p], line = 1, cex.main=1.5)
  }
}
