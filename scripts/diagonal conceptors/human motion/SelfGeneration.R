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
start_z <- reservoir$start_z
n_learn <- reservoir$n_learn
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
n_run <- n_learn
sg_outputs <- list()
for(p in 1:n_pattern){
  pattern <- patterns[[p]]
  print(paste('pattern',p))
  
  output_collector <- matrix(0,M,n_run[p])
  
  # initial state
  z <- start_z[[p]]
  
  # create progress bar
  pb <- txtProgressBar(1,n_run[p], style = 3)
  for(n in 1:(n_run[p])){
    r <- tanh(W %*% z + b)
    z <- (1-leaking_rate) * z + leaking_rate * cs[[p]]*r
    
    output_collector[,n] <- W_out %*% r
    
    # update progress bar
    setTxtProgressBar(pb, n)
  }
  
  sg_outputs[[p]] <- t(output_collector)
}

####################
### PLOT RESULTS ###
####################
# compute the nrmses
target_outputs <- lapply(1:n_pattern, function(p) patterns[[p]][(n_points[p]-n_learn[p]+1):n_points[p],])
nrmses <- nrmse(target_outputs,sg_outputs)

# plot (part of the) patterns
k_pattern <- 15
nrow_plots <- 4
which_patterns_plot <- 1:(nrow_plots^2)
which_patterns_plot <- order(nrmses[[k_pattern]]$all_nrmses, decreasing = T)[1:(nrow_plots^2)]
par(mfrow=c(nrow_plots,nrow_plots), mar=rep(2,4))
for(i in 1:length(which_patterns_plot)){
  # plot the "worst" patterns
  plot(target_outputs[[k_pattern]][,which_patterns_plot[i]], type = 'l', ylim = c(-1,1), ylab = '')
  lines(sg_outputs[[k_pattern]][,which_patterns_plot[i]], col=2)
  nrmse_i <- nrmse(target_outputs[[k_pattern]][,which_patterns_plot[i]], 
                   sg_outputs[[k_pattern]][,which_patterns_plot[i]])
  legend('bottomleft', legend=c(paste0('error=',round(nrmse_i,3)))) 
}

# compute errors
for(p in 1:n_pattern){
  if(p==1) print('Normalized Root Mean Square Error')
  print(paste0('p',p,
               ': ', round(nrmses[[p]]$mean_nrmse,3), 
               ', min=', round(nrmses[[p]]$min_nrmse,3),
               ', max=',round(nrmses[[p]]$max_nrmse,3)))
}
nrmse_df <- t(data.frame('min'=round(sapply(nrmses, function(x) x$min_nrmse),3),
                         'max'=round(sapply(nrmses, function(x) x$max_nrmse),3),
                         'mean'=round(sapply(nrmses, function(x) x$mean_nrmse),3),
                         'std'=round(sqrt(sapply(nrmses, function(x) var(x$all_nrmses))),3)))
colnames(nrmse_df) <- 1:n_pattern

