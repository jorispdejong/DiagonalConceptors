### CLEAR GLOBAL ENVIRONMENT ###
rm(list=ls(all.names = TRUE))

### REPRODUCIBILITY ###
set.seed(6)

### SOURCES ###
source('source/Libraries.R')
source('source/GeneralFunctions.R')
source('source/DiagonalConceptorsFunctions.R')
source('source/PreprocessingFunctions.R')

##################
### DATA SETUP ###
##################
# read pre-generated patterns to avoid randomness
file_path <- 'data/periodic patterns/diagonal conceptors/patterns5000.RData'
if(file.exists(file_path)){
  load(file_path)
}else{
  patterns <- createExamplePatterns(generate_new_data=T, L=1000, 
                                    folder_path='data/periodic patterns/diagonal conceptors/')
  save(patterns, file=file_path) 
}
# set pattern variables
n_pattern <- length(patterns) # number of patterns
n_points <- sapply(patterns, FUN = function(x) dim(x)[1]) # length of each pattern

# plot (part of) the patterns
par(mfrow=c(2,2))
n_plot <- 40
for(p in 1:4) plot(patterns[[p]][1:n_plot], type = 'l', ylim = c(-1,1), ylab = '')

##################
### PARAMETERS ###
##################
# scaling
W_in_scaling <- 2
b_scaling <- 2
W_star_scaling <- 1

# leaking rate
leaking_rate <- 1

# apertures
apertures <- rep(8,n_pattern)
apertures_2 <- apertures^-2

# regularization
reg_out <- 0
reg_W <- 1e-3
compute_D <- F
if(compute_D) reg_D <- 1e-3

# other parameters
verbose <- T

#######################
### RESERVOIR SETUP ###
#######################
# dimensions
M <- ncol(patterns[[1]]) # dimension of output 
N <- 100 # reservoir size

# create raw (not scaled yet) weights
W_in_raw <- createInitialInputWeights() # dim=NxM (never modified after initialization, only scaled)
b_raw <- createInitialBias() # dim=Nx1(never modified after initialization, only scaled)
W_star_raw <- createInitialReservoirWeights(spectral_radius=1) # spectral radius = 1

# scale weights
W_in <- W_in_raw * W_in_scaling
b <- b_raw * b_scaling
W_star <- W_star_raw * W_star_scaling

# period lengths
n_washout <- rep(100,n_pattern)
n_adapt <- rep(500,n_pattern)
n_learn <- n_points - n_adapt - n_washout

### DRIVE RESERVOIR ###
# diagonal conceptors (initial randomly)
cs <- lapply(1:n_pattern, FUN=function(p) matrix(runif(N),N,1))

# initialize empty collectors
all_training_r <- c()
all_training_z_old <- c()
all_training_D_target <- c()
all_training_W_target <- c()
all_training_output <- c()
# collect data from driving the reservoir with different drivers
for(p in 1:n_pattern){
  pattern <- patterns[[p]]
  if(verbose) print(paste('pattern',p))
  
  # adapt period collectors
  z_collector <- matrix(0,N,n_adapt[p])
  
  # learn period collectors
  r_collector <- matrix(0,N,n_learn[p])
  z_old_collector <- matrix(0,N,n_learn[p])
  D_target_collector <- matrix(0,N,n_learn[p])
  W_target_collector <- matrix(0,N,n_learn[p])
  output_collector <- matrix(0,M,n_learn[p])
  
  # initial state
  r <- matrix(0,N,1)
  z <- cs[[p]] * r
  
  # create progress bar
  if(verbose) pb <- txtProgressBar(1,n_washout[p]+n_adapt[p]+n_learn[p])
  for(n in 1:(n_washout[p]+n_adapt[p]+n_learn[p])){
    u <- pattern[n,]
    z_old <- z
    D_target <- W_in %*% u
    W_target <- W_star %*% z_old + D_target
    r <- tanh(W_target + b)
    z <- (1-leaking_rate) * z_old + leaking_rate * cs[[p]]*r
    
    if(n > n_washout[p] && n <= n_washout[p]+n_adapt[p]){
      z_collector[,n-n_washout[p]] <- z
      
      # compute diagonal conceptors
      if(n==n_washout[p]+n_adapt[p]){
        R_adj <- rowMeans(z_collector^2)
        cs[[p]] <- R_adj*(R_adj+apertures_2[p])^-1
      } 
    }
    
    if(n > n_washout[p]+n_adapt[p]){
      # collect states, output and G_target for computing W_out and G later
      r_collector[,n-n_washout[p]-n_adapt[p]] <- r
      z_old_collector[,n-n_washout[p]-n_adapt[p]] <- z_old
      D_target_collector[,n-n_washout[p]-n_adapt[p]] <- D_target
      W_target_collector[,n-n_washout[p]-n_adapt[p]] <- W_target
      output_collector[,n-n_washout[p]-n_adapt[p]] <- u
    }
    
    # update progress bar
    if(verbose) setTxtProgressBar(pb, n)
  }
  
  # construct concatenated matrices
  all_training_r <- cbind(all_training_r, r_collector)
  all_training_z_old <- cbind(all_training_z_old, z_old_collector)
  all_training_output <- cbind(all_training_output, output_collector)
  all_training_W_target <- cbind(all_training_W_target, W_target_collector)
  all_training_D_target <- cbind(all_training_D_target, D_target_collector)
}

# split training outputs matrix into a list
target_outputs <- splitMatrixToList(t(all_training_output), lengths = n_learn, splitAlong = 'rows')

# plot some neuron states
plotNeuronStates(X=all_training_r, n=4, l=50)

#############################
### (RE)COMPUTING WEIGHTS ###
#############################
# compute the output weights
cp_W_out <- computeOutputWeights(all_training_states = all_training_r, 
                                 all_training_output = all_training_output,
                                 reg_out = reg_out, verbose = verbose)
W_out <- cp_W_out$W_out

# recompute reservoir weights (loading the patterns into the reservoir)
if(compute_D){
  cp_D <- recomputeReservoirWeights(all_training_states_old = all_training_z_old,
                                    all_training_W_target = all_training_D_target,
                                    reg_W = reg_D, verbose = verbose)
  D <- cp_D$W
  W <- W_star + D
}else{
  cp_W <- recomputeReservoirWeights(all_training_states_old = all_training_z_old, 
                                    all_training_W_target = all_training_W_target,
                                    reg_W = reg_W, verbose = verbose)
  W <- cp_W$W
}

##############################
### SELF-GENERATE PATTERNS ###
##############################
# let reservoir self-generate with diagonal conceptors
n_run <- n_points
n_washout2 <- 2*n_washout
sg_outputs <- list()
for(p in 1:n_pattern){
  pattern <- patterns[[p]]
  if(verbose) print(paste('pattern',p))
  
  output_collector <- matrix(0,M,n_run[p])
  
  # initial state
  r <- matrix(runif(N),N,1)
  z <- cs[[p]] * r
  
  # create progress bar
  if(verbose) pb <- txtProgressBar(1,n_run[p]+n_washout2[p], style = 3)
  for(n in 1:(n_run[p]+n_washout2[p])){
    r <- tanh(W %*% z + b)
    z <- (1-leaking_rate) * z + leaking_rate * cs[[p]]*r
    
    if(n > n_washout2[p]){
      output_collector[,n-n_washout2[p]] <- W_out %*% r
    }
    
    # update progress bar
    if(verbose) setTxtProgressBar(pb, n)
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
par(mfrow=c(n_pattern,4), mar=c(1,3,2.5,1), oma=c(2,0,2,1))
n_plot <- 30
nrmses <- rep(NA, n_pattern)
st <- c(8,8,3,3)
tr_st <- splitMatrixToList(all_training_r, n_learn)
R <- lapply(tr_st, function(X) X%*%t(X)/nrow(X))
singular_values_R <- lapply(R, function(x) svd(x)$d)
for(p in 1:n_pattern){
  # observed and target pattern
  nrmses[p] <- shifted_patterns[[p]]$nrmse
  if(verbose) print(paste0('pattern ',p,' nrmse : ', round(nrmses[p],10)))
  
  # plot outputs
  plot(shifted_patterns[[p]]$target_pattern[st[p]:(st[p]+n_plot)], type = 'l', 
       xlab = '', ylab = '', xaxt='n', yaxt='n',
       main = if(p==1) 'driver and y', cex.main=1.5,
       ylim = c(-1,1))
  axis(side = 1, at = c(0,15,30), labels = c('0','15','30'), cex.axis = 1.5)
  axis(side = 2, at = c(-1,0,1), labels = c('-1','0','1'), cex.axis = 1.5)
  lines(shifted_patterns[[p]]$pattern[st[p]:(st[p]+n_plot)], col=2, lwd=3, lty=2)
  legend('bottomleft', legend=c(round(nrmses[p],5)), x.intersp = 0, cex=1.3)
  
  # random reservoir states
  n_random_states <- 4
  random_neurons <- sample(1:N, n_random_states)
  matplot(t(tr_st[[p]][random_neurons,1:n_plot]), type='l', 
          xlab='', ylab = '', xaxt='n', yaxt='n',
          main = if(p==1) 'reservoir states', cex.main=1.5,
          ylim = c(-1,1))
  axis(side = 1, at = c(0,15,30), labels = c('0','15','30'), cex.axis = 1.5)
  axis(side = 2, at = c(-1,0,1), labels = c('-1','0','1'), cex.axis = 1.5)
  
  # log10 PC energy
  plot(log10(singular_values_R[[p]]), type = 'l', lwd=2,
       xlab='', ylab='', xaxt='n', yaxt='n', ylim=c(-20,10),
       main = if(p==1) 'log10 PC energy', cex.main=1.5)
  axis(side = 1, at = c(0,50,100), labels = c('0','50','100'), cex.axis = 1.5)
  axis(side = 2, at = c(-20,-10,0,10), labels = c('-20','-10','0','10'), cex.axis = 1.5)
  
  # leading PC energy
  plot(singular_values_R[[p]][1:10], type = 'l',lwd=2,
       xlab='', ylab='', xaxt='n', yaxt='n', ylim = c(0,1200),
       main = if(p==1) 'leading PC energy', cex.main=1.5)
  axis(side = 1, at = c(0,5,10), labels = c('0','5','10'), cex.axis = 1.5)
  axis(side = 2, at = c(0,600,1200), labels = c('0','600','1200'), cex.axis = 1.5)
}

####################
### SAVE RESULTS ###
####################
if(save_results){
  if(verbose) print('Saving results...')
  reservoir <- list('patterns'=patterns, 
                    'cs'=cs,
                    'W'=W,
                    'b'=b,
                    'W_out'=W_out,
                    'n_washout'=n_washout2,
                    'leaking_rate'=leaking_rate)
  save(reservoir, file = 'output/reservoirs/reservoir_dc_pp.RData')
  if(verbose) print('Done!')
}
