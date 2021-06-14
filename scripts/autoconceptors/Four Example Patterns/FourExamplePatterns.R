### CLEAR GLOBAL ENVIRONMENT ###
rm(list=ls(all.names = TRUE))

### REPRODUCIBILITY ###
set.seed(6)

### FOLDER PATH ###
folder_path1 <- 'Autoconceptors/'
folder_path2 <- 'Four Example Patterns/'
folder_path <- paste0(folder_path1, folder_path2)

### SOURCES ###
source('Source/Libraries.R')
source('Source/GeneralFunctions.R')
source('Source/PreprocessingFunctions.R')

##################
### DATA SETUP ###
##################
# read pre-generated patterns to avoid randomness
file_path <- paste0(folder_path, 'Results/RData/patterns.RData')
if(file.exists(file_path)){
  load(file_path)
}else{
  patterns <- createExamplePatterns(1:2000)
  save(patterns, file=file_path) 
}
# set pattern variables
n_pattern <- length(patterns) # number of patterns
n_points <- sapply(patterns, FUN = function(x) dim(x)[1]) # length of each pattern

# plot (part of) the patterns
par(mfrow=c(2,2), mar=c(5.1, 4.1, 4.1, 2.1))
n_plot <- 40
for(p in 1:4) plot(patterns[[p]][1:n_plot], type = 'l', ylim = c(-1,1), ylab = '')

##################
### PARAMETERS ###
##################
# scaling
W_in_scaling <- 1.5
b_scaling <- 0.3
W_star_scaling <- 1.4

# leaking rate
leaking_rate <- 1

# regularization parameters
reg_out <- 1e-3
reg_W <- 1e-6

# period lengths
n_washout <- rep(50,n_pattern)
n_load <- rep(50,n_pattern)
n_cue <- rep(15,n_pattern)
n_recall <- rep(100,n_pattern)

# apertures
apertures <- c(20,20,10,10)
apertures_2 <- apertures^-2

# adapting rates
C_adapt_rate_cue <- 0.02
C_adapt_rate_recall <- 0.01

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

###############
### LOADING ###
###############
# initialize empty collectors
all_training_r <- c()
all_training_r_old <- c()
all_training_W_target <- c()
all_training_output <- c()
# collect data from driving the reservoir with different drivers
for(p in 1:n_pattern){
  pattern <- patterns[[p]]
  print(paste('pattern',p))
  
  # learn period collectors
  r_collector <- matrix(0,N,n_load[p])
  r_old_collector <- matrix(0,N,n_load[p])
  W_target_collector <- matrix(0,N,n_load[p])
  output_collector <- matrix(0,M,n_load[p])
  
  # initial state
  r <- matrix(0,N,1)

  # create progress bar
  pb <- txtProgressBar(1,n_washout[p]+n_load[p])
  for(n in 1:(n_washout[p]+n_load[p])){
    u <- pattern[n,]
    r_old <- r
    W_target <- W_star %*% r_old + W_in %*% u
    r <- (1-leaking_rate) * r_old + leaking_rate * tanh(W_target + b)
    
    if(n > n_washout[p]){
      # collect states, output and W_target for computing W_out and W later
      r_collector[,n-n_washout[p]] <- r
      r_old_collector[,n-n_washout[p]] <- r_old
      W_target_collector[,n-n_washout[p]] <- W_target
      output_collector[,n-n_washout[p]] <- u
    }
    
    # update progress bar
    setTxtProgressBar(pb, n)
  }
  
  # construct concatenated matrices
  all_training_r <- cbind(all_training_r, r_collector)
  all_training_r_old <- cbind(all_training_r_old, r_old_collector)
  all_training_output <- cbind(all_training_output, output_collector)
  all_training_W_target <- cbind(all_training_W_target, W_target_collector)
}

# plot some of the reservoir states
plotNeuronStates(X=all_training_r, n=20, l=50)

# compute the output weights
cp_W_out <- computeOutputWeights(all_training_states = all_training_r, 
                                 all_training_output = all_training_output,
                                 reg_out = reg_out, verbose = T)
W_out <- cp_W_out$W_out

# recompute reservoir weights (loading the patterns into the reservoir)
cp_W <- recomputeReservoirWeights(all_training_states_old = all_training_r_old, 
                                  all_training_W_target = all_training_W_target,
                                  reg_W = reg_W, verbose = T)
W <- cp_W$W

##############
### RECALL ###
##############
# step 1 (Initial washout): starting from a zero network state, 
#                           the reservoir is driven with pj for n_washout steps, 
#                           in order to obtain a task-related reservoir state
# step 2 (Cueing): The reservoir is continued to be driven with pj 
#                  for another n_cue steps according to:
#                  r(n+1)=tanh(W*r(n)+Win pj(n)+b) and
#                  Cj(n+1)=Cj(n)+lamdba_cue*(r(n)r'(n)-Cj(n)r(n)r'(n)-a^-2*Cj(n))
#                  where Cj(0)=0.
#                  This leads to a matrix Cj_cue.
# step 3 (Autonomous recall): The network run is continued for another n_recall steps 
#                             in a mode where the input was switched off and replaced 
#                             by the newly computed reservoir weights W. 
#                             The matrix Cj_cue is continued to be adapted autonomously via:
#                             z(n+1)=Cj(n)tanh(W*z(n)+b) and
#                             Cj(n+1)=Cj(n)+lamdba_recall*(z(n)z'(n)-Cj(n)z(n)z'(n)-a^-2*Cj(n)).
#                             This leads to a matric Cj_recall.

# list that collects the C_cue and C_recall
C_cues <- list()
C_recalls <- list()
for(p in 1:n_pattern){
  print(paste0('pattern ',p))
  pattern <- patterns[[p]]
  
  # initial conceptor matrix
  C <- matrix(0,N,N)
  
  # initial washout
  r <- matrix(0,N,1)

  # create progress bar
  pb <- txtProgressBar(1,n_washout[p]+n_cue[p]+n_recall[p])
  for(n in 1:(n_washout[p]+n_cue[p]+n_recall[p])){
    if(n <= n_washout[p]){
      u <- pattern[n,]
      r <- (1 - leaking_rate)*r + leaking_rate * tanh(W_star %*% r + W_in %*% u + b)
    }else if(n > n_washout[p] && n <= n_washout[p]+ n_cue[p]){
      u <- pattern[n,]
      r <- (1 - leaking_rate)*r + leaking_rate * tanh(W_star %*% r + W_in %*% u + b)
      C <- C + C_adapt_rate_cue * ((r - C %*% r) %*% t(r) - apertures_2[p]*C)
      
      if(n == n_washout[p]+n_cue[p]) C_cues[[p]] <- C
    }else{
      # initial z
      if(n == n_washout[p]+n_cue[p]+1) z <- C %*% r
      u <- pattern[n,]
      r <- tanh(W %*% z + b)
      z <- (1 - leaking_rate) * z + leaking_rate * (C %*% r)
      C <- C + C_adapt_rate_recall * ((z - C %*% z) %*% t(z) - apertures_2[p]*C)
      
      if(n == n_washout[p] + n_cue[p] + n_recall[p]) C_recalls[[p]] <- C
    }
    
    # update progress bar
    setTxtProgressBar(pb,n)
  }
}

#######################################
### MEASURING QUALITY OF CONCEPTORS ###
#######################################
# period lengths
n_run <- n_points

sg_outputs <- list()
for(p in 1:n_pattern){
  print(paste0('pattern ',p))
  pattern <- patterns[[p]]
  
  # empty collectors
  output_collector <- matrix(0,n_run[p],M)
  
  # conceptor matrix
  C <- C_recalls[[p]]
  
  # initial washout
  z <- matrix(0,N,1)
  
  # create progress bar
  pb <- txtProgressBar(1,n_washout[p]+n_run[p])
  for(n in 1:(n_washout[p]+n_run[p])){
    r <- tanh(W %*% z + b)
    z <- C %*% r
    
    if(n>n_washout[p]) output_collector[n-n_washout[p],] <- W_out %*% r
    
    # update progress bar
    setTxtProgressBar(pb,n)
  }
  
  sg_outputs[[p]] <- output_collector
}

# shift the output of the sg_outputs for comparison
shifted_patterns <- lapply(1:n_pattern, function(p){
  shiftPattern(sg_outputs[[p]], 
               matrix(patterns[[p]][1:n_run[p],],n_run[p],1),
               max_shift = 10,
               shift_along = 'rows')
})

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
  legend('bottomleft', legend=c(paste0('error=',round(nrmses[p],3))))
  
}
