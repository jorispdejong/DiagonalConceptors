### CLEAR GLOBAL ENVIRONMENT ###
rm(list=ls(all.names = TRUE))

### REPRODUCIBILITY ###
set.seed(6)

### SOURCES ###
source('source/Libraries.R')
source('source/GeneralFunctions.R')
source('source/ConceptorsFunctions.R')
source('source/PreprocessingFunctions.R')

##################
### DATA SETUP ###
##################
# read pre-generated patterns to avoid randomness
file_path <- 'data/chaotic attractors/patterns.RData'
if(file.exists(file_path)){
  load(file_path)
}else{
  patterns <- getChaoticAttractorsPatterns(generate_new_data = T, L=1500)
  save(patterns, file=file_path) 
}
# set pattern variables
n_pattern <- length(patterns) # number of patterns
n_points <- sapply(patterns, FUN = function(x) dim(x)[1]) # length of each pattern

# plot the patterns
par(mfrow=c(2,2), mar=c(2.1, 3, 2.1, 2.1))
for(p in 1:4) plot(patterns[[p]][,1], patterns[[p]][,2], type = ifelse(p==4, 'p', 'l'), 
                   xlim = c(0,1), ylim = c(0,1), 
                   xlab = '', ylab = '')

##################
### PARAMETERS ###
##################
# scaling
W_in_scaling <- 1.4
b_scaling <- 1
W_star_scaling <- 1.5

# leaking rate
leaking_rate <- 1

# period lengths
n_washout <- rep(100, n_pattern)
n_learn <- n_points - n_washout

# other parameters
verbose <- T
save_results <- F

#######################
### RESERVOIR SETUP ###
#######################
# dimensions
M <- ncol(patterns[[1]]) # dimension of output 
N <- 500 # reservoir size

# create raw (not scaled yet) weights
W_in_raw <- createInitialInputWeights() # dim=NxM (never modified after initialization, only scaled)
b_raw <- createInitialBias() # dim=Nx1(never modified after initialization, only scaled)
W_star_raw <- createInitialReservoirWeights(spectral_radius=1) # spectral radius = 1

# scale weights
W_in <- W_in_scaling * W_in_raw
b <- b_scaling * b_raw
W_star <- W_star_scaling * W_star_raw

#######################
### DRIVE RESERVOIR ###
#######################
all_training_states <- c()
all_training_states_old <- c()
all_training_W_target <- c()
all_training_output <- c()
R <- list()
# collect data from driving the reservoir with different drivers
for(p in 1:n_pattern){
  pattern <- patterns[[p]]
  if(verbose) print(paste('pattern',p))
  
  # empty collectors
  states_collector <- matrix(0,N,n_learn[p])
  states_old_collector <- matrix(0,N,n_learn[p])
  W_target_collector <- matrix(0,N,n_learn[p])
  output_collector <- matrix(0,M,n_learn[p])
  
  # starting state
  x <- matrix(0,N,1)
  
  # create progress bar
  if(verbose) pb <- txtProgressBar(1,n_washout[p]+n_learn[p])
  for(n in 1:(n_washout[p]+n_learn[p])){
    u <- pattern[n,]
    x_old <- x
    W_target <- W_star %*% x_old + W_in %*% u
    x <- (1-leaking_rate) * x_old + leaking_rate * tanh(W_target+b)
    
    if(n > n_washout[p]){
      # collect states, output and W_target for computing W_out and W later
      states_collector[,n-n_washout[p]] <- x
      states_old_collector[,n-n_washout[p]] <- x_old
      W_target_collector[,n-n_washout[p]] <- W_target
      output_collector[,n-n_washout[p]] <- u
    }
    
    # update progress bar
    if(verbose) setTxtProgressBar(pb, n)
  }
  
  # correlation matrix
  R[[p]] <- (states_collector %*%t(states_collector))/n_learn[p]
  
  # construct concatenated matrices
  all_training_states <- cbind(all_training_states, states_collector)
  all_training_states_old <- cbind(all_training_states_old, states_old_collector)
  all_training_output <- cbind(all_training_output, output_collector)
  all_training_W_target <- cbind(all_training_W_target, W_target_collector)
}

# split training outputs matrix into a list
target_outputs <- splitMatrixToList(t(all_training_output), lengths = n_learn, splitAlong = 'rows')

# plot some neuron states
plotNeuronStates(X=all_training_states, n=20, l=20)

#############################
### (RE)COMPUTING WEIGHTS ###
#############################
# regularization constant
reg_out <- 0
reg_W <- 0.1

# compute the output weights
cp_W_out <- computeOutputWeights(all_training_states = all_training_states, 
                                 all_training_output = all_training_output,
                                 reg_out = reg_out, verbose = verbose)
W_out <- cp_W_out$W_out

# recompute reservoir weights (loading the patterns into the reservoir)
cp_W <- recomputeReservoirWeights(all_training_states_old = all_training_states_old, 
                                  all_training_W_target = all_training_W_target,
                                  reg_W = reg_W, verbose = verbose)
W <- cp_W$W

##########################
### COMPUTE CONCEPTORS ###
##########################
# apertures
apertures <- c(140,35,60,20)

# compute conceptors
Cs <- computeConceptors(R, apertures)
singular_values_R <- lapply(R, function(x) svd(x)$d)

##############################
### SELF-GENERATE PATTERNS ###
##############################
# let reservoir self-generate with diagonal conceptors
n_run <- n_learn
sg_outputs <- list()
for(p in 1:n_pattern){
  pattern <- patterns[[p]]
  if(verbose) print(paste('pattern',p))
  
  output_collector <- matrix(0,M,n_run[p])
  
  # initial state
  r <- matrix(runif(N),N,1)
  z <- Cs[[p]] %*% r
  
  # create progress bar
  if(verbose) pb <- txtProgressBar(1,n_run[p]+n_washout[p], style = 3)
  for(n in 1:(n_run[p]+n_washout[p])){
    r <- tanh(W %*% z + b)
    z <- (1-leaking_rate) * z + leaking_rate * (Cs[[p]] %*% r)
    
    if(n > n_washout[p]){
      output_collector[,n-n_washout[p]] <- W_out %*% r
    }
    
    # update progress bar
    if(verbose) setTxtProgressBar(pb, n)
  }
  
  sg_outputs[[p]] <- t(output_collector)
}

####################
### PLOT RESULTS ###
####################
# compare target outputs and self-generated outputs (plotting and Normalized Root Mean Square Error (NRMSE))
par(mfrow=c(2,4), mar=rep(1,n_pattern))
for(p in 1:n_pattern){
  # plot outputs
  plot(patterns[[p]][,1], patterns[[p]][,2], type = ifelse(p==4, 'p', 'l'), 
       xlim = c(0,1), ylim = c(0,1),
       xlab = '', ylab = '', xaxt='n', yaxt='n')
  plot(sg_outputs[[p]][,1], sg_outputs[[p]][,2], type = ifelse(p==4, 'p', 'l'),
       col = 2,
       xlim = c(0,1), ylim = c(0,1),
       xlab = '', ylab = '', xaxt='n', yaxt='n')
}

####################
### SAVE RESULTS ###
####################
if(save_results){
  if(verbose) print('Saving results...')
  reservoir <- list('patterns'=patterns,
                    'Cs'=Cs,
                    'W'=W,
                    'b'=b,
                    'W_out'=W_out,
                    'n_washout'=n_washout,
                    'leaking_rate'=leaking_rate)
  save(reservoir, file = 'output/reservoirs/reservoir_c_ca.RData')
  if(verbose) print('Done!')
}

