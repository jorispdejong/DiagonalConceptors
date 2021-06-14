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
par(mfrow=c(2,2), mar=rep(1.8,4), oma=c(0,0,1,0))
p_names <- c('Rössler', 'Lorenz', 'Mackey-Glass', 'Hénon')
for(p in 1:4) plot(patterns[[p]][,1], patterns[[p]][,2], type = ifelse(p==4, 'p', 'l'), 
                   xlim = c(0,1), ylim = c(0,1), xaxt='n', yaxt='n',
                   xlab = '', ylab = '', main = p_names[p])

##################
### PARAMETERS ###  
##################
# scaling
W_in_scaling <- 1.5
b_scaling <- 1
W_star_scaling <- 1.5

# leaking rate
leaking_rate <- 1

# apertures
apertures <- c(10,6,9,5)
apertures_2 <- apertures^-2

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
W_in <- W_in_raw * W_in_scaling
b <- b_raw * b_scaling
W_star <- W_star_raw * W_star_scaling

# period lengths
n_washout <- rep(100,n_pattern)
n_adapt <- rep(400,n_pattern)
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
  z_collector <- c()
  
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
    
    if(n > n_washout[p] && n <= n_washout[p] + n_adapt[p]){
      z_collector <- cbind(z_collector, z)
      
      R_adj <- rowMeans(z_collector^2)
      cs[[p]] <- R_adj*(R_adj+apertures_2[p])^-1
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
plotNeuronStates(X=all_training_r, n=20, l=50)

#############################
### (RE)COMPUTING WEIGHTS ###
#############################
# regularization
reg_out <- 0
reg_W <- 0.01
compute_D <- F
if(compute_D) reg_D <- 1e-2

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
    
    if(n > n_washout2[p]) output_collector[,n-n_washout2[p]] <- W_out %*% r
    
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
                    'cs'=cs,
                    'W'=W,
                    'b'=b,
                    'W_out'=W_out,
                    'n_washout'=n_washout2,
                    'leaking_rate'=leaking_rate)
  save(reservoir, file = 'output/reservoirs/reservoir_dc_ca.RData')
  if(verbose) print('Done!')
}

