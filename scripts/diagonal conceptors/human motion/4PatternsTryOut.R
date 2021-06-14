### CLEAR GLOBAL ENVIRONMENT ###
rm(list=ls(all.names = TRUE))

### REPRODUCIBILITY ###
set.seed(6)

### FOLDER PATH ###
folder_path1 <- 'Diagonal Conceptors/'
folder_path2 <- 'Stick Person/'
folder_path <- paste0(folder_path1, folder_path2)

### SOURCES ###
source('Source/Libraries.R')
source('Source/GeneralFunctions.R')
source('Source/DiagonalConceptorsFunctions.R')
source('Source/PreprocessingFunctions.R')

### HERBERT'S DATA OR NOT ###
H <- T

##################
### DATA SETUP ###
##################
which_patterns <- 1:4
# read pre-generated patterns to avoid randomness
if(H){
  file_path_preprocessed_patterns <- paste0(folder_path, 'Results/RData/preprocessed_4patterns_H.RData')
}else{
  file_path_preprocessed_patterns <- paste0(folder_path, 'Results/RData/preprocessed_4patterns.RData')
} 
if(file.exists(file_path_preprocessed_patterns)){
  load(file_path_preprocessed_patterns)
}else{
  # get patterns from folder
  raw_patterns <- getPatterns(H = H, which_patterns = which_patterns)
  
  # preprocess patterns and save the scalings parameters
  preprocessed_patterns <- preprocessPatterns(raw_patterns, H = H)
  
  # save preprocessed patterns
  save(preprocessed_patterns, file=file_path_preprocessed_patterns) 
}
# get patterns
removed_dimensions <- preprocessed_patterns$removed_dimensions
scalings <- preprocessed_patterns$scalings
patterns <- preprocessed_patterns$patterns

# set pattern variables
n_pattern <- length(patterns) # number of patterns
n_points <- sapply(patterns, FUN = function(x) dim(x)[1]) # length of each pattern

# plot (part of the) patterns
k_pattern <- 2
nrow_plots <- 3
par(mfrow=c(nrow_plots,nrow_plots), mar=rep(3,4))
for(i in 1:nrow_plots^2) plot(patterns[[k_pattern]][,i], type = 'l', ylim = c(-1,1), ylab = '')

##################
### PARAMETERS ###
##################
# scaling
W_in_scaling <- 0.2
b_scaling <- 0.2
W_star_scaling <- 1

# leaking rate
leaking_rate <- 0.1

# other parameters
show_plots <- T
verbose <- T
save_results <- F

#######################
### RESERVOIR SETUP ###
#######################
# dimensions
M <- ncol(patterns[[1]]) # dimension of output 
N <- 1000 # reservoir size

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
n_learn <- n_points - n_washout

# diagonal conceptors (initial randomly)
cs_initial <- lapply(1:n_pattern, FUN=function(p) matrix(runif(N),N,1))

# initialize empty collectors
all_training_r <- c()
all_training_z_old <- c()
all_training_D_target <- c()
all_training_W_target <- c()
all_training_output <- c()
R_adj <- list()
start_z <- list()
# collect data from driving the reservoir with different drivers
for(p in 1:n_pattern){
  pattern <- patterns[[p]]
  if(verbose) print(paste('pattern',p))
  
  # learn period collectors
  r_collector <- matrix(0,N,n_learn[p])
  z_collector <- matrix(0,N,n_learn[p])
  z_old_collector <- matrix(0,N,n_learn[p])
  D_target_collector <- matrix(0,N,n_learn[p])
  W_target_collector <- matrix(0,N,n_learn[p])
  output_collector <- matrix(0,M,n_learn[p])
  
  # initial state
  r <- matrix(0,N,1)
  z <- cs_initial[[p]] * r
  
  # create progress bar
  if(verbose) pb <- txtProgressBar(1,n_washout[p]+n_learn[p])
  for(n in 1:(n_washout[p]+n_learn[p])){
    u <- pattern[n,]
    z_old <- z
    D_target <- W_in %*% u
    W_target <- W_star %*% z_old + D_target
    r <- tanh(W_target + b)
    z <- (1-leaking_rate) * z_old + leaking_rate * cs_initial[[p]]*r
    
    if(n==n_washout[p]) start_z[[p]] <- z
    
    if(n > n_washout[p]){
      # collect states 
      z_collector[,n-n_washout[p]] <- z
      r_collector[,n-n_washout[p]] <- r
      z_old_collector[,n-n_washout[p]] <- z_old
      D_target_collector[,n-n_washout[p]] <- D_target
      W_target_collector[,n-n_washout[p]] <- W_target
      output_collector[,n-n_washout[p]] <- u
    }
    
    # update progress bar
    if(verbose) setTxtProgressBar(pb, n)
  }
  
  R_adj[[p]] <- rowMeans(z_collector^2)
  
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
if(show_plots) plotNeuronStates(X=all_training_r, n=10, l=n_learn[1]-20)

#############################
### (RE)COMPUTING WEIGHTS ###
#############################
# regularization
reg_out <- 1e-4
reg_W <- 1e-4
compute_D <- T
if(compute_D) reg_D <- 1e-4

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

###################################
### COMPUTE DIAGONAL CONCEPTORS ###
###################################
# apertures
apertures <- rep(35,n_pattern)
apertures_2 <- apertures^-2

# diagonal conceptors
cs <- lapply(1:n_pattern, function(p) R_adj[[p]]*(R_adj[[p]]+apertures_2[p])^-1)

#################################
### SELF-GENERATING RESERVOIR ###
#################################
n_run <- n_learn
sg_outputs <- list()
for(p in 1:n_pattern){
  pattern <- patterns[[p]]
  if(verbose) print(paste('pattern',p))
  
  # output collector
  output_collector <- matrix(0,M,n_run[p])
  
  # initial state
  z <- start_z[[p]]
  
  # create progress bar
  if(verbose) pb <- txtProgressBar(1,n_washout[p]+n_run[p])
  for(n in 1:(n_washout[p]+n_run[p])){
    u <- pattern[n,]
    r <- tanh(W %*% z + b)
    z <- (1-leaking_rate) * z + leaking_rate * cs[[p]] * cs_initial[[p]] * r
    
    if(n > n_washout[p]){
      # collect output
      output_collector[,n-n_washout[p]] <- W_out %*% r
    }
    
    # update progress bar
    if(verbose) setTxtProgressBar(pb, n)
  }
  
  sg_outputs[[p]] <- output_collector
}

####################
### PLOT RESULTS ###
####################
# compare target outputs and self-generated outputs (plotting and Normalized Root Mean Square Error (NRMSE))
par(mfrow=c(2,2), mar=c(1.5,3.5,1,1))
n_plot <- 20
for(p in 1:n_pattern){
  # plot outputs
  if(show_plots) plot(target_outputs[1:n_plot], type = 'l', 
                      xlab = '', ylab = '', xaxt = 'n',
                      ylim = c(-1,1))
  if(show_plots) lines(shifted_patterns[[p]]$pattern[1:n_plot], col=2, lwd=3, lty=2)
  if(show_plots) legend('bottomleft', legend=c(paste0('error = ', 
                                                      round(shifted_patterns[[p]]$nrmse,5))))
  
}

nrmses <- sapply(shifted_patterns, function(x) x$nrmse)
for(p in 1:n_pattern){
  if(verbose){
    if(p==1) print('Normalized Root Mean Square Error')
    print(paste0('p',p,
                 ': ', round(nrmses[p],5)))
  }
}

####################
### SAVE RESULTS ###
####################
if(save_results){
  reservoir <- list('cs'=cs,
                    'W'=W,
                    'b'=b,
                    'W_out'=W_out,
                    'start_cz'=start_x,
                    'leaking_rate'=leaking_rate)
  save(reservoir, file = paste0(folder_path,'Results/RData/reservoir.RData'))
}

