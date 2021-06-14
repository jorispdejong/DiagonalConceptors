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
which_patterns <- 1:10
# read pre-generated patterns to avoid randomness
if(H){
  file_path_preprocessed_patterns <- paste0(folder_path, 'Results/RData/preprocessed_10patterns_H.RData')
}else{
  file_path_preprocessed_patterns <- paste0(folder_path, 'Results/RData/preprocessed_10patterns.RData')
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
W_in_scaling <- if(H) 0.1 else 0.1
b_scaling <- if(H) 0.2 else 0.2
W_star_scaling <- if(H) 1 else 1

# leaking rate
leaking_rate <- if(H) 0.15 else 0.12

# other parameters
iterative <- T
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
n_washout <- rep(50,n_pattern)
n_adapt <- ceiling(n_points * (1/3))
n_learn <- n_points - n_adapt - n_washout

# apertures
if(H){
  apertures <- rep(35,n_pattern)
  apertures[1] <- 10
  apertures[2] <- 10
  apertures[8] <- 8
}else{
  apertures <- rep(35,n_pattern)
  apertures[2] <- 9.5
  apertures[3] <- 10
  apertures[5] <- 10
  apertures[6] <- 8
}
apertures_2 <- apertures^-2

# learning rate vector
if(iterative) learning_rate <- rep(0.5,N)

### DRIVE RESERVOIR ###
# diagonal conceptors (initial randomly)
cs <- lapply(1:n_pattern, FUN=function(p) matrix(runif(N),N,1))

# initialize empty collectors
all_training_r <- c()
all_training_z_old <- c()
all_training_D_target <- c()
all_training_W_target <- c()
all_training_output <- c()
start_z <- list()
# collect data from driving the reservoir with different drivers
for(p in 1:n_pattern){
  pattern <- patterns[[p]]
  if(verbose) print(paste('pattern',p))
  
  # adapt period collectors
  if(!iterative) z_collector <- matrix(0,N,n_adapt[p])
  
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
      if(iterative){
        # online adaptation 
        cs[[p]] <- cs[[p]] + learning_rate * ((1 - cs[[p]]) * z^2 - apertures_2[p] * cs[[p]])
      }else{
        z_collector[,n-n_washout[p]] <- z
        
        # compute diagonal conceptors
        if(n==n_washout[p]+n_adapt[p]){
          R_adj <- rowMeans(z_collector^2)
          cs[[p]] <- R_adj*(R_adj+apertures_2[p])^-1
        } 
      }
      
      # save starting state
      if(n==n_washout[p]+n_adapt[p]) start_z[[p]] <- z
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
if(show_plots) plotNeuronStates(X=all_training_r, n=10, l=n_learn[1]-20)

#############################
### (RE)COMPUTING WEIGHTS ###
#############################
# regularization
reg_out <- if(H) 0.1 else 0.1
reg_W <- 0.1
compute_D <- T
if(compute_D) reg_D <- if(H) 0.1 else 0.01

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
sg_outputs <- list()
for(p in 1:n_pattern){
  if(verbose) print(paste0('Pattern ',p))
  # empty collectors
  p_collector <- matrix(0,M,n_learn[p])
  
  # initial states
  z <- start_z[[p]]
  
  # create progress bar
  if(verbose) pb <- txtProgressBar(1,n_learn[p], style = 3)
  for(n in 1:n_learn[p]){
    # update equation
    z_old <- z
    r <- tanh(W %*% z_old + b)
    z <- (1-leaking_rate) * z_old + leaking_rate * cs[[p]]*r
    
    # compute output
    p_collector[,n] <- W_out %*% r
    
    # update progress bar
    if(verbose) setTxtProgressBar(pb, n)
  }
  sg_outputs[[p]] <- t(p_collector)
}

####################
### PLOT RESULTS ###
####################
# plot (part of the) patterns
k_pattern <- 2
nrow_plots <- 4
which_patterns_plot <- 1:(nrow_plots^2)
# which_patterns_plot <- c(1:(nrow_plots^2-1),20)
par(mfrow=c(nrow_plots,nrow_plots), mar=rep(2,4))
for(i in 1:length(which_patterns_plot)){
  if(show_plots){
    plot(target_outputs[[k_pattern]][,which_patterns_plot[i]], type = 'l', ylim = c(-1,1), ylab = '')
    lines(sg_outputs[[k_pattern]][,which_patterns_plot[i]], col=2)
    nrmse_i <- nrmse(target_outputs[[k_pattern]][,which_patterns_plot[i]], 
                     sg_outputs[[k_pattern]][,which_patterns_plot[i]])
    legend('bottomleft', legend=c(paste0('error=',round(nrmse_i,3)))) 
  }
}

# compute errors
nrmses <- nrmse(target_outputs,sg_outputs)
for(p in 1:n_pattern){
  if(p==1) if(verbose) print('Normalized Root Mean Square Error')
  if(verbose) print(paste0('p',p,
                           ': ', round(nrmses[[p]]$mean_nrmse,3), 
                           ', min=', round(nrmses[[p]]$min_nrmse,3),
                           ', max=',round(nrmses[[p]]$max_nrmse,3)))
}
# 
# # show the worst patterns
# which_max_nrmses <- sapply(nrmses, function(x) which.max(x$all_nrmses))
# par(mfrow=c(3,4), mar=rep(2,4))
# for(p in 1:n_pattern){
#   if(show_plots){
#     plot(target_outputs[[k_pattern]][,which_max_nrmses[p]], type = 'l',
#          ylim = c(-1,1), ylab = '', main =  paste0('pattern ',p))
#     lines(sg_outputs[[k_pattern]][,which_max_nrmses[p]], col=2)
#     nrmse_i <- nrmse(target_outputs[[k_pattern]][,which_max_nrmses[p]], 
#                      sg_outputs[[k_pattern]][,which_max_nrmses[p]])
#     legend('bottomleft', legend=c(paste0('error=',round(nrmse_i,3)))) 
#   }
# }
