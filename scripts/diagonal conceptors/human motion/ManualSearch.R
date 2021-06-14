### CLEAR GLOBAL ENVIRONMENT ###
rm(list=ls(all.names = TRUE))

### REPRODUCIBILITY ###
set.seed(6)

### SOURCES ###
source('source/Libraries.R')
source('source/GeneralFunctions.R')
source('source/DiagonalConceptorsFunctions.R')
source('source/PreprocessingFunctions.R')

### HERBERT'S DATA OR NOT ###
H <- T

##################
### DATA SETUP ###
##################
# read pre-generated patterns to avoid randomness
if(H){
  file_path_preprocessed_patterns <- 'data/human motion/preprocessed/herbert/patterns.RData'
}else{
  file_path_preprocessed_patterns <- 'data/human motion/preprocessed/joris/patterns.RData'
} 
if(file.exists(file_path_preprocessed_patterns)){
  load(file_path_preprocessed_patterns)
}else{
  # get patterns from folder
  raw_patterns <- getPatterns(H = H)
  
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
k_patterns <- c(2,8,9,12)
names_p <- c('boxing', 'sitting down', 'jogging', 'standing up')
n_plots <- 4
axis_show <- matrix(c(350,700,200,400,75,150,100,200),2,4)
par(mfrow=c(n_plots,length(k_patterns)), mar=c(3.1,2.6,2.5,1.1), oma=c(0,0,2.1,0))
for(i in 1:n_plots){
  for(p in 1:length(k_patterns)){
    plot(patterns[[k_patterns[p]]][,i], type = 'l', ylim = c(-1,1), 
         xlab = '', ylab = '', xaxt='n', yaxt = 'n')
    axis(side = 1, at = c(0,axis_show[1,p],axis_show[2,p]), labels = c(0,axis_show[1,p],axis_show[2,p]))
    axis(side = 2, at = c(-1,0,1), labels = c(-1,0,1))
    title(main = if(i==1) names_p[p], line = 0.7)
  }
  title('First dimension of a few patterns', outer = T, cex.main = 2)
}

##################
### PARAMETERS ###
##################
# scaling
W_in_scaling <- if(H) 0.2 else 0.1
b_scaling <- if(H) 0.8 else 0.2
W_star_scaling <- if(H) 1 else 1

# leaking rate
leaking_rate <- if(H) 0.3 else 0.12

# other parameters
verbose <- T
save_results <- F

#######################
### RESERVOIR SETUP ###
#######################
# dimensions
M <- ncol(patterns[[1]]) # dimension of output 
N <- 600 # reservoir size

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
  # apertures[1] <- 8
  # apertures[2] <- 8
  # apertures[4] <- 8
  # apertures[8] <- 8
  # apertures[15] <- 8
}else{
  apertures <- rep(100,n_pattern)
  apertures[15] <- 5
}
apertures_2 <- apertures^-2

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
plotNeuronStates(X=all_training_r, n=10, l=n_learn[1]-20)

#############################
### (RE)COMPUTING WEIGHTS ###
#############################
# regularization
reg_out <- if(H) 0.01 else 0.1
reg_W <- 0.1
compute_D <- F
if(compute_D) reg_D <- if(H) 0.1 else 1

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
k_pattern <- 15
nrow_plots <- 4
which_patterns_plot <- 1:(nrow_plots^2)
# which_patterns_plot <- c(1:(nrow_plots^2-1),20)
par(mfrow=c(nrow_plots,nrow_plots), mar=rep(2,4))
for(i in 1:length(which_patterns_plot)){
  plot(target_outputs[[k_pattern]][,which_patterns_plot[i]], type = 'l', ylim = c(-1,1), ylab = '')
  lines(sg_outputs[[k_pattern]][,which_patterns_plot[i]], col=2)
  nrmse_i <- nrmse(target_outputs[[k_pattern]][,which_patterns_plot[i]], 
                   sg_outputs[[k_pattern]][,which_patterns_plot[i]])
  legend('bottomleft', legend=c(paste0('error=',round(nrmse_i,3)))) 
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


