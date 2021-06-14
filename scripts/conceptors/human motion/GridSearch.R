### CLEAR GLOBAL ENVIRONMENT ###
rm(list=ls(all.names = TRUE))

### REPRODUCIBILITY ###
set.seed(6)

### SOURCES ###
source('source/Libraries.R')
source('source/GeneralFunctions.R')
source('source/ConceptorsFunctions.R')
source('source/PreprocessingFunctions.R')

### HERBERT'S DATA OR NOT ###
H <- F

##################
### DATA SETUP ###
##################
# read pre-generated patterns to avoid randomness
if(H){
  file_path_preprocessed_patterns <- 'data/stick person/preprocessed/herbert/patterns.RData'
}else{
  file_path_preprocessed_patterns <- 'data/stick person/preprocessed/joris/patterns.RData'
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
k_pattern <- 1
nrow_plots <- 3
par(mfrow=c(nrow_plots,nrow_plots))
for(i in 1:nrow_plots^2) plot(patterns[[k_pattern]][,i], type = 'l', ylim = c(-1,1), ylab = '')

##################
### PARAMETERS ###
##################
# leaking rate
if(H){
  leaking_rate <- 0.3
}else{
  leaking_rate <- 0.3
}

# period lengths
n_washout <- rep(70,n_pattern)
n_learn <- n_points - n_washout

# other parameters
show_plots <- T
verbose <- T
save_results <- F

#######################
### RESERVOIR SETUP ###
#######################
# dimensions
M <- ncol(patterns[[1]]) #dimension of output 
N <- 600 #reservoir size

# create raw (not scaled yet) weights
W_in_raw <- createInitialInputWeights() # dim=NxM (never modified after initialization, only scaled)
b_raw <- createInitialBias() # dim=Nx1(never modified after initialization, only scaled)
W_star_raw <- createInitialReservoirWeights(spectral_radius=1) # spectral radius = 1

# scaling
if(H){
  W_in_scaling <- 0.2
  b_scaling <- 0.8
  W_star_scaling <- 1
}else{
  W_in_scaling <- 0.2
  b_scaling <- 0.8
  W_star_scaling <- 1
}

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
start_x <- list()
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
      
      # collect starting state
      if(n == n_washout[p] + 1) start_x[[p]] <- x
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
if(show_plots) plotNeuronStates(X=all_training_states, n=20, l=n_learn[1]-20)


#############################
### (RE)COMPUTING WEIGHTS ###
#############################
# regularization
if(H){
  reg_out <- 1
  reg_W <- 1e-3
}else{
  reg_out <- 1
  reg_W <- 1e-3
}
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

############################
### GRID SEARCH APERTURE ###
############################
which_patterns <- 1:15 # which patterns will be tried out
apertures_seq <- seq(1,100,8) # apertures to loop over
nrmses <- matrix(0,length(apertures_seq),length(which_patterns)) # nrmse collector

if(verbose) pb <- txtProgressBar(1,length(apertures_seq), style = 3) # progress bar
if(verbose) print('Computing errors for different apertures...')
for(a in 1:length(apertures_seq)){
  apertures <- rep(apertures_seq[a], length(R))
  
  # compute conceptors
  Cs <- list()
  for(p in 1:length(R)){
    if(p %in% which_patterns){
      Cs[[p]] <- computeConceptor(R[[p]],apertures[p])
    }
    else{
      Cs[[p]] <- NA
    }
  }
  
  # let reservoir self-generate with conceptors
  sg_outputs <- list()
  for(p in 1:n_pattern){
    if(p %in% which_patterns){
      pattern <- patterns[[p]]

      # output collector
      output_collector <- matrix(0,M,n_learn[p])
      
      # initial state
      z <- start_x[[p]]
      
      for(n in 1:n_learn[p]){
        r <- tanh(W %*% z + b)
        z <- (1-leaking_rate) * z + leaking_rate * (Cs[[p]] %*% r)
        output_collector[,n] <- W_out %*% z
      }
      
      sg_outputs[[p]] <- t(output_collector) 
    }else{
      sg_outputs[[p]] <- NA
    }
  }
  
  # compute nrmse
  nrmses[a,] <- sapply(which_patterns, 
                       function(p) nrmse(target_outputs[[p]],sg_outputs[[p]])$mean_nrmse)
  
  # set progress bar
  setTxtProgressBar(pb,a)
}

# plot nrmses
par(mfrow=c(round(sqrt(length(which_patterns))),ceiling(sqrt(length(which_patterns)))))
for(p in 1:length(which_patterns)){
  plot(apertures_seq, nrmses[,p], type = 'l', 
       xlab = 'apertures', ylab = 'nrmse', 
       main = paste0('pattern: ', which_patterns[p]))
}

#######################
### FOUND APERTURES ###
#######################
best_apertures <- apertures_seq[apply(nrmses, 2, which.min)]

# compute the best conceptors and let reservoir self-generate
Cs <- computeConceptors(R,best_apertures, verbose = verbose)

# let reservoir self-generate with conceptors
sg_outputs <- list()
for(p in 1:n_pattern){
  pattern <- patterns[[p]]
  if(verbose) print(paste('pattern',p))
  
  output_collector <- matrix(0,M,n_learn[p])
  
  # initial state
  z <- start_x[[p]]
  
  # create progress bar
  if(verbose) pb <- txtProgressBar(1,n_learn[p], style = 3)
  for(n in 1:n_learn[p]){
    r <- tanh(W %*% z + b)
    z <- (1-leaking_rate) * z + leaking_rate * (Cs[[p]] %*% r)
    
    output_collector[,n] <- W_out %*% z 
    
    # update progress bar
    if(verbose) setTxtProgressBar(pb, n)
  }
  
  sg_outputs[[p]] <- t(output_collector)
}

####################
### PLOT RESULTS ###
####################
k_pattern <- 2
nrow_plots <- 3
par(mfrow=c(nrow_plots,nrow_plots))
for(i in 1:nrow_plots^2){
  plot(sg_outputs[[k_pattern]][,i], type = 'l', ylim = c(-1,1), ylab = '')
  lines(target_outputs[[k_pattern]][,i], col=2, lty=2, lwd=2)
}

# compare target outputs and self-generated outputs (Normalized Root Mean Square Error (NRMSE))
nrmses <- nrmse(target_outputs, sg_outputs)
for(p in 1:n_pattern){
  if(p==1) if(verbose) print('Normalized Root Mean Square Error')
  if(verbose) print(paste0('p',p,
                           ': ', round(nrmses[[p]]$mean_nrmse,3), 
                           ', min=', round(nrmses[[p]]$min_nrmse,3),
                           ', max=',round(nrmses[[p]]$max_nrmse,3)))
}

####################
### SAVE RESULTS ###
####################
if(save_results){
  reservoir <- list('preprocessed_patterns'=preprocessed_patterns,
                    'Cs'=Cs,
                    'W'=W,
                    'b'=b,
                    'W_out'=W_out,
                    'start_z'=start_x,
                    'leaking_rate'=leaking_rate)
  if(H){
    file_name_reservoir <- 'output/conceptors/stick person/herbert/reservoir.RData'
  }else{
    file_name_reservoir <- 'output/conceptors/stick person/joris/reservoir.RData'
  }
  if(verbose) print('saving reservoir...')
  save(reservoir, file = file_name_reservoir)
  if(verbose) print('done!')
}

