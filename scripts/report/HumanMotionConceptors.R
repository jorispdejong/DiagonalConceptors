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
file_path_preprocessed_patterns <- 'data/human motion/preprocessed/herbert/patterns.RData'
if(file.exists(file_path_preprocessed_patterns)){
  load(file_path_preprocessed_patterns)
}else{
  # get patterns from folder
  raw_patterns <- getPatterns(H = T)
  
  # preprocess patterns and save the scalings parameters
  preprocessed_patterns <- preprocessPatterns(raw_patterns, H = T)
  
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
k_patterns <- c(1,13,4,3)
names_p <- c('boxing 1', 'standing up', 'cart wheel', 'boxing 2')
n_plots <- 1
axis_show <- matrix(c(350,700,200,400,75,150,100,200),2,4)
par(mfrow=c(n_plots,length(k_patterns)), mar=c(3.1,2.6,2.5,1.1), oma=c(0,0,0,0))
for(i in 1:n_plots){
  for(p in 1:length(k_patterns)){
    plot(patterns[[k_patterns[p]]][,i], type = 'l', ylim = c(-1,1), 
         xlab = '', ylab = '', yaxt = 'n', lwd=2, cex.axis=1.3)
    #axis(side = 1, at = c(0,axis_show[1,p],axis_show[2,p]), labels = c(0,axis_show[1,p],axis_show[2,p]))
    axis(side = 2, at = c(-1,0,1), labels = c(-1,0,1), cex.axis=1.3)
    title(main = if(i==1) names_p[p], line = 0.7)
  }
}

##################
### PARAMETERS ###
##################
# scaling
W_in_scaling <- 0.2
b_scaling <- 0.8
W_star_scaling <- 1

# leaking rate
leaking_rate <- 0.3

# period lengths
n_washout <- rep(50,n_pattern)
n_learn <- n_points - n_washout

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
plotNeuronStates(X=all_training_states, n=20, l=n_learn[1]-20)

#############################
### (RE)COMPUTING WEIGHTS ###
#############################
# regularization
reg_out <- 0
reg_W <- 1e-3

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
apertures <- c(25,55,17,49,89,49,9,65,49,9,25,25,17,9,7)

# compute conceptors
Cs <- computeConceptors(R, apertures, verbose = verbose)

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
  z <- start_x[[p]]
  
  # create progress bar
  if(verbose) pb <- txtProgressBar(1,n_run[p], style = 3)
  for(n in 1:n_run[p]){
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
# compute the nrmses
nrmses <- nrmse(target_outputs,sg_outputs)

# plot (part of the) patterns
k_pattern <- 15
nrow_plots <- 3
which_patterns_plot <- 1:(nrow_plots^2)
par(mfrow=c(nrow_plots,nrow_plots), mar=rep(2,4))
for(i in 1:length(which_patterns_plot)){
  # plot the "worst" patterns
  which_patterns_plot <- order(nrmses[[k_pattern]]$all_nrmses, decreasing = T)[1:(nrow_plots^2)]
  plot(target_outputs[[k_pattern]][,which_patterns_plot[i]], type = 'l', 
       ylim = c(-1,1), ylab = '')
  lines(sg_outputs[[k_pattern]][,which_patterns_plot[i]], col=2, lwd=2)
  nrmse_i <- nrmse(target_outputs[[k_pattern]][,which_patterns_plot[i]], 
                   sg_outputs[[k_pattern]][,which_patterns_plot[i]])
  legend('bottomleft', legend=c(paste0('error=',round(nrmse_i,3)))) 
}

# compute errors
for(p in 1:n_pattern){
  if(p==1) if(verbose) print('Normalized Root Mean Square Error')
  if(verbose) print(paste0('p',p,
                           ': ', round(nrmses[[p]]$mean_nrmse,3), 
                           ', min=', round(nrmses[[p]]$min_nrmse,3),
                           ', max=',round(nrmses[[p]]$max_nrmse,3)))
}
nrmse_df <- t(data.frame('min'=round(sapply(nrmses, function(x) x$min_nrmse),3),
                         'max'=round(sapply(nrmses, function(x) x$max_nrmse),3),
                         'mean'=round(sapply(nrmses, function(x) x$mean_nrmse),3),
                         'std'=round(sqrt(sapply(nrmses, function(x) var(x$all_nrmses))),3)))
colnames(nrmse_df) <- 1:n_pattern


generate_plot_for_report <- F
if(generate_plot_for_report){
  ##############################
  ### SELF-GENERATE PATTERNS ###
  ##############################
  # let reservoir self-generate with diagonal conceptors
  n_run <- 4*n_learn
  sg_outputs <- list()
  for(p in 1:n_pattern){
    pattern <- patterns[[p]]
    if(verbose) print(paste('pattern',p))
    
    output_collector <- matrix(0,M,n_run[p])
    
    # initial state
    z <- start_x[[p]]
    
    # create progress bar
    if(verbose) pb <- txtProgressBar(1,n_run[p], style = 3)
    for(n in 1:n_run[p]){
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
  # plot (part of the) patterns
  k_patterns <- c(1,13,4,3)
  names_p <- c('boxing 1', 'standing up', 'cart wheel', 'boxing 2')
  n_plots <- 1
  par(mfrow=c(n_plots,length(k_patterns)), mar=c(3.1,2.6,2.5,1.1), oma=c(0,0,0,0))
  for(i in 1:n_plots){
    for(p in 1:length(k_patterns)){
      plot(sg_outputs[[k_patterns[p]]][,i], type = 'l', ylim = c(-1,1), 
           xlab = '', ylab = '', yaxt = 'n', col='#FF000088', lwd=3, cex.axis=1.3)
      lines(patterns[[k_patterns[p]]][n_washout[k_patterns[p]]:n_points[k_patterns[p]],i], lwd=2)
      axis(side = 2, at = c(-1,0,1), labels = c(-1,0,1), cex.axis=1.3)
      title(main = if(i==1) names_p[p], line = 1, cex.main=1.5)
    }
  }
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
                    'start_z'=start_x,
                    'n_learn'=n_learn,
                    'leaking_rate'=leaking_rate)
  save(reservoir, file = 'output/reservoirs/reservoir_c_hm_H.RData')
  if(verbose) print('Done!')
}

