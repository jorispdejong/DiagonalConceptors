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
# read pre-generated patterns to avoid randomness
if(H){
  file_path_preprocessed_patterns <- paste0(folder_path, 'Results/RData/preprocessed_patterns_H.RData')
}else{
  file_path_preprocessed_patterns <- paste0(folder_path, 'Results/RData/preprocessed_patterns.RData')
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
k_pattern <- 6
nrow_plots <- 3
par(mfrow=c(nrow_plots,nrow_plots), mar=c(5.1, 4.1, 4.1, 2.1))
for(i in 1:nrow_plots^2) plot(patterns[[k_pattern]][,i], type = 'l', ylim = c(-1,1), ylab = '')

#######################
### RESERVOIR SETUP ###
#######################
# dimensions
M <- ncol(patterns[[1]]) # dimension of output 
N <- 1000 # reservoir size

# create weights
W_in_raw <- createInitialInputWeights() # dim=NxM (never modified after initialization)
b_raw <- createInitialBias() # dim=Nx1(never modified after initialization)

leaking_rate <- 0.4 # leaking rate
input_scaling <- 0.15 # scaling of the input weights
bias_scaling <- 0.2 # scaling of the bias vector

Win <- W_in_raw * input_scaling
b <- b_raw * bias_scaling

learning_rate <- 0.5 # learning rate 
aperture <- rep(4.5,n_pattern) # aperture
#aperture[15] <- 10
reg_out <- 3.5 # regularizer for  W_out training
reg_D <- 1 # regularizer for  D training
reg_G <- 3.5 # regularizer for  G training
compute_G <- T # recompute G or D

# how many first training points we discard
n_washout <- rep(50, n_pattern)
n_adapt <- rep(50, n_pattern)
n_learn <- n_points - n_adapt - n_washout

### DRIVE RESERVOIR ###
if(length(aperture)==1){
  apertures <- rep(aperture, n_pattern)
}else{
  if(length(aperture)!=n_pattern) stop('length of R and alpha must be equal')
  apertures <- aperture
}
apertures_2 <- apertures^-2

# empty collectors
all_training_z <- c()
all_training_cz_old <- c()
all_training_output <- c()
all_training_D_target <- c()
all_training_t <- c()
start_z <- list()
c_collectors <- list()
Cs <- list()
# collect data from driving the reservoir with different drivers
for(p in 1:n_pattern){
  print(paste0('Pattern ',p))
  pattern <- patterns[[p]]
  
  # empty collectors
  z_collector <- matrix(0,N,n_learn[p])
  cz_old_collector <- matrix(0,N,n_learn[p])
  p_collector <- matrix(0,M,n_learn[p])
  c_collector <- matrix(0,N,n_adapt[p])
  D_target_collector <- matrix(0,N,n_learn[p])
  t_collector <- matrix(0,N,n_learn[p])
  
  # initial states
  z <- matrix(0,N,1)
  c <- matrix(1,N,1)
  cz <- matrix(0,N,1)
  
  # create progress bar
  pb <- txtProgressBar(min = 1, max = n_washout[p] + n_adapt[p] + n_learn[p], style = 3)
  
  for(n in 1:(n_washout[p] + n_adapt[p] + n_learn[p])){
    u <- pattern[n,] 
    cz_old <- cz
    D_target <- Win %*% u
    t <- cz + D_target
    z <- tanh(t + b)
    cz <- c * z
    if(n > n_washout[p] && n <= n_adapt[p] + n_washout[p]){
      c <- c + learning_rate * ((cz - c*cz) * cz - apertures_2[p]*c)
      c_collector[,n - n_washout[p]] <- c
    }
    if(n == n_adapt[p] + n_washout[p]){
      Cs[[p]] <- c
      start_z[[p]] <- cz
    }
    if(n > n_washout[p] + n_adapt[p]){
      z_collector[,n - n_washout[p] - n_adapt[p]] <- z
      cz_old_collector[,n - n_washout[p] - n_adapt[p]] <- cz_old
      p_collector[,n - n_washout[p] - n_adapt[p]] <- u
      t_collector[,n - n_washout[p] - n_adapt[p]] <- t
      D_target_collector[,n - n_washout[p] - n_adapt[p]] <- D_target
    }
    # update progress bar
    setTxtProgressBar(pb, n)
  }
  
  c_collectors[[p]] <- c_collector
  all_training_z <- cbind(all_training_z, z_collector)
  all_training_cz_old <- cbind(all_training_cz_old, cz_old_collector)
  all_training_output <-cbind(all_training_output, p_collector)
  all_training_t <- cbind(all_training_t, t_collector)
  all_training_D_target <- cbind(all_training_D_target, D_target_collector)
}

### COMPUTE W OUT
print('computing W_out...')
W_out_training <- ridgeRegression(X=all_training_z, Y=all_training_output, reg = reg_out)
print('done!')
W_out <- W_out_training[[1]]
training_outputs_matrix <- W_out %*% all_training_z
print(paste0('mean NRMSE W_out: ', round(W_out_training[[2]],5)))

# split training outputs matrix into a list
training_outputs <- splitMatrixToList(t(training_outputs_matrix), lengths = n_learn, splitAlong = 'rows')
original_outputs <- splitMatrixToList(t(all_training_output), lengths = n_learn, splitAlong = 'rows')

### RECOMPUTE D or G
if(compute_G){
  print('recomputing G...')
  G_training <- ridgeRegression(X=all_training_cz_old, Y=all_training_t, reg = reg_G)
  print('done!')
  G <- G_training[[1]]
  print(paste0('mean NRMSE G: ', round(G_training[[2]],5)))
  print(paste0('mean abs G: ', round(mean(colMeans(G)),5)))
}else{
  print('computing D...')
  D_training <- ridgeRegression(X=all_training_cz_old, Y=all_training_D_target, reg = reg_D)
  print('done!')
  D <- D_training[[1]]
  print(paste0('mean NRMSE D: ', round(D_training[[2]],5)))
  print(paste0('mean abs D: ', round(mean(colMeans(D)),5)))
  G <- diag(1,N,N) + D 
}

# plot some neuron states
par(mfrow=c(1,1))
plotNeuronStates(X=all_training_z, n=10, l=300)

# plot the original patterns and the trained output
k_pattern <- 1
nrow_plots <- 2
par(mfrow=c(nrow_plots,nrow_plots))
for(i in 1:nrow_plots^2){
  plot(original_outputs[[k_pattern]][,i], type = 'l', ylim = c(-1,1), ylab = '')
  lines(training_outputs[[k_pattern]][,i], col=2)
}

############
### RFCs ###
############
# empty collectors
RFC_outputs <- list()
for(p in 1:n_pattern){
  print(paste0('Pattern ',p))
  
  # empty collectors
  p_collector <- matrix(0,M,n_learn[p])
  
  # initial states
  c <- Cs[[p]]
  cz <- start_z[[p]]
  # create progress bar
  pb <- txtProgressBar(min = 1, max = n_learn[p], style = 3)
  for(n in 1:n_learn[p]){
    z <- tanh(G%*%cz + b)
    cz <- c * z
    p_collector[,n] <- W_out %*% z
    
    # update progress bar
    setTxtProgressBar(pb, n)
  }
  RFC_outputs[[p]] <- t(p_collector)
}

# NRMSE's
nrmses <- rep(NA, n_pattern)
for(p in 1:n_pattern){
  nrmses[p] <- nrmse(original_outputs[[p]], RFC_outputs[[p]])
  print(paste0('pattern ',p,' nrmse : ', round(nrmses[p],5)))
}

# plot the original patterns, the trained output, and the conceptor patterns
k_pattern <- 1
nrow_plots <- 3
par(mfrow=c(nrow_plots,nrow_plots), mar=c(1,1,1,1))
for(i in 1:nrow_plots^2){
  plot(original_outputs[[k_pattern]][,i], 
       type = 'l', 
       ylim = c(-1,1),
       ylab = '', xlab = '', 
       xaxt='n', yaxt='n')
  #lines(training_outputs[[k_pattern]][,i], col=2)
  lines(RFC_outputs[[k_pattern]][,i], col=2, lty=2, lwd=4)
  # Add a legend
  #op <- par(cex = 0.8) # font size in legend
  #legend(1, 1, legend=c(paste0("NRMSE: ",round(nrmse(original_outputs[[k_pattern]][,i], RFC_outputs[[k_pattern]][,i]),5))))
}

