### CLEAR GLOBAL ENVIRONMENT ###
rm(list=ls(all.names = TRUE))

### REPRODUCIBILITY ###
set.seed(6)

### SOURCES ###
source('./Libraries.R') 
source('./Functions.R')
source('./Data Processing/PreprocessData.R')

### HERBERT'S DATA OR NOT
H <- T

##################
### DATA SETUP ###
##################
# get patterns from folder
raw_patterns <- if(H) getPatternsH(which_patterns = 1:4) else getPatterns(which_patterns = 1:4)
# preprocess patterns and save the scalings parameters
preprocess_patterns <- preprocessData(raw_patterns)
removed_dimensions <- preprocess_patterns$removed_dimensions
scalings <- preprocess_patterns$scalings
patterns <- preprocess_patterns$patterns
# set pattern variables
n_pattern <- length(patterns) # number of patterns
n_points <- sapply(patterns, FUN = function(x) dim(x)[1]) # length of each pattern

# plot (part of the) patterns
k_pattern <- 1
nrow_plots <- 3
par(mfrow=c(nrow_plots,nrow_plots))
for(i in 1:nrow_plots^2) plot(patterns[[k_pattern]][,i], type = 'l', ylim = c(-1,1), ylab = '')

#######################
### RESERVOIR SETUP ###
#######################
# dimensions
M <- ncol(patterns[[1]]) # dimension of output 
N <- 1000 # reservoir size
K <- 2000 # feature space dimension

# create weights
W_in_raw <- createInitialInputWeights() # dim=NxM (never modified after initialization)
b_raw <- createInitialBias() # dim=Nx1(never modified after initialization)
G_star_raw <- createInitialBackProjectionWeights() # dim=NxK 
G_star_raw <- diag(1,N,K)
F_raw <- createInitialForwardProjectionWeights() # dim=KxN 
F_raw <- diag(1,K,N)

leaking_rate <- 0.4 # leaking rate
input_scaling <- 0.15 # scaling of the input weights
bias_scaling <- 0.2 # scaling of the bias vector
spectral_radius_scaling <- 0.9 # scaling of the spectral radius of GF

Win <- W_in_raw * input_scaling
b <- b_raw * bias_scaling
spectral_radius_G_F <- spectralRadius(G_star_raw %*% F_raw)
F <- (F_raw / sqrt(spectral_radius_G_F)) * sqrt(spectral_radius_scaling)
G_star <- (G_star_raw / sqrt(spectral_radius_G_F)) * sqrt(spectral_radius_scaling)

learning_rate <- 0.5 # learning rate 
aperture <- c(5, 8, 5, 5) # aperture
reg_out <- 1 # regularizer for  W_out training
reg_D <- 0.01 # regularizer for D training

# how many first training points we discard
n_washout <- rep(50, n_pattern)
n_adapt <- rep(200, n_pattern)
n_learn <- n_points - n_adapt - n_washout

dr <- driveReservoirTrainingRFC(patterns = patterns, 
                                W_in = W_in, b = b, G_star = G_star, F=F, 
                                leaking_rate = leaking_rate,
                                reg_out = reg_out, reg_D = reg_D, 
                                aperture=aperture, learning_rate=learning_rate,
                                n_pattern = n_pattern, n_adapt = n_adapt, 
                                n_washout = n_washout, n_learn = n_learn)
# assign values from the training
all_training_r <- dr$all_training_r
original_outputs <- dr$original_outputs
training_outputs <- dr$training_outputs
W_out <- dr$W_out
D <- dr$D
G <- G_star + D
start_z <- dr$start_z
Cs <- dr$Cs

# plot some neuron states
par(mfrow=c(1,1))
plotNeuronStates(X=all_training_r, n=10, l=300)

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
drc <- driveReservoirWithRFC(Cs=Cs, 
                             G=G_star, b=b, F=F, W_out=W_out,
                             start_z=start_z,
                             n_pattern=n_pattern, n_learn = n_learn)
RFC_outputs <- drc$outputs

# NRMSE's
nrmses <- rep(NA, n_pattern)
for(p in 1:n_pattern){
  nrmses[p] <- nrmse(original_outputs[[p]], RFC_outputs[[p]])
  print(paste0('pattern ',p,' nrmse : ', round(nrmses[p],5)))
}

# plot the original patterns, the trained output, and the conceptor patterns
k_pattern <- 2
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
  op <- par(cex = 0.8) # fontsize in legend
  legend(1, 1, legend=c(paste0("NRMSE: ",round(nrmse(original_outputs[[k_pattern]][,i], RFC_outputs[[k_pattern]][,i]),5))))
}
