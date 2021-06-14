### CLEAR GLOBAL ENVIRONMENT ###
rm(list=ls(all.names = TRUE))

### REPRODUCIBILITY ###
set.seed(6)

### SOURCES ###
source('Libraries.R') 
source('Functions.R')

##################
### DATA SETUP ###
##################
# example patterns setup
patterns <- createExamplePatterns(1:1000)
n_pattern <- length(patterns) # number of patterns
n_points <- sapply(patterns, FUN = function(x) dim(x)[1]) # length of each pattern

# plot the patterns
par(mfrow=c(2,2))
plot_length <- 30
for(i in 1:4) plot(patterns[[i]][1:plot_length], type = 'l', ylim = c(-1,1), ylab = '')

#######################
### RESERVOIR SETUP ###
#######################
# dimensions
M <- 1 # dimension of output
N <- 100 # reservoir size
K <- 300 # feature space dimension

# create weights
W_in_raw <- createInitialInputWeights() # dim=NxM (never modified after initialization)
b_raw <- createInitialBias() # dim=Nx1(never modified after initialization)
G_star_raw <- createInitialBackProjectionWeights() # dim=NxK 
F_raw <- createInitialForwardProjectionWeights() # dim=KxN 

leaking_rate <- 0.3 # leaking rate
input_scaling <- 1.2 # scaling of the input weights
bias_scaling <- 0.2 # scaling of the bias vector
spectral_radius_scaling <- 1.4 # scaling of the spectral radius of GF

Win <- W_in_raw * input_scaling
b <- b_raw * bias_scaling
spectral_radius_G_F <- spectralRadius(G_star_raw %*% F_raw)
F <- (F_raw / sqrt(spectral_radius_G_F)) * sqrt(spectral_radius_scaling)
G_star <- (G_star_raw / sqrt(spectral_radius_G_F)) * sqrt(spectral_radius_scaling)

learning_rate <- 0.5 # learning rate 
aperture <- 4.6 # aperture
reg_out <- 1e-4 # regularizer for  W_out training
reg_D <- 1e-4 # regularizer for  G training

# how many first training points we discard
n_washout <- rep(100, n_pattern)
n_adapt <- rep(500, n_pattern)
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
plotNeuronStates(X=all_training_r, n=10, l=plot_length)

# plot the original patterns and the trained output
par(mfrow=c(2,2))
plot_length <- 30
for(i in 1:4){
  plot(original_outputs[[i]][1:plot_length], type = 'l', ylim = c(-1,1), ylab = '')
  lines(training_outputs[[i]][1:plot_length], col=2)
}

############
### RFCs ###
############
drc <- driveReservoirWithRFC(Cs=Cs, 
                             G=G, b=b, F=F, W_out=W_out,
                             start_z=start_z,
                             n_pattern=n_pattern, n_learn=n_learn)
RFC_outputs <- drc$outputs

# NRMSE's
nrmses <- rep(NA, n_pattern)
for(p in 1:n_pattern){
  nrmses[p] <- nrmse(original_outputs[[p]], RFC_outputs[[p]])
  print(paste0('pattern ',p,' nrmse : ', round(nrmses[p],5)))
}

# plot the original patterns, the trained output, and the conceptor patterns
par(mfrow=c(2,2), mar=c(1,1,1,1))
for(p in 1:4){
  if(p==1 || p==2) plot_length <- 40
  if(p==3 || p==4) plot_length <- 20
  plot(original_outputs[[p]][1:plot_length], type = 'l', 
       ylim = c(-1,1),
       ylab = '', xlab = '', 
       xaxt='n', yaxt='n')
  #lines(training_outputs[[p]][1:plot_length], col=3)
  lines(RFC_outputs[[p]][1:plot_length], col=2, lty=2, lwd=4) 
  # Add a legend
  op <- par(cex = 0.9) # fontsize in legend
  legend(1, 1, legend=c(paste0("NRMSE: ",round(nrmses[p],5))))
}
