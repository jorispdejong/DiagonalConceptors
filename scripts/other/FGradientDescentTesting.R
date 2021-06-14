### CLEAR GLOBAL ENVIRONMENT ###
rm(list=ls(all.names = TRUE))

### REPRODUCIBILITY ###
set.seed(123)

### SOURCES ###
source('./Libraries.R') 
source('./Functions.R')

##################
### DATA SETUP ###
##################
pattern <- matrix(sin(2*pi*(1:1000)/8.8342522), ncol=1)
# set pattern variables
n_points <- length(pattern) # length of each pattern

# plot (part of the) patterns
par(mar = c(5.1, 4.1, 4.1, 2.1)) # Set the margin to default
par(mfrow=c(1,1))
plot_length <- 40
plot(pattern[1:plot_length], type = 'l', ylim = c(-1,1), ylab = '')

#######################
### RESERVOIR SETUP ###
#######################
# dimensions
M <- 1 #dimension of output 
N <- 20 #reservoir size

# create weights
W_in_raw <- createInitialInputWeights() # dim=NxM (never modified after initialization)
b_raw <- createInitialBias() # dim=Nx1(never modified after initialization)
W_star_raw <- createInitialReservoirWeights() # dim=NxN

leaking_rate <- 1
input_scaling <- 0.8
bias_scaling <- 0.1
spectral_radius_scaling <- 1

W_in <- W_in_raw * input_scaling
b <- b_raw * bias_scaling
W_star <- W_star_raw * spectral_radius_scaling

# how many first training points we discard
n_washout <- 100
n_learn <- n_points - n_washout

### LEARNING ###
start_X <- matrix(0,N,1)
X_collector <- matrix(0,N,n_learn)
X_old_collector <- matrix(0,N,n_learn)
W_target_collector <- matrix(0,N,n_learn)
output_collector <- matrix(0,M,n_learn)
x <- matrix(0,N,1)
for(n in 1:(n_learn+n_washout)){
  u <- pattern[n,]
  x_old <- x
  W_target <- W_star %*% x_old + W_in %*% u
  x <- (1-leaking_rate)*x_old + leaking_rate*tanh(W_target+b)
  
  if(n==n_washout) start_X[,1] <- x
  
  if(n > n_washout){
    X_collector[,n-n_washout] <- x
    X_old_collector[,n-n_washout] <- x_old
    W_target_collector[,n-n_washout] <- W_target
    output_collector[,n-n_washout] <- u
  }
}

# correlation matrix
R <- (X_collector %*%t(X_collector))/n_learn


# compute output matrix and pattern readout
W_out_training <- ridgeRegression(X=X_collector, Y=output_collector, reg = 1e-8)
W_out <- W_out_training[[1]]
training_outputs_matrix <- W_out %*% X_collector
print(paste0('mean NRMSE W_out: ', round(W_out_training[[2]],5)))
#print(paste0('mean abs W_out: ', round(mean(colMeans(W_out)),5)))

# compute W
W_training <- ridgeRegression(X=X_old_collector,Y=W_target_collector,reg=1e-8)
W <- W_training[[1]]
print(paste0('mean NRMSE W: ', round(W_training[[2]],5)))
print(paste0('Spectral radius W: ', round(spectralRadius(W),5)))

# assign values from the training
all_training_states <- X_collector
original_outputs <- output_collector
training_outputs <- training_outputs_matrix

# plot some neuron states
par(mfrow=c(1,1))
plotNeuronStates(X=all_training_states, n=10, l=plot_length)

# plot the original patterns and the trained output
plot(original_outputs[1:plot_length], type = 'l', ylim = c(-1,1), ylab = '')
lines(training_outputs[1:plot_length], col=2, lwd=4, lty=2)


#################
### CONCEPTOR ###
#################
aperture <- 50
C <- computeConceptor(R, aperture)

output_C_collector <- matrix(0,M,n_learn)
x <- start_X
for(n in (n_washout+1):(n_learn+n_washout)){
  u <- pattern[n,]
  x_old <- x
  x <- tanh(W_star%*%x_old + W_in %*% u + b)
  x <- C %*% x
  output_C_collector[,n-n_washout] <- W_out %*% x
}
conceptors_output <- output_C_collector

# NRMSE
error <- nrmse(original_outputs, conceptors_output)
print(paste0('nrmse : ', round(error,5)))

# plot the original patterns, the trained output, and the conceptor patterns
plot(original_outputs[1:plot_length], type = 'l', ylim = c(-1,1), ylab = '')
lines(conceptors_output[1:plot_length], col=2, lty=2, lwd=4)

#####################
### BEETJE KLOTEN ###
#####################
svd_C <- svd(C)
plot(svd_C$d, type = 'l')
lines(cumsum(svd_C$d)/sum(svd_C$d), col=2)
lines(c(1,N), rep(0.95,2), col=3, lty=2)

##################################
### COMPUTE F GRADIENT DESCENT ###
##################################
apertures_2 <- aperture^-2
K <- 100
n_washout <- 100
n_adapt <- 500
n_learn <- n_points - n_adapt - n_washout
learning_rate_D <- 0.5
learning_rate_F <- 0.1
spectral_radius_scaling <- 1

# initial states
z <- matrix(0,K,1)
D <- diag(rep(1,K))
F_raw <- createInitialForwardProjectionWeights() # dim=KxN 
spectral_radius_FT_F <- spectralRadius(t(F_raw) %*% F_raw)
F <- (F_raw / sqrt(spectral_radius_FT_F)) * sqrt(spectral_radius_scaling)

Ft <- t(F)
a <- norm(F, type='F')^-2
FtDF <- Ft%*%D%*%F
aDF <- a*D%*%F
# create progress bar
pb <- txtProgressBar(min = 1, max = n_washout, style = 3)
for(n in 1:n_washout){
  u <- pattern[n,]
  z_old <- z
  r <- tanh(W_star %*% Ft %*% z_old + W_in %*% u + b)
  z <- aDF %*% r
  # update progress bar
  setTxtProgressBar(pb, n)
}

loss <- rep(NA,n_adapt)
# create progress bar
pb <- txtProgressBar(min = (n_washout+1), max = (n_adapt+n_washout), style = 3)
for(n in (n_washout+1):(n_adapt+n_washout)){
  u <- pattern[n,]
  z_old <- z
  Ft <- t(F)
  a <- norm(F, type='F')^-2
  FtDF <- Ft%*%D%*%F
  loss[n - n_washout] <- norm(C - a*FtDF, type = 'F')^2
  aDF <- a*D%*%F
  r <- tanh(W_star %*% Ft %*% z_old + W_in %*% u + b)
  z <- aDF %*% r
  #plot(z, type = 'l')
  #Sys.sleep(0.1)
  Dz2 <- diag(as.vector(z*z))
  D_adaptation <- learning_rate_D * (Dz2 - D%*%Dz2 - apertures_2*D)
  D <- D + D_adaptation
  F_adaptation <- 4 * learning_rate_F * a * (D%*%F%*%C - a*(tr(FtDF%*%C)*F + D%*%F%*%FtDF) + a^2 * tr(FtDF%*%FtDF)*F)
  F <- F + F_adaptation
  if(n==(n_adapt+n_washout)){
    D_final <- D
    F_final <- F
  }
  if(sum(as.integer(is.na(F)))>0){
    print(n)
    break
  }
  # update progress bar
  setTxtProgressBar(pb, n)
}
plot(loss, type = 'l')

Ft_final <- t(F_final)
a_final <- norm(F_final, type='F')^-2
aDF_final <- a_final * D_final %*% F_final

output_F_collector <- matrix(0,M,n_learn)
compare_output_collector <- matrix(0,M,n_learn)
# create progress bar
pb <- txtProgressBar(min = 1, max = n_learn, style = 3)
for(n in 1:n_learn){
  u <- pattern[n_washout+n_adapt+n,]
  z_old <- z
  r <- tanh(W_star %*% Ft_final %*% z_old + W_in %*% u + b)
  z <- aDF_final %*% r
  output_F_collector[,n] <- W_out %*% r
  compare_output_collector[,n] <- u
  # update progress bar
  setTxtProgressBar(pb, n)
}

# plot the original patterns, the trained output, and the conceptor patterns
plot(compare_output_collector[1:plot_length], type = 'l', ylim = c(-1,1), ylab = '')
lines(output_F_collector[1:plot_length], type = 'l', col=2, lty=2, lwd=4) 

# NRMSE
error <- nrmse(compare_output_collector, output_F_collector)
print(paste0('nrmse : ', round(error,5)))
