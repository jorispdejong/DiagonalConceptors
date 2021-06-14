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
file_path <- 'data/periodic patterns/diagonal conceptors/patterns.RData'
if(file.exists(file_path)){
  load(file_path)
}else{
  patterns <- createExamplePatterns(generate_new_data=T, L=1000, 
                                    folder_path='data/periodic patterns/diagonal conceptors/')
  save(patterns, file=file_path) 
}
# set pattern variables
n_pattern <- length(patterns) # number of patterns
n_points <- sapply(patterns, FUN = function(x) dim(x)[1]) # length of each pattern

# plot (part of) the patterns
par(mfrow=c(2,2))
n_plot <- 40
for(p in 1:4) plot(patterns[[p]][1:n_plot], type = 'l', ylim = c(-1,1), ylab = '')

##################
### PARAMETERS ###
##################
# scaling
W_in_scaling <- 1
b_scaling <- 0.2
W_star_scaling <- 1

# leaking rate
leaking_rate <- 1

# apertures
apertures <- rep(8,n_pattern)
apertures_2 <- apertures^-2

# regularization
reg_out <- 0
reg_W <- 1e-3
compute_D <- F
if(compute_D) reg_D <- 1e-3

# other parameters
verbose <- T

#######################
### RESERVOIR SETUP ###
#######################
# dimensions
M <- ncol(patterns[[1]]) # dimension of output 
N <- 100 # reservoir size

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
n_adapt <- rep(500,n_pattern)
n_learn <- rep(500,n_pattern)

#######################
### DRIVE RESERVOIR ###
#######################
results <- list()
for(sim in 1:2){
  # diagonal conceptors (sim=1: initialized randomly, sim=2: identity)
  if(sim == 1){
    cs <- lapply(1:n_pattern, FUN=function(p) matrix(1,N,1))
  }else{
    cs <- lapply(1:n_pattern, FUN=function(p) matrix(runif(N),N,1))
  }
  
  # initialize empty collectors
  all_training_r <- c()
  all_z_adapt <- list()
  all_z_learn <- list()
  all_training_z_old <- c()
  all_training_D_target <- c()
  all_training_W_target <- c()
  all_training_output <- c()
  # collect data from driving the reservoir with different drivers
  for(p in 1:n_pattern){
    pattern <- patterns[[p]]
    if(verbose) print(paste('pattern',p))
    
    # adapt period collectors
    z_collector_adapt <- matrix(0,N,n_adapt[p])
    z_collector_learn <- matrix(0,N,n_learn[p])
    
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
        z_collector_adapt[,n-n_washout[p]] <- z
        
        # compute diagonal conceptors
        if(n==n_washout[p]+n_adapt[p]){
          R_adj <- rowMeans(z_collector_adapt^2)
          cs[[p]] <- R_adj*(R_adj+apertures_2[p])^-1
        } 
      }
      
      if(n > n_washout[p]+n_adapt[p]){
        # collect states, output and G_target for computing W_out and G later
        z_old_collector[,n-n_washout[p]-n_adapt[p]] <- z_old
        W_target_collector[,n-n_washout[p]-n_adapt[p]] <- W_target

        # collect state for later comparison
        z_collector_learn[,n-n_washout[p]-n_adapt[p]] <- z
      }
      
      # update progress bar
      if(verbose) setTxtProgressBar(pb, n)
    }
    
    # construct concatenated matrices
    all_training_z_old <- cbind(all_training_z_old, z_old_collector)
    all_training_W_target <- cbind(all_training_W_target, W_target_collector)

    # states
    all_z_adapt[[p]] <- z_collector_adapt
    all_z_learn[[p]] <- z_collector_learn
  }
  
  #############################
  ### (RE)COMPUTING WEIGHTS ###
  #############################
  
  cp_W <- recomputeReservoirWeights(all_training_states_old = all_training_z_old, 
                                    all_training_W_target = all_training_W_target,
                                    reg_W = reg_W, verbose = verbose)
  W <- cp_W$W
  
  ##############################
  ### SELF-GENERATE PATTERNS ###
  ##############################
  # let reservoir self-generate with diagonal conceptors
  n_run <- n_learn
  n_washout2 <- 2*n_washout
  all_z_run <- list()
  for(p in 1:n_pattern){
    pattern <- patterns[[p]]
    if(verbose) print(paste('pattern',p))
    
    z_collector <- matrix(0,N,n_run[p])
    
    # initial state
    r <- matrix(0,N,1)
    z <- cs[[p]] * r
    
    # create progress bar
    if(verbose) pb <- txtProgressBar(1,n_run[p]+n_washout2[p], style = 3)
    for(n in 1:(n_run[p]+n_washout2[p])){
      r <- tanh(W %*% z + b)
      z <- (1-leaking_rate) * z + leaking_rate * cs[[p]]*r
      
      if(n > n_washout2[p]) z_collector[,n-n_washout2[p]] <- z
      
      # update progress bar
      if(verbose) setTxtProgressBar(pb, n)
    }
    
    all_z_run[[p]] <- z_collector
  }
  
  results[[sim]] <- list(all_z_adapt=all_z_adapt,
                         all_z_learn=all_z_learn,
                         all_z_run=all_z_run)
}


####################
### PLOT RESULTS ###
####################
# plot states
par(mfrow=c(1,2), mar=c(1,1,3,1), oma=rep(0,4))
dims <- c(10,45)
cex <- 1.2
main_titles <- c(expression('D'[0]*' - identity'), expression('D'[0]*' - random'))

for(sim in 1:2){
  z_adapt <- results[[sim]]$all_z_adapt
  z_learn <- results[[sim]]$all_z_learn
  
  ### D0 - Identity
  # empty plot
  plot(NA,NA, 
       xlim = c(-1,1), ylim = c(-1,1), 
       xaxt='n', yaxt='n',
       main = main_titles[sim], cex.main=1.5)
  # pattern 1
  ellips <- sortXY(z_adapt[[1]][dims[1],], z_adapt[[1]][dims[2],])
  points(ellips$x, ellips$y, type = 'b', lwd=5)

  # pattern 2
  ellips <- sortXY(z_adapt[[2]][dims[1],], z_adapt[[2]][dims[2],])
  points(ellips$x, ellips$y, type = 'b', col=2, lwd=5)
  
  legend('topright', legend = c(expression('D'[0]*'r'^1),
                                expression('D'[0]*'r'^2)), 
         col=1:2, pch=16, cex=1.2)
}




