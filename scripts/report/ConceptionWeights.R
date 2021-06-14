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
file_path <- 'data/periodic patterns/diagonal conceptors/patterns5000.RData'
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

# apertures
apertures <- rep(8,n_pattern)
apertures_2 <- apertures^-2

# leaking rate
leaking_rate <- 1

# regularization
reg_out <- 0
reg_W <- 1e-3

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
n_adapt <- rep(3000,n_pattern)
n_learn <- n_points - n_washout - n_adapt

#######################
### DRIVE RESERVOIR ###
#######################
# diagonal conceptors initialized randomly
cs_init <- lapply(1:n_pattern, FUN=function(p) matrix(runif(N),N,1))

results <- list()
for(sim in 1:2){
  if(sim == 1){
    iterative <- F
  }else{
    iterative <- T
  }
  
  # initialize empty collectors
  all_z_adapt <- list()
  all_z_learn <- list()
  all_training_r <- c()
  all_training_z_old <- c()
  all_training_D_target <- c()
  all_training_W_target <- c()
  all_training_output <- c()
  
  if(iterative){
    learning_rate <- 0.5
    all_cs_collector <- list()
  }
  # diagonal conceptors
  cs <- cs_init
  # collect data from driving the reservoir with different drivers
  for(p in 1:n_pattern){
    pattern <- patterns[[p]]
    if(verbose) print(paste('pattern',p))
    
    # adapt period collectors
    z_collector_adapt <- matrix(0,N,n_adapt[p])
    if(iterative) cs_collector <- matrix(0,N,n_adapt[p])
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
      z <- (1-leaking_rate) * z_old + leaking_rate * cs[[p]] * r
      
      if(n > n_washout[p] && n <= n_washout[p]+n_adapt[p]){
        # collect z states
        z_collector_adapt[,n-n_washout[p]] <- z
        
        if(iterative){
          if(n > n_washout[p] + 0.5*n_adapt[p]){
            #learning_rate <- learning_rate * (1-((n-n_washout[p]-ceiling(0.5*n_adapt[p])) / ceiling(0.5*n_adapt[p])))
          }
          
          cs[[p]] <- cs[[p]] + learning_rate * ((1-cs[[p]])*z^2 - apertures_2[p]*cs[[p]])
          cs_collector[,n-n_washout[p]] <- cs[[p]]
          
        }else{
          # compute diagonal conceptors
          if(n==n_washout[p]+n_adapt[p]){
            R_adj <- rowMeans(z_collector_adapt^2)
            cs[[p]] <- R_adj*(R_adj+apertures_2[p])^-1
          } 
        }
      }
      
      if(n > n_washout[p]+n_adapt[p]){
        # collect states, output and G_target for computing W_out and G later
        r_collector[,n-n_washout[p]-n_adapt[p]] <- r
        z_old_collector[,n-n_washout[p]-n_adapt[p]] <- z_old
        D_target_collector[,n-n_washout[p]-n_adapt[p]] <- D_target
        W_target_collector[,n-n_washout[p]-n_adapt[p]] <- W_target
        output_collector[,n-n_washout[p]-n_adapt[p]] <- u
        
        # collect state for later comparison
        z_collector_learn[,n-n_washout[p]-n_adapt[p]] <- z
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
    
    # states
    all_z_adapt[[p]] <- z_collector_adapt
    all_z_learn[[p]] <- z_collector_learn
    
    # cs collector
    if(iterative) all_cs_collector[[p]] <- cs_collector
  }
  
  # split training outputs matrix into a list
  target_outputs <- splitMatrixToList(t(all_training_output), lengths = n_learn, splitAlong = 'rows')
  
  #############################
  ### (RE)COMPUTING WEIGHTS ###
  #############################
  # compute the output weights
  cp_W_out <- computeOutputWeights(all_training_states = all_training_r, 
                                   all_training_output = all_training_output,
                                   reg_out = reg_out, verbose = verbose)
  W_out <- cp_W_out$W_out
  
  # recompute reservoir weights (storing the patterns into the reservoir)
  cp_W <- recomputeReservoirWeights(all_training_states_old = all_training_z_old, 
                                    all_training_W_target = all_training_W_target,
                                    reg_W = reg_W, verbose = verbose)
  W <- cp_W$W
  
  ##############################
  ### SELF-GENERATE PATTERNS ###
  ##############################
  # let reservoir self-generate with diagonal conceptors
  n_run <- n_points
  n_washout2 <- 2*n_washout
  all_z_run <- list()
  all_z_run_washout <- list()
  sg_outputs <- list()
  for(p in 1:n_pattern){
    pattern <- patterns[[p]]
    if(verbose) print(paste('pattern',p))
    
    output_collector <- matrix(0,M,n_run[p])
    z_washout_collector <- matrix(0,N,n_washout2[p])
    z_collector <- matrix(0,N,n_run[p])
    
    # initial state
    r <- matrix(1,N,1)
    z <- cs[[p]] * r
    
    # create progress bar
    if(verbose) pb <- txtProgressBar(1,n_run[p]+n_washout2[p], style = 3)
    for(n in 1:(n_run[p]+n_washout2[p])){
      r <- tanh(W %*% z + b)
      z <- (1-leaking_rate) * z + leaking_rate * cs[[p]] * r
      
      if(n > n_washout2[p]){
        z_collector[,n-n_washout2[p]] <- z
        output_collector[,n-n_washout2[p]] <- W_out %*% r
      }else{
        z_washout_collector[,n] <- z
      }
      # update progress bar
      if(verbose) setTxtProgressBar(pb, n)
    }
    
    all_z_run[[p]] <- z_collector
    all_z_run_washout[[p]] <- z_washout_collector
    sg_outputs[[p]] <- t(output_collector)
  }
  
  # shift the output of the sg_outputs for comparison
  shifted_patterns <- lapply(1:n_pattern, function(p) shiftPattern(sg_outputs[[p]], 
                                                                   patterns[[p]],
                                                                   max_shift = 100,
                                                                   shift_along = 'rows'))
  
  if(iterative){
    results[[sim]] <- list(all_z_adapt=all_z_adapt,
                           all_z_learn=all_z_learn,
                           all_z_run=all_z_run, 
                           all_z_run_washout=all_z_run_washout,
                           all_cs_collector=all_cs_collector,
                           cs=cs,
                           shifted_patterns=shifted_patterns)
  }else{
    results[[sim]] <- list(all_z_adapt=all_z_adapt,
                           all_z_learn=all_z_learn,
                           all_z_run=all_z_run, 
                           all_z_run_washout=all_z_run_washout,
                           cs=cs,
                           shifted_patterns=shifted_patterns) 
  }
}

# show nrmses
nrmses <- sapply(results, function(y) lapply(y$shifted_patterns, function(x) x$nrmse))
nrmses

####################
### PLOT RESULTS ###
####################
show_iterative_states <- F
if(show_iterative_states){
  # plot states
  dims <- c(10,45)
  cex <- 1.2
  main_titles <- c('Explicit', 'Iterative')
  n_plot <- 500
  
  par(mfrow=c(1,2), mar=c(1,1,3,1), oma=rep(0,4))
  for(sim in 1:2){
    z_learn <- results[[sim]]$all_z_learn
    
    ### D0 - Identity
    # empty plot
    plot(NA,NA, 
         xlim = c(-1,1), ylim = c(-1,1), 
         xaxt='n', yaxt='n',
         main = main_titles[sim], cex.main=1.5)
    # pattern 1
    ellips <- sortXY(z_learn[[1]][dims[1],1:n_plot], z_learn[[1]][dims[2],1:n_plot])
    points(ellips$x, ellips$y, type = 'b', lwd=5)
    
    # pattern 2
    ellips <- sortXY(z_learn[[2]][dims[1],1:n_plot], z_learn[[2]][dims[2],1:n_plot])
    points(ellips$x, ellips$y, type = 'b', col=2, lwd=5)
    
    legend('topright', legend = c(expression('D'[0]*'r'^1),
                                  expression('D'[0]*'r'^2)), 
           col=1:2, pch=16, cex=1.2)
  }
  
}

# plot convergence
show_convergence <- T
if(show_convergence){
  colfunc<-colorRampPalette(c("gray","black"))
  x <- t(results[[2]]$all_cs_collector[[1]])
  
  par(mfrow=c(1,1), mar=c(3,3,3,1), oma=c(0,0,0,0))
  matplot(x[,order(x[nrow(x),])], type = 'l',
          xlab = '', ylab = '', yaxt = 'n',
          col = colfunc(N),
          main = 'Convergence of conception weights', 
          cex.main=1.5, cex.axis=1.2)
  axis(side = 2, at = c(0,0.5,1), cex.axis=1.2)
}

# plot conception weights
show_conception_weights <- T
if(show_conception_weights){
  # conception weights
  par(mfrow=c(1,2), mar=c(3,3,2.5,1), oma=c(0,0,1.5,0))
  for(sim in 1:2){
    plot(results[[sim]]$cs[[1]], pch=16, cex=1.2,
         ylim = c(0,1), xaxt='n', yaxt='n', xlab='', ylab='',
         main = ifelse(sim==1, 'Explicit','Iterative'))
    axis(side = 2, at=c(0,0.5,1), cex.axis=1.2)
    mtext(expression('c'['i']), side=1, line=1.5, cex=1.8)
    abline(h=0, col=2, lty=2, lwd=2)
    abline(h=0.5, col=2, lty=2, lwd=2)
    abline(h=1, col=2, lty=2, lwd=2)
  }
  title('Conception weights', outer = T, cex.main = 2, line=-0.3) 
}
