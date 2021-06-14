#############################
### DRIVING THE RESERVOIR ###
#############################
### WASHOUT PERIOD ###
driveReservoirWashoutPeriod <- function(patterns,
                                        cs_init,
                                        W_in, b, W_star, leaking_rate,
                                        n_washout_max=100, tol=1e-5, 
                                        state_init1=matrix(0,N,1), state_init2=matrix(1,N,1),
                                        verbose=F){
  ### INPUT: patterns = list of patterns
  #          cs_init = list of initial diagonal conceptors
  #          W_in = input weights
  #          b = bias weights
  #          W_star = reservoir weights
  #          leaking_rate = leaking rate
  #          n_washout_max = integer or vector of integers containing the max length of the washout period
  #          tol = tolerance that determines when the two initial state vectors have converged
  #          state_init1 = initial state vector 1
  #          state_init2 = initial state vector 2 (should be different from initial state vector 1)
  #
  ### OUTPUT: n_washout = vector containing the found washout length after which the state vector had converged
  #           final_cz_washout = the final state of the system (this will be used for continuing driving the reservoir)
  
  
  ### CHECKS ###
  # patterns is a list containing matrices?
  if(!is.list(patterns)){
    stop('patterns must be a list of patterns, which are matrices of dimension LxM')
  }else{
    if(length(unique(sapply(patterns, function(x) ncol(x))))!=1) stop('patterns must all have the same number of columns')
  }
  # n_pattern is defined?
  if(!exists('n_pattern')) n_pattern <- length(patterns)
  # n_washout_max is an array?
  if(length(n_washout_max)==1){
    n_washout_max <- rep(n_washout_max,n_pattern) 
  }else{
    if(length(n_washout_max)!=n_pattern) stop(paste0('length of n_washout_max does not agree with the correct length, ', n_pattern))
  }
  # state_init1 and state_init2 are equal?
  if(all(state_init1==state_init2)){
    stop('state_init1 and starting_state2 are equal, but they must be different otherwise the algorithm is useless')
  }else{
    if(nrow(state_init1)!=N || ncol(state_init1)!=1) stop('state_init1 does not have the correct dimensions. It must be a matrix of dimension Nx1')
    if(nrow(state_init2)!=N || ncol(state_init2)!=1) stop('state_init2 does not have the correct dimensions. It must be a matrix of dimension Nx1')
  }
  # cs_init is a list of n_pattern diagonal conceptors?
  if(length(cs_init)!=n_pattern){
    stop(paste0('cs_init must be a list of ',n_pattern, ' matrices of dimension Nx1')) 
  }else{
    if(!all(sapply(cs_init, function(x) nrow(x)==N && ncol(x)==1))) stop(paste0('cs_init must be a list of ',n_pattern, ' matrices of dimension Nx1'))
  }
  
  # how many first training points we discard maximally
  n_washout <- rep(NA, n_pattern)
  
  # initialize empty collectors
  final_cz_washout <- list()
  # collect data from driving the reservoir with different drivers
  for(p in 1:n_pattern){
    pattern <- patterns[[p]]
    if(verbose) print(paste('pattern',p))
    
    # two different initial state
    z1 <- state_init1
    z2 <- state_init2
    cz1 <- cs_init[[p]] * z1
    cz2 <- cs_init[[p]] * z2
    
    # create progress bar
    if(verbose) pb <- txtProgressBar(1,n_washout_max[p])
    for(n in 1:n_washout_max[p]){
      u <- pattern[n,] 
      
      cz1_old <- cz1
      cz2_old <- cz2
      
      D <-  W_in %*% u + b
      z1 <- tanh(W_star %*% cz1_old + D)
      z2 <- tanh(W_star %*% cz2_old + D)
      
      cz1 <- (1-leaking_rate) * cz1_old + leaking_rate * cs_init[[p]]*z1
      cz2 <- (1-leaking_rate) * cz2_old + leaking_rate * cs_init[[p]]*z2
      
      distance <- sqrt(sum((cz1-cz2)^2))
      if(distance < tol){
        final_cz_washout[[p]] <- cz1
        n_washout[p] <- n
        break
      }
      if(n==n_washout_max[p]){
        final_cz_washout[[p]] <- cz1
        n_washout[p] <- n
        if(verbose) print(paste0('State vectors not converged within ', n_washout_max[p], 
                                 ' steps. Difference was ', round(distance,7)))
      }
      
      # set progress bar
      if(verbose) setTxtProgressBar(pb,n)
    }
  }
  
  return(list('n_washout'=n_washout, 'final_cz_washout'=final_cz_washout))
}
### ADAPTATION PERIOD ###
driveReservoirAdaptPeriod <- function(patterns,
                                      cs_init,
                                      W_in, b, W_star, leaking_rate,
                                      apertures,
                                      start_cz, n_washout,
                                      n_adapt_max=500, tol=1e-2,
                                      verbose=F){
  ### INPUT: patterns = list of patterns
  #          cs_init = list of initial diagonal conceptors
  #          W_in = input weights
  #          b = bias weights
  #          W_star = reservoir weights
  #          leaking_rate = leaking rate
  #          apertures = number or vector containing the apertures
  #          start_cz = list with the starting states of the reservoir, which are the final states of the washout period
  #          n_washout = vector with the washout lengths found in the washout period
  #          n_adapt_max = integer or vector of integers containing the max length of the adaption period
  #          tol = tolerance that determines when the diagonal conceptor stop changing
  #
  ### OUTPUT: n_adapt = vector containing the found adaptation lengths after which the diagonal conceptors had stopped changing
  #           final_cz_adapt = the final state of the system (this will be used for continuing driving the reservoir)
  
  ### CHECKS ###
  # patterns is a list containing matrices?
  if(!is.list(patterns)){
    stop('patterns must be a list of patterns, which are matrices of dimension LxM')
  }else{
    if(length(unique(sapply(patterns, function(x) ncol(x))))!=1) stop('patterns must all have the same number of columns')
  }
  # n_pattern is defined?
  if(!exists('n_pattern')) n_pattern <- length(patterns)
  # n_adapt_max is an array?
  if(length(n_adapt_max)==1){
    n_adapt_max <- rep(n_adapt_max,n_pattern) 
  }else{
    if(length(n_adapt_max)!=n_pattern) stop(paste0('length of n_adapt_max does not agree with the correct length, ', n_pattern))
  }
  # apertures is an array?
  if(length(apertures)==1){
    apertures <- rep(apertures,n_pattern) 
  }else{
    if(length(apertures)!=n_pattern) stop(paste0('length of apertures does not agree with the correct length, ', n_pattern))
  }
  # cs_init is a list of n_pattern diagonal conceptors?
  if(length(cs_init)!=n_pattern){
    stop(paste0('cs_init must be a list of ',n_pattern, ' matrices of dimension Nx1')) 
  }else{
    if(!all(sapply(cs_init, function(x) nrow(x)==N && ncol(x)==1))) stop(paste0('cs_init must be a list of ',n_pattern, ' matrices of dimension Nx1'))
  }
  
  # how many first training points we discard maximally
  n_adapt <- rep(NA, n_pattern)
  
  # initialize empty collectors
  cs <- list()
  final_cz_adapt <- list()
  # collect data from driving the reservoir with different drivers
  for(p in 1:n_pattern){
    pattern <- patterns[[p]]
    if(verbose) print(paste('pattern',p))
    
    # set starting cz
    cz <- start_cz[[p]]
    
    # collector that will collect the computed conceptor vectors
    cs_adapt_collector <- cs_init[[p]]
    # collector that collects the cz states that are used to compute the conceptor vectors
    cz_collector_adapt <- matrix(cz,N,1)
    
    # create progress bar
    if(verbose) pb <- txtProgressBar(n_washout[p]+1,n_washout[p]+n_adapt_max[p])
    for(n in (n_washout[p]+1):(n_washout[p]+n_adapt_max[p])){
      u <- pattern[n,] 
      cz_old <- cz
      z <- tanh(W_star %*% cz_old + W_in %*% u + b)
      cz <- (1-leaking_rate) * cz_old + leaking_rate * cs_init[[p]]*z
      
      # collect cz
      cz_collector_adapt <- cbind(cz_collector_adapt, cz)
      
      if(n > n_washout[p]+1){
        R_adj <- rowMeans(cz_collector_adapt^2)
        cs_adapt_collector <- cbind(cs_adapt_collector, R_adj*(R_adj+apertures[p]^-2)^-1)
        
        # compute difference with last conceptor vector
        distance <- sqrt(sum((cs_adapt_collector[,ncol(cs_adapt_collector)] - cs_adapt_collector[,ncol(cs_adapt_collector)-1])^2))
        if(distance < tol){
          cs[[p]] <- as.matrix(cs_adapt_collector[,ncol(cs_adapt_collector)])
          final_cz_adapt[[p]] <- cz
          n_adapt[p] <- ncol(cs_adapt_collector)
          break
        }
        if(n==n_washout[p]+n_adapt_max[p]){
          cs[[p]] <- as.matrix(cs_adapt_collector[,ncol(cs_adapt_collector)])
          final_cz_adapt[[p]] <- cz
          n_adapt[p] <- ncol(cs_adapt_collector)
          if(verbose) print(paste0('Conceptor vectors not converged within ', n_adapt_max[p], 
                                   ' steps. Difference was ', round(distance,5)))
        }
      }
      
      # set progress bar
      if(verbose) setTxtProgressBar(pb,n)
    }
  }
  
  return(list('cs'=cs, 'n_adapt'=n_adapt, 'final_cz_adapt'=final_cz_adapt))
}
### LEARNING PERIOD ###
driveReservoirLearnPeriod <- function(patterns,
                                      cs,
                                      W_in, b, W_star, leaking_rate,
                                      start_cz,
                                      n_washout, n_adapt, n_learn,
                                      verbose=F){
  ### INPUT: patterns = list of patterns
  #          cs = list of diagonal conceptors
  #          W_in = input weights
  #          b = bias weights
  #          W_star = reservoir weights
  #          leaking_rate = leaking rate
  #          start_cz = list with the starting states of the reservoir which are the final states of the adaptation period
  #          n_learn = number or vector with the learning lengths
  #
  ### OUTPUT: all_training_z = concatenated matrix with the state vectors of the run
  #           all_training_cz_old = concatenated matrix with the conceptor multiplied state vectors of the run, delayed by 1 time step
  #           all_training_output = concatenated matrix with the outputs of the run
  #           all_training_W_target = concatenated matrix with the target vectors for computing the matrix W later
  
  ### CHECKS ###
  # patterns is a list containing matrices?
  if(!is.list(patterns)){
    stop('patterns must be a list of patterns, which are matrices of dimension LxM')
  }else{
    if(length(unique(sapply(patterns, function(x) ncol(x))))!=1) stop('patterns must all have the same number of columns')
  }
  # n_pattern is defined?
  if(!exists('n_pattern')) n_pattern <- length(patterns)
  # n_learn is an array?
  if(length(n_learn)==1){
    n_learn <- rep(n_learn,n_pattern) 
  }else{
    if(length(n_learn)!=n_pattern) stop(paste0('length of n_learn does not agree with the correct length, ', n_pattern))
  }
  # cs is a list of n_pattern diagonal conceptors?
  if(length(cs)!=n_pattern){
    stop(paste0('cs must be a list of ',n_pattern, ' matrices of dimension Nx1')) 
  }else{
    if(!all(sapply(cs, function(x) nrow(x)==N && ncol(x)==1))) stop(paste0('cs must be a list of ',n_pattern, ' matrices of dimension Nx1'))
  }
  
  
  # initialize empty collectors
  all_training_z <- c()
  all_training_cz_old <- c()
  all_training_W_target <- c()
  all_training_output <- c()
  # collect data from driving the reservoir with different drivers
  for(p in 1:n_pattern){
    pattern <- patterns[[p]]
    if(verbose) print(paste('pattern',p))
    
    # empty collectors
    z_collector <- matrix(0,N,n_learn[p])
    cz_old_collector <- matrix(0,N,n_learn[p])
    W_target_collector <- matrix(0,N,n_learn[p])
    output_collector <- matrix(0,M,n_learn[p])
    
    # initial state
    cz <- start_cz[[p]]
    
    # create progress bar
    if(verbose) pb <- txtProgressBar(n_washout[p]+1+n_adapt[p]+1,n_washout[p]+n_adapt[p]+n_learn[p])
    for(n in (n_washout[p]+n_adapt[p]+1):(n_washout[p]+n_adapt[p]+n_learn[p])){
      u <- pattern[n,]
      cz_old <- cz
      W_target <- W_star %*% cz_old + W_in %*% u
      z <- tanh(W_target + b)
      cz <- (1-leaking_rate) * cz_old + leaking_rate * cs[[p]]*z
      
      # collect states, output and G_target for computing W_out and G later
      z_collector[,n-n_washout[p]-n_adapt[p]] <- z
      cz_old_collector[,n-n_washout[p]-n_adapt[p]] <- cz_old
      W_target_collector[,n-n_washout[p]-n_adapt[p]] <- W_target
      output_collector[,n-n_washout[p]-n_adapt[p]] <- u
      
      # update progress bar
      if(verbose) setTxtProgressBar(pb, n)
    }
    
    # construct concatenated matrices
    all_training_z <- cbind(all_training_z, z_collector)
    all_training_cz_old <- cbind(all_training_cz_old, cz_old_collector)
    all_training_output <- cbind(all_training_output, output_collector)
    all_training_W_target <- cbind(all_training_W_target, W_target_collector)
  }
  
  # split training outputs matrix into a list
  target_outputs <- splitMatrixToList(t(all_training_output), lengths = n_learn, splitAlong = 'rows')
  
  return(list('all_training_z'=all_training_z, 'all_training_cz_old'=all_training_cz_old,
              'all_training_output'=all_training_output, 'all_training_W_target'=all_training_W_target,
              'target_outputs'=target_outputs))
}

#################################
### SELF GENERATING RESERVOIR ###
#################################
selfGeneratingReservoir <- function(cs,
                                    start_cz=NULL,
                                    n_washout=100,
                                    W, W_out, b, leaking_rate,
                                    n_run, 
                                    verbose=F){
  ### INPUT: cs = list of diagonal conceptors
  #          start_cz = list containing the starting state of the reservoir
  #          W = recomputed reservoir weights
  #          W_out = output weights
  #          b = bias weights
  #          leaking_rate = leaking rate
  #          n_run = integer or vector of integers containing the length of self-generation period
  #
  ### OUTPUT: DC_outputs = list containing the outputs of the self-generated run of the reservoir
  #
  
  ### CHECKS ###
  # n_pattern is defined?
  if(!exists('n_pattern')) n_pattern <- length(cs)
  # n_run is an array?
  if(length(n_run)==1){
    n_run <- rep(n_run,n_pattern) 
  }else{
    if(length(n_run)!=n_pattern) stop(paste0('length of n_run does not agree with the correct length, ', n_pattern))
  }
  # start_cz is NULL? then we must do a washout period
  n_wash <- rep(0, n_pattern)
  if(is.null(start_cz)){
    if(length(n_washout)==1){
      n_wash <- rep(n_washout, n_pattern)
    }else if(length(n_washout)!=n_pattern){
      stop(paste0('length of n_washout is not correct, it should be ',n_pattern))
    }else{
      n_wash <- n_washout
    }
    start_cz <- lapply(1:n_pattern, function(x) matrix(0,N,1))
  }
  
  sg_outputs <- list()
  for(p in 1:n_pattern){
    if(verbose) print(paste0('Pattern ',p))
    # empty collectors
    p_collector <- matrix(0,M,n_run[p])
    
    # initial states
    cz <- start_cz[[p]]
    
    # create progress bar
    if(verbose) pb <- txtProgressBar(1,n_run[p], style = 3)
    for(n in 1:(n_wash[p] + n_run[p])){
      # update equation
      cz_old <- cz
      z <- tanh(W %*% cz_old + b)
      cz <- (1-leaking_rate) * cz_old + leaking_rate * cs[[p]]*z
      
      # compute output
      if(n > n_wash[p]) p_collector[,n-n_wash[p]] <- W_out %*% z
      
      # update progress bar
      if(verbose) setTxtProgressBar(pb, n)
    }
    sg_outputs[[p]] <- t(p_collector)
  }
  
  return(list('sg_outputs'=sg_outputs))
}


#######################
### FULL SIMULATION ###
#######################
# store the patterns in a reservoir
storePatternsInReservoir <- function(patterns,
                                     cs_init,
                                     W_in_raw, b_raw, W_star_raw, 
                                     W_in_scaling, b_scaling, W_star_scaling,
                                     leaking_rate,
                                     apertures,
                                     reg_out, reg_W,
                                     n_washout_max=100, n_adapt_max=500,
                                     verbose=F,
                                     show_plots=T){
  # scale weights
  scaled_weights <- scaleWeights(W_in_raw = W_in_raw, W_in_scaling = W_in_scaling,
                                 b_raw = b_raw, b_scaling = b_scaling,
                                 W_star_raw = W_star_raw, W_star_scaling = W_star_scaling)
  W_in <- scaled_weights$W_in
  b <- scaled_weights$b
  W_star <- scaled_weights$W_star
  
  # drive reservoir for max n_washout_max time steps, but stop when two initial state vectors converge
  dr_washout <- driveReservoirWashoutPeriod(patterns = patterns,
                                            cs_init = cs_init,
                                            W_in = W_in, b = b, W_star = W_star, leaking_rate = leaking_rate,
                                            n_washout_max = n_washout_max,
                                            verbose = verbose)
  n_washout <- dr_washout$n_washout
  final_cz_washout <- dr_washout$final_cz_washout
  
  # drive reservoir for max n_adapt_max time steps, but stop when conceptors do not change any more
  dr_adapt <- driveReservoirAdaptPeriod(patterns = patterns,
                                        cs_init = cs_init,
                                        W_in = W_in, b = b, W_star = W_star, leaking_rate = leaking_rate,
                                        apertures = apertures,
                                        start_cz = final_cz_washout, n_washout = n_washout,
                                        n_adapt_max = n_adapt_max,
                                        verbose = verbose)
  cs <- dr_adapt$cs
  n_adapt <- dr_adapt$n_adapt
  start_cz <- dr_adapt$final_cz_adapt
  
  # drive reservoir for the n_learn time step and collect states
  n_learn <- n_points - n_washout - n_adapt
  dr_learn <- driveReservoirLearnPeriod(patterns = patterns,
                                        cs = cs,
                                        W_in = W_in, b = b, W_star = W_star, leaking_rate = leaking_rate,
                                        start_cz = start_cz,
                                        n_washout = n_washout, n_adapt = n_adapt, n_learn = n_learn,
                                        verbose = verbose)
  all_training_z <- dr_learn$all_training_z
  all_training_cz_old <- dr_learn$all_training_cz_old
  all_training_output <- dr_learn$all_training_output
  all_training_W_target <- dr_learn$all_training_W_target
  target_outputs <- lapply(dr_learn$target_outputs, as.matrix)
  
  # plot some neuron states
  if(show_plots) plotNeuronStates(X=all_training_z, n=20, l=50)
  
  # compute the output weights
  cp_W_out <- computeOutputWeights(all_training_states = all_training_z, 
                                   all_training_output = all_training_output,
                                   reg_out = reg_out, verbose = verbose)
  W_out <- cp_W_out$W_out
  
  # recompute reservoir weights (loading the patterns into the reservoir)
  cp_W <- recomputeReservoirWeights(all_training_states_old = all_training_cz_old, 
                                    all_training_W_target = all_training_W_target,
                                    reg_W = reg_W, verbose = verbose)
  W <- cp_W$W
  
  # let reservoir self-generate with diagonal conceptors
  n_run <- n_learn
  sg <- selfGeneratingReservoir(cs = cs,
                                n_washout = n_washout,
                                W = W, W_out = W_out, b = b, leaking_rate =  leaking_rate,
                                n_run = n_run, verbose = verbose)
  sg_outputs <- sg$sg_outputs
  
  # shift the output of the sg_outputs for comparison
  shifted_patterns <- lapply(1:n_pattern, function(p) shiftPattern(sg_outputs[[p]], 
                                                                   target_outputs[[p]], 
                                                                   shift_along = 'rows'))
  
  # compare target outputs and self-generated outputs (plotting and Normalized Root Mean Square Error (NRMSE))
  par(mfrow=c(2,2), mar=rep(3,4))
  n_plot <- 50
  nrmses <- rep(NA, n_pattern)
  for(p in 1:n_pattern){
    # plot outputs
    if(show_plots) plot(shifted_patterns[[p]]$target_pattern[1:n_plot], type = 'l', xlab = '', ylab = '', ylim = c(-1,1))
    if(show_plots) lines(shifted_patterns[[p]]$pattern[1:n_plot], col=2, lwd=3, lty=2)
    
    nrmses[p] <- shifted_patterns[[p]]$nrmse
    if(verbose) print(paste0('pattern ',p,' nrmse : ', round(nrmses[p],5)))
  }
  
  return(list('nrmses'=nrmses, 
              'reservoir'=list('cs'=cs,
                               'W'=W,
                               'b'=b,
                               'W_out'=W_out,
                               'start_cz'=start_cz,
                               'leaking_rate'=leaking_rate)))
}






