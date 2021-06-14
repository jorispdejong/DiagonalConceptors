#############################
### DRIVING THE RESERVOIR ###
#############################
### WASHOUT PERIOD ###
driveReservoirWashoutPeriod <- function(patterns,
                                        W_in, b, W_star, leaking_rate,
                                        n_washout_max=100, tol=1e-5, 
                                        state_init1=matrix(0,N,1), state_init2=matrix(1,N,1),
                                        verbose=F){
  ### INPUT: patterns = list of patterns
  #          W_in = input weights
  #          b = bias weights
  #          W = reservoir weights
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
  
  # how many first training points we discard maximally
  n_washout <- rep(NA, n_pattern)
  
  # initialize empty collectors
  final_x_washout <- list()
  # collect data from driving the reservoir with different drivers
  for(p in 1:n_pattern){
    pattern <- patterns[[p]]
    if(verbose) print(paste('pattern',p))
    
    # two different initial state
    x1 <- state_init1
    x2 <- state_init2
    
    # create progress bar
    if(verbose) pb <- txtProgressBar(1,n_washout_max[p], style = 3)
    for(n in 1:n_washout_max[p]){
      u <- pattern[n,]
      
      x1_old <- x1
      x2_old <- x2
      
      D <- W_in %*% u + b
      x1 <- (1-leaking_rate)*x1_old + leaking_rate*tanh(W_star %*% x1_old + D)
      x2 <- (1-leaking_rate)*x2_old + leaking_rate*tanh(W_star %*% x2_old + D)
      
      distance <- sqrt(sum((x1-x2)^2))
      if(distance < tol){
        final_x_washout[[p]] <- x1
        n_washout[p] <- n
        break
      }
      if(n==n_washout_max[p]){
        final_x_washout[[p]] <- x1
        n_washout[p] <- n
        print(paste0('State vectors not converged within ', n_washout_max[p], ' steps. Difference was ', round(distance,5)))
      }
      
      # set progress bar
      if(verbose) setTxtProgressBar(pb,n)
    }
  }
  
  return(list('n_washout'=n_washout, 'final_x_washout'=final_x_washout))
}
### LEARNING PERIOD ###
driveReservoirLearnPeriod <- function(patterns,
                                      W_in, b, W_star, 
                                      leaking_rate,
                                      start_x,
                                      n_washout, n_learn,
                                      verbose=F){
  ### INPUT: patterns = list of patterns
  #          W_in = input weights
  #          b = bias weights
  #          W = reservoir weights
  #          leaking_rate = leaking rate
  #          start_x = list with the starting states of the reservoir which are the final states of the washout period
  #          n_learn = number or vector with the learning lengths
  #
  ### OUTPUT: all_training_x = concatenated matrix with the state vectors of the run
  #           all_training_x_old = concatenated matrix with the conceptor multiplied state vectors of the run, delayed by 1 time step
  #           all_training_output = concatenated matrix with the outputs of the run
  #           all_training_W_target = concatenated matrix with the target vectors for computing the matrix G later
  
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
  
  all_training_x <- c()
  all_training_x_old <- c()
  all_training_W_target <- c()
  all_training_output <- c()
  R <- list()
  # collect data from driving the reservoir with different drivers
  for(p in 1:n_pattern){
    pattern <- patterns[[p]]
    if(verbose) print(paste('pattern',p))
    
    # empty collectors
    x_collector <- matrix(0,N,n_learn[p])
    x_old_collector <- matrix(0,N,n_learn[p])
    W_target_collector <- matrix(0,N,n_learn[p])
    output_collector <- matrix(0,M,n_learn[p])
    
    # starting state
    x <- start_x[[p]]
    
    # create progress bar
    if(verbose) pb <- txtProgressBar(n_washout[p]+1,n_washout[p]+n_learn[p])
    for(n in (n_washout[p]+1):(n_washout[p]+n_learn[p])){
      u <- pattern[n,]
      x_old <- x
      W_target <- W_star %*% x_old + W_in %*% u
      x <- (1-leaking_rate)*x_old + leaking_rate*tanh(W_target+b)
      
      # collect states, output and W_target for computing W_out and W later
      x_collector[,n-n_washout[p]] <- x
      x_old_collector[,n-n_washout[p]] <- x_old
      W_target_collector[,n-n_washout[p]] <- W_target
      output_collector[,n-n_washout[p]] <- u
      
      # update progress bar
      if(verbose) setTxtProgressBar(pb, n)
    }
    
    # correlation matrix
    R[[p]] <- (x_collector %*%t(x_collector))/n_learn[p]
    
    # construct concatenated matrices
    all_training_x <- cbind(all_training_x, x_collector)
    all_training_x_old <- cbind(all_training_x_old, x_old_collector)
    all_training_output <- cbind(all_training_output, output_collector)
    all_training_W_target <- cbind(all_training_W_target, W_target_collector)
  }
  
  # split training outputs matrix into a list
  target_outputs <- splitMatrixToList(t(all_training_output), lengths = n_learn, splitAlong = 'rows')
  
  return(list('all_training_x'=all_training_x, 'all_training_x_old'=all_training_x_old,
              'all_training_output'=all_training_output, 'all_training_W_target'=all_training_W_target,
              'target_outputs'=target_outputs, 'R'=R))
}

#################################
### SELF GENERATING RESERVOIR ###
#################################
selfGeneratingReservoir <- function(Cs,
                                    start_x=NULL,
                                    n_washout=100,
                                    W, W_out, b, leaking_rate,
                                    n_run,
                                    which_patterns=NA,
                                    verbose=F){
  ### INPUT: Cs = list of conceptor matrices
  #          start_x = list containing the starting state of the reservoir
  #          W = recomputed reservoir weights
  #          W_out = output weights
  #          b = bias weights
  #          leaking_rate = leaking rate
  #          n_run = integer or vector of integers containing the length of self-generation period
  #
  ### OUTPUT: C_outputs = list containing the outputs of the self-generated run of the reservoir
  #
  
  ### CHECKS ###
  # n_pattern is defined?
  if(!exists('n_pattern')) n_pattern <- length(Cs)
  # which_patterns is given?
  if(is.na(which_patterns)) which_patterns <- 1:n_pattern
  # n_run is an array?
  if(length(n_run)==1){
    n_run <- rep(n_run,n_pattern) 
  }else{
    if(length(n_run)!=n_pattern) stop(paste0('length of n_run does not agree with the correct length, ', n_pattern))
  }
  # start_cz is NULL? then we must do a washout period
  n_wash <- rep(0, n_pattern)
  if(is.null(start_x)){
    if(length(n_washout)==1){
      n_wash <- rep(n_washout, n_pattern)
    }else if(length(n_washout)!=n_pattern){
      stop(paste0('length of n_washout is not correct, it should be ',n_pattern))
    }else{
      n_wash <- n_washout
    }
    start_x <- lapply(1:n_pattern, function(x) matrix(runif(N),N,1))
  }
  
  sg_outputs <- list()
  for(p in 1:n_pattern){
    if(p %in% which_patterns){
      if(verbose) print(paste0('Pattern ',p))
      # empty collectors
      output_collector <- matrix(0,M,n_run[p])
      
      # initial states
      x <- start_x[[p]]
      
      # create progress bar
      if(verbose) pb <- txtProgressBar(1,n_run[p], style = 3)
      for(n in 1:(n_wash[p] + n_run[p])){
        x_old <- x
        x <- (1-leaking_rate)*x_old + leaking_rate * Cs[[p]] %*% tanh(W %*% x + b)
        
        # compute output
        if(n > n_wash[p]) output_collector[,n-n_wash[p]] <- W_out %*% x
        
        # update progress bar
        if(verbose) setTxtProgressBar(pb, n)
      }
      sg_outputs[[p]] <- t(output_collector) 
    }else{
      sg_outputs[[p]] <- NA
    }
  }
  
  return(list('sg_outputs'=sg_outputs))
}
selfGeneratingReservoirConceptors <- function(R,
                                              apertures,
                                              start_x=NULL,
                                              n_washout=100,
                                              W, W_out, b, leaking_rate,
                                              n_run,
                                              which_patterns=NA,
                                              verbose=F){
  # which_patterns is given?
  if(length(which_patterns)==1){
    if(is.na(which_patterns)) which_patterns <- 1:n_pattern 
  }
  
  
  if(length(apertures)==1){
    apertures <- rep(apertures, length(R))
  }else{
    if(length(apertures)!=length(R)) stop('length of R and apertures must be equal')
  }
  
  # compute conceptors
  Cs <- list()
  if(verbose) pb <- txtProgressBar(1,length(R), style = 3)
  if(verbose) print('Computing conceptors...')
  for(p in 1:length(R)){
    if(p %in% which_patterns) Cs[[p]] <- computeConceptor(R[[p]],apertures[p])
    else Cs[[p]] <- NA
    
    if(verbose) setTxtProgressBar(pb,p)
  }
  
  
  # let reservoir self-generate with diagonal conceptors
  n_run <- n_learn
  sg <- selfGeneratingReservoir(Cs = Cs,
                                start_x = start_x,
                                n_washout = n_washout,
                                W = W, W_out = W_out, b = b, leaking_rate = leaking_rate,
                                n_run = n_run, 
                                which_patterns = which_patterns,
                                verbose = verbose)
  
  return(list('sg_outputs'=sg$sg_outputs))
}

#######################
### FULL SIMULATION ###
#######################
# store the patterns in a reservoir
storePatternsInReservoir <- function(W_in_raw, b_raw, W_star_raw, 
                                     W_in_scaling=1, b_scaling=1, W_star_scaling=1,
                                     leaking_rate=1,
                                     apertures=20,
                                     reg_out=1e-3, reg_W=1e-3,
                                     n_washout_max=100,
                                     verbose=F,
                                     show_plots=T){
  # scale weights
  scaled_weights <- scaleWeights(W_in_raw = W_in_raw, 
                                 W_in_scaling = W_in_scaling,
                                 b_raw = b_raw, 
                                 b_scaling = b_scaling,
                                 W_star_raw = W_star_raw, 
                                 W_star_scaling = W_star_scaling)
  W_in <- scaled_weights$W_in
  b <- scaled_weights$b
  W_star <- scaled_weights$W_star
  
  # drive reservoir for max n_washout_max time steps, but stop when two initial state vectors converge
  dr_washout <- driveReservoirWashoutPeriod(patterns = patterns,
                                            W_in = W_in, b = b, W_star = W_star, 
                                            leaking_rate = leaking_rate,
                                            n_washout_max = n_washout_max,
                                            verbose = verbose)
  n_washout <- dr_washout$n_washout
  start_x <- dr_washout$final_x_washout
  
  # drive reservoir for the n_learn time step and collect states
  n_learn <- n_points - n_washout
  dr_learn <- driveReservoirLearnPeriod(patterns = patterns,
                                        W_in = W_in, b = b, W_star = W_star, 
                                        leaking_rate = leaking_rate,
                                        start_x = start_x,
                                        n_washout = n_washout, n_learn = n_learn,
                                        verbose = verbose)
  all_training_x <- dr_learn$all_training_x
  all_training_x_old <- dr_learn$all_training_x_old
  all_training_output <- dr_learn$all_training_output
  all_training_W_target <- dr_learn$all_training_W_target
  target_outputs <- lapply(dr_learn$target_outputs, as.matrix)
  R <- dr_learn$R
  
  # plot some neuron states
  if(show_plots) plotNeuronStates(X=all_training_x, n=20, l=50)
  
  # compute the output weights
  cp_W_out <- computeOutputWeights(all_training_states = all_training_x, 
                                   all_training_output = all_training_output,
                                   reg_out = reg_out, verbose = verbose)
  W_out <- cp_W_out$W_out
  
  # recompute reservoir weights (loading the patterns into the reservoir)
  cp_W <- recomputeReservoirWeights(all_training_states_old = all_training_x_old, 
                                    all_training_W_target = all_training_W_target,
                                    reg_W = reg_W, verbose = verbose)
  W <- cp_W$W
  
  # compute conceptors
  Cs <- computeConceptors(R, apertures)
  
  # let reservoir self-generate with diagonal conceptors
  n_run <- n_learn
  sg <- selfGeneratingReservoir(Cs = Cs,
                                n_washout = n_washout,
                                W = W, W_out = W_out, b = b, 
                                leaking_rate =  leaking_rate,
                                n_run = n_run, verbose = verbose)
  sg_outputs <- sg$sg_outputs
  
  # shift the output of the sg_outputs for comparison
  shifted_patterns <- lapply(1:n_pattern, function(p) shiftPattern(sg_outputs[[p]], 
                                                                   target_outputs[[p]], 
                                                                   shift_along = 'rows'))
  
  # compare target outputs and self-generated outputs (plotting and Normalized Root Mean Square Error (NRMSE))
  par(mfrow=c(2,2), mar=rep(3,4))
  n_plot <- 40
  nrmses <- rep(NA, n_pattern)
  for(p in 1:n_pattern){
    # plot outputs
    if(show_plots) plot(shifted_patterns[[p]]$target_pattern[1:n_plot], type = 'l', xlab = '', ylab = '', ylim = c(-1,1))
    if(show_plots) lines(shifted_patterns[[p]]$pattern[1:n_plot], col=2, lwd=3, lty=2)
    
    nrmses[p] <- shifted_patterns[[p]]$nrmse
    if(verbose) print(paste0('pattern ',p,' nrmse : ', round(nrmses[p],5)))
    if(show_plots) legend(1, 1, legend=c(paste0("NRMSE: ",round(nrmses[p],5))))
  }
  
  reservoir <- list('patterns'= patterns,
                    'Cs'=Cs,
                    'W'=W,
                    'b'=b,
                    'W_out'=W_out,
                    'n_washout'=n_washout,
                    'leaking_rate'=leaking_rate)
  
  return(list('nrmses'=nrmses, 
              'reservoir'=reservoir))
}

