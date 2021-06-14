#########################
### ERROR COMPUTATION ###
#########################
# Normalized Root-Mean-Square Error (NRMSE)
nrmse <- function(obs, pred){
  if(is.matrix(obs) && is.matrix(pred)){
    max_nrmse <- 0
    min_nrmse <- 1e6
    all_nrmses <- rep(NA, ncol(obs))
    for(i in 1:ncol(obs)){
      all_nrmses[i] <- sqrt(mean((obs[,i]-pred[,i])^2)/mean((obs[,i]-mean(obs[,i]))^2))
      if(all_nrmses[i] > max_nrmse) max_nrmse <- all_nrmses[i]
      if(all_nrmses[i] < min_nrmse) min_nrmse <- all_nrmses[i]
    }
    mean_nrmse <- sqrt(mean((obs-pred)^2)/mean((obs-mean(obs))^2))
    return(list('mean_nrmse'=mean_nrmse, 
                'min_nrmse'=min_nrmse,
                'max_nrmse'=max_nrmse, 
                'all_nrmses'=all_nrmses))
  }else if(is.list(obs) && is.list(pred)){
    nrmse_list <- list()
    for(j in 1:length(obs)){
      max_nrmse <- 0
      min_nrmse <- 1e6
      all_nrmses <- rep(NA, ncol(obs[[j]]))
      for(i in 1:ncol(obs[[j]])){
        all_nrmses[i] <- sqrt(mean((obs[[j]][,i]-pred[[j]][,i])^2)/mean((obs[[j]][,i]-mean(obs[[j]][,i]))^2))
        if(all_nrmses[i] > max_nrmse) max_nrmse <- all_nrmses[i]
        if(all_nrmses[i] < min_nrmse) min_nrmse <- all_nrmses[i]
      }
      mean_nrmse <- sqrt(mean((obs[[j]]-pred[[j]])^2)/mean((obs[[j]]-mean(obs[[j]]))^2))
      nrmse_list[[j]] <- list('mean_nrmse'=mean_nrmse, 
                              'min_nrmse'=min_nrmse,
                              'max_nrmse'=max_nrmse, 
                              'all_nrmses'=all_nrmses)
    }
    return(nrmse_list)
  }else{
    return(sqrt(mean((obs-pred)^2)/mean((obs-mean(obs))^2)))
  }
}
# Root-Mean-Square Error (RMSE)
rmse <- function(obs, pred){
  if(is.list(obs)){
    rmse <- 0
    for(i in 1:length(obs)){
      rmse <- rmse + (1/length(obs))*sqrt(mean((obs[[i]]-pred[[i]])^2))
    }
    return(rmse)
  }else{
    return(sqrt(mean((obs-pred)^2)))
  }
}

########################
### RIDGE REGRESSION ###
########################
# solving b=((X'X)^-1)X'y
ridgeRegression <- function(X,Y,reg){
  output <- t((solve(X %*% t(X) + diag(reg,dim(X)[1])) %*% X) %*% t(Y))
  nrmse <- nrmse(output %*% X, Y)$mean_nrmse
  return(list(output,nrmse))
}

#############################
### CONCEPTOR COMPUTATION ###
#############################
# computer conceptor matrix
computeConceptor <- function(R, aperture){
  R %*% solve(R + diag(aperture^(-2),N))
}
# compute all conceptor matrices in a list
computeConceptors <- function(R, apertures, verbose=F){
  if(!is.list(R)) stop('R must be a list')
  
  if(length(apertures)==1){
    apertures <- rep(apertures, length(R))
  }else{
    if(length(apertures)!=length(R)) stop('length of R and apertures must be equal')
  }
  
  C <- list()
  if(verbose) pb <- txtProgressBar(1,length(R),style = 3)
  for(i in 1:length(R)){
    C[[i]] <- computeConceptor(R[[i]], apertures[i])
    
    if(verbose) setTxtProgressBar(pb,i)
  }
  
  return(C)
}

#####################################################
### MATRIX COMPUTATION, MANIPULATION AND CREATION ###
#####################################################
# compute trace of a matrix
tr <- function(A) sum(diag(A))
# create normally distributed sparse matrix
createSparseMatrix <- function(nrow, ncol, density, dist = 'norm', range=c(-0.5,0.5)){
  n_items <- ceiling(nrow * ncol * density)
  non_zero_items <- sample(1:(nrow*ncol), n_items)
  vec <- rep(0, nrow * ncol)
  if(dist == 'norm'){
    vec[non_zero_items] <- rnorm(n_items)
  }else if(dist == 'uniform'){
    vec[non_zero_items] <- runif(n_items, range[1], range[2])
  }else{
    stop('distribution must be norm or uniform')
  }
  m <- matrix(vec, nrow = nrow)
  return(m)
}
# compute spectral radius of matrix
spectralRadius <- function(M){
  return(abs(eigen(M,only.values=TRUE)$values[1]))
}
# split vector into chunks of length L[i] and put into list
splitMatrixToList <- function(M,lengths,splitAlong='cols'){
  if(splitAlong=='rows'){
    if(sum(lengths)!=dim(M)[1]) stop('matrix cannot be split evenly')
    
    cs <- cumsum(lengths)
    l <- list()
    l[[1]] <- M[1:cs[1],]
    for(i in 2:length(lengths)){
      begin <- cs[i-1] + 1
      end <- cs[i]
      l[[i]] <- M[begin:end,]
    }
  }else{
    if(sum(lengths)!=dim(M)[2]) stop('x cannot be split evenly')
    
    cs <- cumsum(lengths)
    l <- list()
    l[[1]] <- M[,1:cs[1]]
    for(i in 2:length(lengths)){
      begin <- cs[i-1] + 1
      end <- cs[i]
      l[[i]] <- M[,begin:end]
    }
  }
  return(l)
}
# rbind a matrix multiple times
multipleRbind <- function(m, n=0, target_length=0){
  if(n==0 && target_length==0) stop('either n or target_length must be bigger than 0')
  
  if(target_length > 0){
    almost_full <- do.call("rbind", replicate(floor(target_length/nrow(m)), m, simplify = FALSE))
    mr <- rbind(almost_full, m[1:(target_length-nrow(almost_full)),])
  }else{
    mr <- do.call("rbind", replicate(n, m, simplify = FALSE))
  }
  return(mr)
}

#######################
### RESERVOIR SETUP ###
#######################
# create initial reservoir weights
createInitialReservoirWeights <- function(spectral_radius=1){
  if(N<20) 
    density <- 1
  else 
    density <- 10/N
  
  # create sparse matrix with values samples from normal distribution
  W_raw <- createSparseMatrix(N,N,density)
  # scale such that spectral radius is 1
  W <- (W_raw/spectralRadius(W_raw))*spectral_radius
  return(W)
}
# create initial input weights
createInitialInputWeights <- function() matrix(rnorm(N*M),N,M)
# create initial bias
createInitialBias <- function() matrix(rnorm(N*M),N,1)
# scale the input weight, bias weights, reservoir weights
scaleWeights <- function(W_in_raw, W_in_scaling,  
                         b_raw, b_scaling,
                         W_star_raw, W_star_scaling){
  W_in <- W_in_raw * W_in_scaling
  b <- b_raw * b_scaling
  W_star <- W_star_raw * W_star_scaling
  return(list('W_in'=W_in, 'b'=b, 'W_star'=W_star))
}

###########################
### RECOMPUTING WEIGHTS ###
###########################
computeOutputWeights <- function(all_training_states, all_training_output, reg_out, verbose=F){
  ### compute W_out
  if(verbose) print('computing W_out...')
  W_out_training <- ridgeRegression(X=all_training_states, Y=all_training_output, reg = reg_out)
  if(verbose) print('done!')
  W_out <- W_out_training[[1]]
  if(verbose) print(paste0('mean NRMSE W_out: ', round(W_out_training[[2]],5)))
  return(list('W_out'=W_out))
}
recomputeReservoirWeights <- function(all_training_states_old, all_training_W_target, reg_W, verbose=F){
  ### recompute W
  if(verbose) print('recomputing W...')
  W_training <- ridgeRegression(X=all_training_states_old, Y=all_training_W_target, reg = reg_W)
  if(verbose) print('done!')
  W <- W_training[[1]]
  if(verbose) print(paste0('mean NRMSE W: ', round(W_training[[2]],5)))
  return(list('W'=W))
}

#####################
### NORMALIZATION ###
#####################
# get min and max of each column of a matrix
getScalings <- function(patterns){
  l <- lapply(patterns, FUN = function(y){
    minmax <- apply(y, 2, FUN = function(x) c(min(x),max(x)))
    rownames(minmax) <- c('min', 'max')
    return(minmax)})
  return(l)
}
# normalize between values a and b
normalize <- function(x,min,max,a=-1,b=1){
  n <- a + (b-a)*((x-min(x))/(max(x)-min(x)))
  return(n)
}
# loop over all patterns and normalize column-wise using the scaling factors 
normalizePatterns <- function(patterns, scalings){
  norm_patterns <- lapply(matrix(1:length(scalings), nrow=1), 
                          FUN = function(j)
                            apply(matrix(1:ncol(scalings[[j]]), nrow=1), 2, 
                                  FUN = function(i) normalize(patterns[[j]][,i],
                                                              min = scalings[[j]][1,i], 
                                                              max = scalings[[j]][2,i])))
  return(norm_patterns)
}

################
### PLOTTING ###
################
# plot random neuron states
plotNeuronStates <- function(X, n=8, l=50){
  par(mfrow=c(1,1))
  random_neurons <- sample(1:N, n)
  matplot(t(X[random_neurons,]), type='l', xlim = c(0,l), ylim = c(-1,1), ylab = '')
  title('Some neuron states')
}

#########################
### SHIFTING PATTERNS ###
#########################
# find the phase of two vectors
alignVectors <- function(vec1, vec2, chunk_size=50, keep_only_intersection=F){
  if(length(vec1) > length(vec2)){
    big_vec <- vec1
    small_vec <- vec2
  }else{
    big_vec <- vec2
    small_vec <- vec1
  }
  errors <- rep(NA, length(small_vec))
  chunks <- chunk(small_vec, chunk_size)
  
  smallest_error <- NA
  chunk_index <- NA
  best_index <- NA
  for(k in 1:length(chunks)){
    n_runs <- length(big_vec)-length(chunks[[k]])
    errors <- rep(NA, n_runs)
    for(i in 1:n_runs){
      errors[i] <- sqrt(sum( (big_vec[i:(length(chunks[[k]])+i-1)] - chunks[[k]])^2 ))
    }
    if(k==1){
      smallest_error <- min(errors)
      chunk_index <- k
      best_index <- which.min(errors)
    }else{
      if(min(errors) < smallest_error){
        smallest_error <- min(errors)
        chunk_index <- k
        best_index <- which.min(errors)
      }
    }
  }
  if(chunk_index==1){
    prepend_chunk <- rep(NA, best_index)
    append_chunk <- Reduce(c, chunks[(chunk_index+1):length(chunks)])
    if(length(big_vec) > best_index + length(chunks[[chunk_index]]) + length(append_chunk)){
      append_chunk <- c(append_chunk, rep(NA, length(big_vec) - (best_index + length(chunks[[chunk_index]]) + length(append_chunk))))
      output_big_vec <- big_vec
    }else{
      output_big_vec <- c(big_vec, rep(NA, best_index + length(chunks[[chunk_index]]) + length(append_chunk) - length(big_vec)))
    }
    output_small_vec <- c(prepend_chunk, chunks[[chunk_index]], append_chunk)
  }
  
  if(chunk_index==length(chunks)){
    prepend_chunk <- Reduce(c, chunks[1:(chunk_index-1)])
    append_chunk <- rep(NA, length(big_vec)-best_index-length(chunks[[chunk_index]]))
    if(length(prepend_chunk) < best_index){
      prepend_chunk <- c(rep(NA, best_index - length(prepend_chunk)), prepend_chunk)
      output_big_vec <- big_vec
    }else{
      output_big_vec <- c(rep(NA, length(prepend_chunk) - best_index), big_vec)
    }
    output_small_vec <- c(prepend_chunk, chunks[[chunk_index]], append_chunk)
  }
  
  if(chunk_index > 1 && chunk_index < length(chunks)){
    prepend_chunk <- Reduce(c, chunks[1:(chunk_index-1)])
    if(length(prepend_chunk) < best_index){
      prepend_chunk <- c(rep(NA, best_index - length(prepend_chunk)), prepend_chunk)
      output_big_vec <- big_vec
    }else{
      output_big_vec <- c(rep(NA, length(prepend_chunk) - best_index), big_vec)
    }
    
    append_chunk <- Reduce(c, chunks[(chunk_index+1):length(chunks)])
    if(length(big_vec) > best_index + length(chunks[[chunk_index]]) + length(append_chunk)){
      append_chunk <- c(append_chunk, rep(NA, length(big_vec) - (best_index + length(chunks[[chunk_index]]) + length(append_chunk))))
      output_big_vec <- output_big_vec
    }else{
      output_big_vec <- c(output_big_vec, rep(NA, best_index + length(chunks[[chunk_index]]) + length(append_chunk) - length(big_vec)))
    }
    
    output_small_vec <- c(prepend_chunk, chunks[[chunk_index]], append_chunk)
  }
  
  if(keep_only_intersection){
    intersection <- intersect(which(!is.na(output_small_vec)),which(!is.na(output_big_vec)))
    output_small_vec <- output_small_vec[intersection]
    output_big_vec <- output_big_vec[intersection] 
  }
  
  if(length(vec1) > length(vec2)){
    return(list('vec1'=output_big_vec, 'vec2'=output_small_vec))
  }
  
  return(list('vec1'=output_small_vec, 'vec2'=output_big_vec))
}
# shift one vector with a maximum number of steps to find the least error overlap
shiftPattern <- function(pattern, target_pattern, max_shift=10, shift_along='cols'){
  if(dim(pattern)[1]!=dim(target_pattern)[1] || dim(pattern)[2]!=dim(target_pattern)[2]) 
    stop('pattern and target_pattern must have equal dimensions') 
  
  error <- rep(NA, max_shift)
  if(shift_along=='cols'){
    for(i in 1:max_shift){
      error[i] <- sqrt(sum( (target_pattern[,1:(ncol(target_pattern)-i+1)] - pattern[,i:ncol(pattern)])^2 ))
    }
    
    best_i <- which.min(error)
    
    pattern_cropped <- pattern[,best_i:ncol(pattern)]
    target_pattern_cropped <- target_pattern[,1:(ncol(target_pattern)-best_i+1)]
  }else{
    for(i in 1:max_shift){
      error[i] <- sqrt(sum( (target_pattern[1:(nrow(target_pattern)-i+1),] - pattern[i:nrow(pattern),])^2 ))
    }
    
    best_i <- which.min(error)
    
    pattern_cropped <- pattern[best_i:nrow(pattern),]
    target_pattern_cropped <- target_pattern[1:(nrow(target_pattern)-best_i+1),]
  }
  
  nrmse <- nrmse(pattern_cropped, target_pattern_cropped)
  
  l <- list('pattern'=pattern_cropped, 'target_pattern'=target_pattern_cropped, 'nrmse'=nrmse)
  return(l)
}
# phase align patterns. Shift the smaller pattern over the big one and find where the error is the smallest
phaseAlignPatterns <- function(observed, target){
  if(is.list(observed) && is.list(target)){
    return_list <- list()
    pb <- txtProgressBar(1,length(observed),style = 3)
    for(p in 1:length(observed)){
      observed_p <- observed[[p]]
      target_p <- target[[p]]
      
      l_target_p <- nrow(target_p)
      l_observed_p <- nrow(observed_p)
      l <- l_observed_p - l_target_p
      
      if(l==0){
        return_list[[p]] <- list('shift'=best_i,
                                 'observed_pattern'=observed_p,
                                 'target_pattern'=target_p,
                                 'nrmse'=nrmse(observed_p, target_p))
      }else{
        for(i in 1:abs(l)){
          observed_i <- if(l>0) as.matrix(observed_p[i:(i+l_target_p-1),],nrow=l_target_p-1) else observed_p
          target_i <- if(l>0) target_p else as.matrix(target_p[i:(i+l_observed_p-1),],nrow=l_observed_p-1)
          nrmse_i <- nrmse(observed_i,target_i)
          mean_nrmse <- nrmse_i$mean_nrmse
          if(i==1){
            best_i <- i
            best_mean_nrmse <- mean_nrmse
            best_observed <- observed_i
            best_target <- target_i
            best_nrmse <- nrmse_i
          }else{
            if(mean_nrmse < best_mean_nrmse){
              best_mean_nrmse <- mean_nrmse
              best_i <- i
              best_nrmse <- nrmse_i
              best_observed <- observed_i
              best_target <- target_i
            } 
          }
        }
        
        return_list[[p]] <- list('shift'=best_i,
                                 'observed_pattern'=best_observed,
                                 'target_pattern'=best_target,
                                 'nrmse'=best_nrmse)
      }
      setTxtProgressBar(pb,p)
    }
    return(return_list)
  }
  
  if(is.matrix(observed) && is.matrix(target)){
    l_target <- nrow(target)
    l_observed <- nrow(observed)
    l <- l_observed - l_target
    
    if(l==0){
      return(list('observed_pattern'=observed,
                  'target_pattern'=target,
                  'nrmse'=nrmse(observed, target)))
    }else{
      for(i in 1:abs(l)){
        observed_i <- if(l>0) observed[i:(i+l_target-1),] else observed
        target_i <- if(l>0) target else target[i:(i+l_observed-1),]
        nrmse_i <- nrmse(observed_i,target_i)
        mean_nrmse <- nrmse_i$mean_nrmse
        if(i==1){
          best_i <- i
          best_mean_nrmse <- mean_nrmse
          best_observed <- observed_i
          best_target <- target_i
          best_nrmse <- nrmse_i
        }else{
          if(mean_nrmse < best_mean_nrmse){
            best_mean_nrmse <- mean_nrmse
            best_i <- i
            best_nrmse <- nrmse_i
            best_observed <- observed_i
            best_target <- target_i
          } 
        }
      }
      
      return(list('shift'=best_i,
                  'observed_pattern'=best_observed,
                  'target_pattern'=best_target,
                  'nrmse'=best_nrmse))
    }
  }
}

######################
### SORT XY VALUES ###
######################
# sort x and y values that lie on the boundary of some area
# so that a line can be drawn from one point to the adjacent point etc.
sortXY <- function(x,y){
  points_sorted <- c(1)
  x_old <- x[1]
  y_old <- y[1]
  x_sorted <- c(x_old)
  y_sorted <- c(y_old)
  count <- 0
  while(length(x_sorted) < length(x) || count > length(x)){
    remaining_points <- (1:length(x))[-points_sorted]
    next_point <- order(sqrt((x_old - x[remaining_points])^2 + (y_old - y[remaining_points])^2))[1]
    
    x_new <- x[remaining_points][next_point]
    y_new <- y[remaining_points][next_point]
    
    x_sorted <- c(x_sorted, x_new)
    y_sorted <- c(y_sorted, y_new)
    
    x_old <- x_new
    y_old <- y_new
    
    points_sorted <- c(points_sorted, which(x[remaining_points][next_point]==x))
    
    count <- count + 1
  }
  
  return(list(x=x_sorted,y=y_sorted))
}
