### CLEAR GLOBAL ENVIRONMENT ###
rm(list=ls(all.names = TRUE))

### REPRODUCIBILITY ###
set.seed(6)

### SOURCES ###
source('source/Libraries.R')
source('source/GeneralFunctions.R')
source('source/ConceptorsFunctions.R')
source('source/PreprocessingFunctions.R')
source('source/PostprocessingFunctions.R')

### HERBERT'S DATA OR NOT ###
H <- F

#################
### LOAD DATA ###
#################
load(paste0('output/conceptors/stick person/',ifelse(H,'herbert','joris'),'/reservoir.RData'))

#####################
### SET VARIABLES ###
#####################
preprocessed_patterns <- reservoir$preprocessed_patterns
patterns <- preprocessed_patterns$patterns
scalings <- preprocessed_patterns$scalings
removed_dimensions <- preprocessed_patterns$removed_dimensions
Cs <- reservoir$Cs
W <- reservoir$W
b <- reservoir$b
W_out <- reservoir$W_out
start_x <- reservoir$start_z
leaking_rate <- reservoir$leaking_rate

# system dimensions
N <- nrow(W)
M <- nrow(W_out)

# number of patterns
n_pattern <- length(Cs)

###################################
### LET RESERVOIR SELF GENERATE ###
###################################
# let reservoir self-generate with diagonal conceptors
n_run <- sapply(patterns, nrow)
sg_outputs <- list()
for(p in 1:n_pattern){
  pattern <- patterns[[p]]
  print(paste('pattern',p))
  
  output_collector <- matrix(0,M,n_run[p])
  
  # initial state
  z <- start_x[[p]]
  
  # create progress bar
  pb <- txtProgressBar(1,n_run[p], style = 3)
  for(n in 1:n_run[p]){
    r <- tanh(W %*% z + b)
    z <- (1-leaking_rate) * z + leaking_rate * Cs[[p]] %*% r
    
    output_collector[,n] <- W_out %*% r
    
    # update progress bar
    setTxtProgressBar(pb, n)
  }
  
  sg_outputs[[p]] <- t(output_collector)
}
names(sg_outputs) <- names(patterns)

# shift the output of the sg_outputs for comparison
shifted_patterns <- lapply(1:n_pattern, function(p) shiftPattern(patterns[[p]],
                                                                 sg_outputs[[p]],
                                                                 max_shift = 100,
                                                                 shift_along = 'rows'))

# plot (part of the) patterns
k_pattern <- 1
nrow_plots <- 3
par(mfrow=c(nrow_plots,nrow_plots), mar=c(3,3,1.5,3))
for(i in 1:nrow_plots^2){
  plot(shifted_patterns[[k_pattern]]$target_pattern[,i], type = 'l', ylim = c(-1,1), ylab = '')
  lines(shifted_patterns[[k_pattern]]$pattern[,i], col=2, lwd=2)
}

# compare target outputs and self-generated outputs (Normalized Root Mean Square Error (NRMSE))
nrmses <- lapply(shifted_patterns, function(x) x$nrmse)
for(p in 1:n_pattern){
  if(p==1) print('Normalized Root Mean Square Error')
  print(paste0('p',p,
               ': ', round(nrmses[[p]]$mean_nrmse,3), 
               ', min=', round(nrmses[[p]]$min_nrmse,3),
               ', max=',round(nrmses[[p]]$max_nrmse,3)))
}

#######################
### POST PROCESSING ###
#######################
post_process <- T
if(post_process){
  postProcessData(patterns = sg_outputs, 
                  scalings = scalings, 
                  removed_dimensions = removed_dimensions, 
                  output_path = paste0('output/conceptors/stick person/',
                                       ifelse(H,'herbert','joris'),
                                       '/json/'))
}
