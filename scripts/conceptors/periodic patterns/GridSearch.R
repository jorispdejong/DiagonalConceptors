### CLEAR GLOBAL ENVIRONMENT ###
rm(list=ls(all.names = TRUE))

### REPRODUCIBILITY ###
set.seed(6)

### SOURCES ###
source('source/Libraries.R')
source('source/GeneralFunctions.R')
source('source/ConceptorsFunctions.R')
source('source/PreprocessingFunctions.R')

##################
### DATA SETUP ###
##################
# read pre-generated patterns to avoid randomness
file_path <- 'data/periodic patterns/conceptors/patterns.RData'
if(file.exists(file_path)){
  load(file_path)
}else{
  patterns <- createExamplePatterns(generate_new_data=T, L=1000, 
                                    folder_path='data/periodic patterns/conceptors/')
  save(patterns, file=file_path) 
}
# set pattern variables
n_pattern <- length(patterns) # number of patterns
n_points <- sapply(patterns, FUN = function(x) dim(x)[1]) # length of each pattern

# plot (part of) the patterns
par(mfrow=c(2,2))
n_plot <- 40
for(p in 1:4) plot(patterns[[p]][1:n_plot], type = 'l', ylim = c(-1,1), ylab = '')

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

#######################################
### STORE PATTERNS IN THE RESERVOIR ###
#######################################
# other parameters
n_washout_max <- 100
show_plots <- F
verbose <- F
save_results <- T

# parameter sequences
W_in_scaling_s <- 0.8
b_scaling_s <- 0.5
W_star_scaling_s <- 1
leaking_rate_s <- 1
apertures_s <- seq(15,50,5)
reg_out_s <- (1e-1)^(3:6)
reg_W_s <- (1e-1)^(3:6)

# all possible combinations
parameter_combinations <- expand.grid('W_in_scaling'=W_in_scaling_s, 
                                      'b_scaling'=b_scaling_s, 
                                      'W_star_scaling'=W_star_scaling_s, 
                                      'leaking_rate'=leaking_rate_s, 
                                      'apertures'=apertures_s, 
                                      'reg_out'=reg_out_s, 
                                      'reg_W'=reg_W_s)

# Go through all the parameter combinations and save the reservoir that has the lowest error
# create progress bar
pb <- txtProgressBar(1,nrow(parameter_combinations), style=3)
best_i <- c()
err <- 1e6 # starting error
for(i in 1:nrow(parameter_combinations)){
  # store patterns in reservoir
  spir <- storePatternsInReservoir(W_in_raw = W_in_raw, 
                                   b_raw = b_raw, 
                                   W_star_raw = W_star_raw,
                                   W_in_scaling = parameter_combinations$W_in_scaling[i], 
                                   b_scaling = parameter_combinations$b_scaling[i], 
                                   W_star_scaling = parameter_combinations$W_star_scaling[i],
                                   leaking_rate = parameter_combinations$leaking_rate[i], 
                                   apertures = parameter_combinations$apertures[i],
                                   reg_out = parameter_combinations$reg_out[i],
                                   reg_W = parameter_combinations$reg_W[i],
                                   n_washout_max = n_washout_max,
                                   show_plots=show_plots,
                                   verbose = verbose)
  
  # compute overall nrmse over all patterns
  nrmse_norm <- sqrt(sum((spir$nrmses)^2))
  
  # if overall nrmse is lower than previous error then save the reservoir
  if(nrmse_norm < err){
    reservoir <- spir$reservoir
    if(save_results) save(reservoir, file='output/reservoirs/reservoir_c_pp.RData')
    err <- nrmse_norm
    best_i <- c(best_i, i)
    print(paste0('nrmse=',round(err,5)))
  }
  
  # update progress bar
  setTxtProgressBar(pb,i)
}





