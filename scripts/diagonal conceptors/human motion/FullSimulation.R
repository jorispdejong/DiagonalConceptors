### CLEAR GLOBAL ENVIRONMENT ###
rm(list=ls(all.names = TRUE))

### REPRODUCIBILITY ###
set.seed(6)

### FOLDER PATH ###
folder_path1 <- 'Diagonal Conceptors/'
folder_path2 <- 'Stick Person/'
folder_path <- paste0(folder_path1, folder_path2)

### SOURCES ###
source('Source/Libraries.R')
source('Source/GeneralFunctions.R')
source('Source/DiagonalConceptorsFunctions.R')
source('Source/PreprocessingFunctions.R')

### HERBERT'S DATA OR NOT ###
H <- F

##################
### DATA SETUP ###
##################
# read pre-generated patterns to avoid randomness
if(H){
  file_path_preprocessed_patterns <- paste0(folder_path, 'Results/RData/preprocessed_patterns_H.RData')
}else{
  file_path_preprocessed_patterns <- paste0(folder_path, 'Results/RData/preprocessed_patterns.RData')
} 
if(file.exists(file_path_preprocessed_patterns)){
  load(file_path_preprocessed_patterns)
}else{
  # get patterns from folder
  raw_patterns <- getPatterns(H = H)
  
  # preprocess patterns and save the scalings parameters
  preprocessed_patterns <- preprocessPatterns(raw_patterns, H = H)
  
  # save preprocessed patterns
  save(preprocessed_patterns, file=file_path_preprocessed_patterns) 
}
# get patterns
removed_dimensions <- preprocessed_patterns$removed_dimensions
scalings <- preprocessed_patterns$scalings
patterns <- preprocessed_patterns$patterns

# set pattern variables
n_pattern <- length(patterns) # number of patterns
n_points <- sapply(patterns, FUN = function(x) dim(x)[1]) # length of each pattern

# plot (part of the) patterns
k_pattern <- 6
nrow_plots <- 3
par(mfrow=c(nrow_plots,nrow_plots))
for(i in 1:nrow_plots^2) plot(patterns[[k_pattern]][,i], type = 'l', ylim = c(-1,1), ylab = '')

##################
### PARAMETERS ###
##################
# scaling
W_in_scaling <- 0.25
b_scaling <- 0.2
W_star_scaling <- 0.8

# leaking rate
leaking_rate <- 0.4

# apertures
apertures <- 5

# regularization
reg_out <- 3.5
reg_W <- 3.5

# other parameters
n_washout_max <- 50
n_adapt_max <- 50
show_plots <- T
verbose <- T
save_results <- F

#######################
### RESERVOIR SETUP ###
#######################
# dimensions
M <- ncol(patterns[[1]]) # dimension of output 
N <- 700 # reservoir size

# create raw (not scaled yet) weights
W_in_raw <- createInitialInputWeights() # dim=NxM (never modified after initialization, only scaled)
b_raw <- createInitialBias() # dim=Nx1(never modified after initialization, only scaled)
W_star_raw <- createInitialReservoirWeights(spectral_radius=1) # spectral radius = 1

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

# set initial random diagonal conceptors
cs_init <- lapply(1:n_pattern, FUN=function(p) matrix(runif(N),N,1))

######################
### WASHOUT PERIOD ###
######################
# drive reservoir for max n_washout_max time steps, but stop when two initial state vectors converge
dr_washout <- driveReservoirWashoutPeriod(patterns = patterns,
                                          cs_init = cs_init,
                                          W_in = W_in, b = b, W_star = W_star, leaking_rate = leaking_rate,
                                          n_washout_max = n_washout_max, tol = 0,
                                          verbose = verbose)
n_washout <- dr_washout$n_washout
final_cz_washout <- dr_washout$final_cz_washout

#########################
### ADAPTATION PERIOD ###
#########################
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

####################
### LEARN PERIOD ###
####################
# drive reservoir for the n_learn time step and collect states
n_learn <- n_points - (n_washout+1) - (n_adapt+1)
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
if(show_plots) plotNeuronStates(X=all_training_z, n=20, l=150)

#############################
### (RE)COMPUTING WEIGHTS ###
#############################
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

##############################
### SELF-GENERATE PATTERNS ###
##############################
# let reservoir self-generate with diagonal conceptors
n_run <- n_learn
sg <- selfGeneratingReservoir(cs = cs,
                              start_cz = start_cz,
                              W = W, W_out = W_out, b = b, leaking_rate =  leaking_rate,
                              n_run = n_run, verbose = verbose)
sg_outputs <- sg$sg_outputs

####################
### PLOT RESULTS ###
####################
# plot (part of the) patterns
k_pattern <- 1
nrow_plots <- 3
par(mfrow=c(nrow_plots,nrow_plots), mar=c(3,3,1,1))
for(i in 1:nrow_plots^2){
  plot(target_outputs[[k_pattern]][,i], type = 'l', ylim = c(-1,1), ylab = '')
  lines(sg_outputs[[k_pattern]][,i], col=2, lty=2, lwd=2)
}
# compare target outputs and self-generated outputs (Normalized Root Mean Square Error (NRMSE))
nrmses <- rep(NA, n_pattern)
for(p in 1:n_pattern){
  nrmses[p] <- nrmse(target_outputs[[p]], sg_outputs[[p]])
  if(verbose) print(paste0('pattern ',p,' nrmse : ', round(nrmses[p],5)))
}

####################
### SAVE RESULTS ###
####################
if(save_results){
  reservoir <- list('Cs'=Cs,
                    'W'=W,
                    'b'=b,
                    'W_out'=W_out,
                    'start_cz'=start_x,
                    'n_run'=n_run,
                    'leaking_rate'=leaking_rate)
  if(H){
    file_name_reservoir <- paste0(folder_path,'Results/RData/reservoir_H.RData')
  }else{
    file_name_reservoir <- paste0(folder_path,'Results/RData/reservoir.RData')
  }
  save(reservoir, file = file_name_reservoir)
}

