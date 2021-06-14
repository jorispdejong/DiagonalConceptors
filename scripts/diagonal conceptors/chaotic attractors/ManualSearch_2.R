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
file_path <- 'data/chaotic attractors/patterns.RData'
if(file.exists(file_path)){
  load(file_path)
}else{
  patterns <- getChaoticAttractorsPatterns(generate_new_data = T, L=1500)
  save(patterns, file=file_path) 
}
# set pattern variables
n_pattern <- length(patterns) # number of patterns
n_points <- sapply(patterns, FUN = function(x) dim(x)[1]) # length of each pattern

# plot the patterns
par(mfrow=c(2,2), mar=c(2.1, 3, 2.1, 2.1))
for(p in 1:4) plot(patterns[[p]][,1], patterns[[p]][,2], type = ifelse(p==4, 'p', 'l'), 
                   xlim = c(0,1), ylim = c(0,1), 
                   xlab = '', ylab = '')

##################
### PARAMETERS ###
##################
# scaling
W_in_scaling <- 1.1
b_scaling <- 0.8
W_star_scaling <- 1.2

# leaking rate
leaking_rate <- 1

# apertures
apertures <- c(30,30,10,10)

# regularization
reg_out <- 0.1
reg_W <- 0.015

# other parameters
n_washout_max <- 100
n_adapt_max <- 400
show_plots <- T
verbose <- T
save_results <- T

#######################
### RESERVOIR SETUP ###
#######################
# dimensions
M <- ncol(patterns[[1]]) # dimension of output 
N <- 400 # reservoir size

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
                                          n_washout_max = n_washout_max, tol = 1e-5,
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
                              n_washout = n_washout_max,
                              W = W, W_out = W_out, b = b, leaking_rate =  leaking_rate,
                              n_run = n_run, verbose = verbose)
sg_outputs <- sg$sg_outputs

####################
### PLOT RESULTS ###
####################
# plot the original patterns and the trained output
par(mfrow=c(2,4))
for(p in 1:n_pattern){
  # plot outputs
  if(show_plots) plot(target_outputs[[p]][,1],target_outputs[[p]][,2], type = ifelse(p==4, 'p', 'l'), 
                      xlim = c(0,1), ylim = c(0,1), xlab = '', ylab = '')
  if(show_plots) plot(sg_outputs[[p]][,1],sg_outputs[[p]][,2], type = ifelse(p==4, 'p', 'l'), col=2,
                      xlim = c(0,1), ylim = c(0,1), xlab = '', ylab = '')
}

####################
### SAVE RESULTS ###
####################
if(save_results){
  if(verbose) print('Saving results...')
  reservoir <- list('patterns'=patterns, 
                    'cs'=cs,
                    'W'=W,
                    'b'=b,
                    'W_out'=W_out,
                    'n_washout'=n_washout_max,
                    'leaking_rate'=leaking_rate)
  save(reservoir, file = 'output/reservoirs/reservoir_dc_ca_2.RData')
  if(verbose) print('Done!')
}