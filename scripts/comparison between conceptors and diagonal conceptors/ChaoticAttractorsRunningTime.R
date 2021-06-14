### CLEAR GLOBAL ENVIRONMENT ###
rm(list=ls(all.names = TRUE))

### REPRODUCIBILITY ###
set.seed(6)

### SOURCES ###
source('source/Libraries.R')
source('source/GeneralFunctions.R')

#############################
### EXPERIMENT PARAMETERS ###
#############################
n_runs <- c(1000,5000,10000)
running_times <- data.frame('conceptors'=rep(NA,length(n_runs)),
                            'diagonal conceptors'=rep(NA,length(n_runs)))

##################
### CONCEPTORS ###
##################
# load the reservoirs
load('output/conceptors/chaotic attractors/reservoir.RData')
Cs <- reservoir$Cs
W <- reservoir$W
b <- reservoir$b
W_out <- reservoir$W_out
n_washout <- reservoir$n_washout
leaking_rate <- reservoir$leaking_rate

# system dimensions
N <- nrow(W)
M <- nrow(W_out)

# number of patterns
n_pattern <- length(Cs)

# loop over number of runs
pb <- txtProgressBar(1,length(n_runs), style = 3) # create progress bar
for(j in 1:length(n_runs)){
  # set begin time
  begin_time <- Sys.time()
  
  for(p in 1:n_pattern){
    # initial state
    r <- matrix(runif(N),N,1)
    z <- Cs[[p]] %*% r
    
    for(n in 1:(n_runs[j]+n_washout[p])){
      r <- tanh(W %*% z + b)
      z <- (1-leaking_rate) * z + leaking_rate * Cs[[p]] %*% r
    }
  }
  
  # set end time
  end_time <- Sys.time()
  
  # compute and store running time
  diff_time <- end_time - begin_time
  running_times$conceptors[j] <- diff_time[[1]]
  
  # update progress bar
  setTxtProgressBar(pb,j)
}

###########################
### DIAGONAL CONCEPTORS ###
###########################
# load the reservoirs
load('output/diagonal conceptors/chaotic attractors/reservoir.RData')
cs <- reservoir$cs
W <- reservoir$W
b <- reservoir$b
W_out <- reservoir$W_out
n_washout <- reservoir$n_washout
leaking_rate <- reservoir$leaking_rate

# system dimensions
N <- nrow(W)
M <- nrow(W_out)

# number of patterns
n_pattern <- length(cs)

# loop over number of runs
pb <- txtProgressBar(1,length(n_runs), style = 3) # create progress bar
for(j in 1:length(n_runs)){
  # set begin time
  begin_time <- Sys.time()
  
  for(p in 1:n_pattern){
    # initial state
    r <- matrix(runif(N),N,1)
    z <- cs[[p]] * r
    
    for(n in 1:(n_runs[j]+n_washout[p])){
      r <- tanh(W %*% z + b)
      z <- (1-leaking_rate) * z + leaking_rate * cs[[p]]*r
    }
  }
  
  # set end time
  end_time <- Sys.time()
  
  # compute and store running time
  diff_time <- end_time - begin_time
  running_times$diagonal.conceptors[j] <- diff_time[[1]]
  
  # update progress bar
  setTxtProgressBar(pb,j)
}

####################
### PLOT RESULTS ###
####################
par(mfrow=c(1,1))
plot(n_runs, running_times$conceptors, type = 'l',
     ylim = c(0,max(running_times$conceptors)), 
     xlab = 'number of time steps', ylab = 'time in seconds',
     main = 'Conceptors vs. Diagonal Conceptors Running Time', lty=2, lwd=2)
points(n_runs, running_times$conceptors, cex=1.2, pch=16)
lines(n_runs, running_times$diagonal.conceptors, lty=1, lwd=2)
points(n_runs, running_times$diagonal.conceptors, cex=1.2, pch=16)
legend('topleft', legend = c('conceptors', 'diagonal conceptors'), 
       lwd=2, lty=c(2,1))


