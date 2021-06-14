### CLEAR GLOBAL ENVIRONMENT ###
rm(list=ls(all.names = TRUE))

### REPRODUCIBILITY ###
set.seed(6)

### SOURCES ###
source('source/Libraries.R')
source('source/GeneralFunctions.R')

###########################
### EXPERIMENT SETTINGS ###
###########################
# 1 = periodic patterns sine - sine
# 2 = periodic patterns 5-periodic - 5-periodic 
# 3 = periodic patterns sine - 5-periodic
# 4 = periodic patterns all
# 5 = chaotic attractors
# 6 = human motion

experiment <- 6
show_progress_bar <- T

#########################
### PERIODIC PATTERNS ###
#########################
if(experiment %in% c(1,2,3)){
  print(paste0('running experiment ', experiment, ', please wait...'))
  
  # load reservoir
  load('output/reservoirs/reservoir_c_pp.RData')
  
  # set variables
  patterns <- reservoir$patterns
  Cs <- reservoir$Cs
  W <- reservoir$W
  b <- reservoir$b
  W_out <- reservoir$W_out
  leaking_rate <- reservoir$leaking_rate
  
  # system dimensions
  N <- nrow(W)
  M <- nrow(W_out)
  
  # morphing patterns
  if(experiment == 1){
    pattern_i <- 1
    pattern_j <- 2
  }else if(experiment == 2){
    pattern_i <- 3
    pattern_j <- 4
  }else{
    pattern_i <- 3
    pattern_j <- 1
  }
  
  # mu
  if(experiment == 1){
    mu_min <- -2
    mu_max <- 3
  }else if(experiment ==2){
    mu_min <- 0
    mu_max <- 1
  }else{
    mu_min <- 0
    mu_max <- 1
  }
  
  # period lengths
  if(experiment == 1){
    n_washout <- 100
    n_before_morph <- 500
    n_morph <- 200
    n_after_morph <- 500
  }else if(experiment ==2){
    n_washout <- 100
    n_before_morph <- 500
    n_morph <- 50
    n_after_morph <- 500
  }else{
    n_washout <- 100
    n_before_morph <- 500
    n_morph <- 50
    n_after_morph <- 500
  }
  n_run <- n_before_morph + n_morph + n_after_morph
  
  # collectors
  output_collector <- matrix(0,n_run,M)
  mu_collector <- rep(NA,n_run)
  
  # initial state
  r <- matrix(runif(N),N,1)
  z <- Cs[[pattern_i]] %*% r
  
  # create progress bar
  if(show_progress_bar) pb <- txtProgressBar(1,n_run+n_washout, style = 3)
  for(n in 1:(n_run+n_washout)){
    # set mu
    if(n < n_washout + n_before_morph){
      mu <- mu_min
    }else if(n >= n_washout + n_before_morph && n <= n_washout + n_before_morph + n_morph){
      mu <- mu_min + ((n - n_washout - n_before_morph)/n_morph) * (mu_max - mu_min)
      # if(n == n_washout + n_before_morph) print(paste0('begin mu=',mu))
      # if(n == n_washout + n_before_morph + n_morph) print(paste0('end mu=',mu))
    }else{
      mu <- mu_max
    }
    
    # update equation
    r <- tanh(W %*% z + b)
    C_morph <- (1-mu)*Cs[[pattern_i]] + mu*Cs[[pattern_j]]
    z <- (1-leaking_rate) * z + leaking_rate * C_morph %*% r
    
    # collect output
    if(n > n_washout){
      output_collector[n-n_washout,] <- W_out %*% r
      mu_collector[n-n_washout] <- mu
    }
    
    # update progress bar
    if(show_progress_bar) setTxtProgressBar(pb, n)
  }
  
  # plot the patterns
  if(experiment == 1){
    n_plot <- 300
  }else if(experiment ==2){
    n_plot <- 100
  }else{
    n_plot <- 100
  }
  plot_start_point <- n_before_morph - ceiling(0.5*(n_plot - n_morph))
  
  if(experiment == 1){
    main_title <- 'Frequency morphing'
  }else if(experiment ==2){
    main_title <- 'Shape morphing'
  }else{
    main_title <- 'Heterogeneous morphing'
  }
  
  par(mfrow=c(2,1), mar=c(3.1, 3, 2.1, 2.1))
  plot(output_collector[plot_start_point:(plot_start_point+n_plot)], 
       type = 'l', lwd=2,
       xlab = '', ylab = '', yaxt='n',
       main = main_title)
  axis(side = 2, at = c(-1,0,1))
  plot(mu_collector[plot_start_point:(plot_start_point+n_plot)], 
       col=2, lty=2, lwd=2,
       xlab = '', ylab = '', yaxt='n',
       main = paste0(mu_min,' \u2264 μ \u2264 ',mu_max))
  axis(side = 2, at = c(-2,-1,0,1,2,3))
  
  if(experiment == 1){
    # find periods
    sp_before <- spectrum(head(output_collector,n_before_morph), method = 'ar', plot=F)  
    period_before <- 1/sp_before$freq[which.max(sp_before$spec)]
    
    sp_after <- spectrum(tail(output_collector,n_after_morph), method = 'ar', plot=F)  
    period_after <- 1/sp_after$freq[which.max(sp_after$spec)]
    
    print(paste0('period mu=',mu_min,' -> ', round(period_before,3)))
    print(paste0('period mu=',mu_max,' -> ', round(period_after,3))) 
  }
}
if(experiment == 4){
  print(paste0('running experiment ', experiment, ', please wait...'))
  
  # load reservoir
  load('output/reservoirs/reservoir_c_pp.RData')
  
  # set variables
  patterns <- reservoir$patterns
  Cs <- reservoir$Cs
  W <- reservoir$W
  b <- reservoir$b
  W_out <- reservoir$W_out
  leaking_rate <- reservoir$leaking_rate
  
  # system dimensions
  N <- nrow(W)
  M <- nrow(W_out)
  
  # period lengths
  n_washout <- 100
  n_before_morph <- 500
  n_morph <- c(200,50,50)
  n_after_morph <- 500
  
  # outputs collector
  all_outputs <- list()
  for(sim in 1:3){
    # morphing patterns
    if(sim == 1){
      pattern_i <- 1
      pattern_j <- 2
    }else if(sim == 2){
      pattern_i <- 3
      pattern_j <- 4
    }else{
      pattern_i <- 3
      pattern_j <- 1
    }
    
    # mu
    if(sim == 1){
      mu_min <- -2
      mu_max <- 3
    }else if(sim ==2){
      mu_min <- 0
      mu_max <- 1
    }else{
      mu_min <- 0
      mu_max <- 1
    }
    
    n_run <- n_before_morph + n_morph[sim] + n_after_morph
    
    # collectors
    output_collector <- matrix(0,n_run,M)
    mu_collector <- rep(NA,n_run)
    
    # initial state
    r <- matrix(runif(N),N,1)
    z <- Cs[[pattern_i]] %*% r
    
    # create progress bar
    if(show_progress_bar) pb <- txtProgressBar(1,n_run+n_washout, style = 3)
    for(n in 1:(n_run+n_washout)){
      # set mu
      if(n < n_washout + n_before_morph){
        mu <- mu_min
      }else if(n >= n_washout + n_before_morph && n <= n_washout + n_before_morph + n_morph[sim]){
        mu <- mu_min + ((n - n_washout - n_before_morph)/n_morph[sim]) * (mu_max - mu_min)
        # if(n == n_washout + n_before_morph) print(paste0('begin mu=',mu))
        # if(n == n_washout + n_before_morph + n_morph[sim]) print(paste0('end mu=',mu))
      }else{
        mu <- mu_max
      }
      
      # update equation
      r <- tanh(W %*% z + b)
      C_morph <- (1-mu)*Cs[[pattern_i]] + mu*Cs[[pattern_j]]
      z <- (1-leaking_rate) * z + leaking_rate * C_morph %*% r
      
      # collect output
      if(n > n_washout){
        output_collector[n-n_washout,] <- W_out %*% r
        mu_collector[n-n_washout] <- mu
      }
      
      # update progress bar
      if(show_progress_bar) setTxtProgressBar(pb, n)
    }
    all_outputs[[sim]] <- output_collector
  }
  
  # plot the patterns
  n_plots <- c(300,100,100)
  par(mfrow=c(3,1), mar=c(3.1, 3, 2.5, 2.1), oma=c(0,0,3,0))
  for(sim in 1:length(n_plots)){
    plot_start_point <- n_before_morph - ceiling(0.5*(n_plots[sim] - n_morph[sim]))
    
    if(sim == 1){
      main_title <- 'Frequency morphing'
    }else if(sim ==2){
      main_title <- 'Shape morphing'
    }else{
      main_title <- 'Heterogeneous morphing'
    }
    
    plot(all_outputs[[sim]][plot_start_point:(plot_start_point+n_plots[sim])], 
         type = 'l', lwd=2, ylim = if(sim==2) c(-1.3,1.3),
         xlab = '', ylab = '', yaxt='n',
         main = main_title, cex.main=1.5, cex.axis=1.5)
    axis(side = 2, at = c(-1,0,1), cex.axis=1.5)
  }
  title('Conceptors', outer = T, cex.main=2.5)
}

##########################
### CHAOTIC ATTRACTORS ###
##########################
if(experiment %in% c(5)){
  print(paste0('running experiment ', experiment, ', please wait...'))
  
  # load reservoir
  load('output/reservoirs/reservoir_c_ca.RData')
  
  # set variables
  patterns <- reservoir$patterns
  Cs <- reservoir$Cs
  W <- reservoir$W
  b <- reservoir$b
  W_out <- reservoir$W_out
  leaking_rate <- reservoir$leaking_rate
  
  # system dimensions
  N <- nrow(W)
  M <- nrow(W_out)
  
  # patterns
  if(experiment == 5){
    pattern_i <- 1
    pattern_j <- 3
  }
  
  # mu
  if(experiment == 5){
    mu_min <- 0
    mu_max <- 1
  }
  
  # period lengths
  if(experiment == 5){
    n_washout <- 100
    n_before_morph <- 500
    n_morph <- 500
    n_after_morph <- 500
  }
  n_run <- n_before_morph + n_morph + n_after_morph
  
  # output collector
  output_collector <- matrix(0,n_run,M)
  mu_collector <- rep(NA,n_run)
  
  # initial state
  r <- matrix(runif(N),N,1)
  z <- Cs[[pattern_i]] %*% r
  
  # create progress bar
  if(show_progress_bar) pb <- txtProgressBar(1,n_run+n_washout, style = 3)
  for(n in 1:(n_run+n_washout)){
    # set mu
    if(n < n_washout + n_before_morph){
      mu <- mu_min
    }else if(n >= n_washout + n_before_morph && n <= n_washout + n_before_morph + n_morph){
      mu <- mu_min + ((n - n_washout - n_before_morph)/n_morph) * (mu_max - mu_min)
      # if(n == n_washout + n_before_morph) print(paste0('begin mu=',mu))
      # if(n == n_washout + n_before_morph + n_morph) print(paste0('end mu=',mu))
    }else{
      mu <- mu_max
    }
    
    # update equation
    r <- tanh(W %*% z + b)
    C_morph <- (1-mu)*Cs[[pattern_i]] + mu*Cs[[pattern_j]]
    z <- (1-leaking_rate) * z + leaking_rate * C_morph %*% r
    
    # collect output
    if(n > n_washout){
      output_collector[n-n_washout,] <- W_out %*% r
      mu_collector[n-n_washout] <- mu
    }
    
    # update progress bar
    if(show_progress_bar) setTxtProgressBar(pb, n)
  }
  
  # plot a few outputs to see how the figure changes over time
  if(experiment == 5){
    n_outputs <- 9
    n_plot <- 200
  }
  max_plot <- nrow(output_collector)-n_plot
  interval <- max_plot/(n_outputs-1)
  begin_end_points <- t(sapply(1:n_outputs, function(i) round(c(1+(i-1)*interval,(i-1)*interval+n_plot))))
  
  if(experiment == 5){
    main_title <- 'Rössler to Mackey-Glass'
  }
  
  # plot the patterns
  colfunc<-colorRampPalette(c("royalblue","springgreen"))
  par(mfrow=c(ceiling(sqrt(n_outputs)),round(sqrt(n_outputs))), 
      mar=rep(1,4), oma=c(0,0,4,0))
  for(i in 1:n_outputs){
    n <- ceiling(begin_end_points[i,1] + 0.5*(begin_end_points[i,2] - begin_end_points[i,1]))
    if(n < n_before_morph){
      col = colfunc(3)[1]
    }else if(n >= n_before_morph && n < n_before_morph + n_morph){
      col = colfunc(3)[2]
    }else{
      col = colfunc(3)[3]
    }
    
    plot(output_collector[begin_end_points[i,1]:begin_end_points[i,2],1], 
         output_collector[begin_end_points[i,1]:begin_end_points[i,2],2], 
         type = 'l', lwd=2, col=col,
         xlab = '', ylab = '', xaxt='n', yaxt='n')
  }
  # title(main = main_title, outer = T, cex.main=1.5, line = 0)
  title(main = 'Conceptors', outer = T, cex.main=3)
}

####################
### HUMAN MOTION ###
####################
if(experiment %in% c(6)){
  print(paste0('running experiment ', experiment, ', please wait...'))
  
  # load reservoir
  load('output/reservoirs/reservoir_c_hm_H.RData')
  
  # set variables
  patterns <- reservoir$patterns
  Cs <- reservoir$Cs
  W <- reservoir$W
  b <- reservoir$b
  W_out <- reservoir$W_out
  start_z <- reservoir$start_z
  n_learn <- reservoir$n_learn
  leaking_rate <- reservoir$leaking_rate
  
  n_pattern <- length(patterns) # number of patterns
  n_points <- sapply(patterns, FUN = function(x) nrow(x)) # length of each pattern
  
  # system dimensions
  N <- nrow(W)
  M <- nrow(W_out)
  
  # patterns
  if(experiment == 6){
    pattern_i <- 1
    pattern_j <- 3
  }
  
  # mu
  if(experiment == 6){
    mu_min <- 0
    mu_max <- 1
  }
  
  # period lengths
  if(experiment == 6){
    n_before_morph <- n_learn[pattern_i]
    n_morph <- 50
    n_after_morph <- n_learn[pattern_j]+200 # +100 is for the required washout period
  }
  n_run <- n_before_morph + n_morph + n_after_morph
  
  # output collector
  output_collector <- matrix(0,n_run,M)
  mu_collector <- rep(NA,n_run)
  
  # initial state
  z <- start_z[[pattern_i]]
  
  # create progress bar
  pb <- txtProgressBar(1,n_run, style = 3)
  for(n in 1:(n_run)){
    # set mu
    if(n < n_before_morph){
      mu <- mu_min
    }else if(n >= n_before_morph && n <= n_before_morph + n_morph){
      mu <- mu_min + ((n - n_before_morph)/n_morph) * (mu_max - mu_min)
      # if(n == n_before_morph) print(paste0('begin mu=',mu))
      # if(n == n_before_morph + n_morph) print(paste0('end mu=',mu))
    }else{
      mu <- mu_max
    }
    
    # update equation
    r <- tanh(W %*% z + b)
    C_morph <- (1-mu)*Cs[[pattern_i]] + mu*Cs[[pattern_j]]
    z <- (1-leaking_rate) * z + leaking_rate * C_morph %*% r
    
    # collect output
    output_collector[n,] <- W_out %*% r
    mu_collector[n] <- mu
    
    # update progress bar
    setTxtProgressBar(pb, n)
  }
  
  # plot a few outputs to see how the figure changes over time
  if(experiment == 6){
    which_dims <- c(3,4,5)
  }
  
  if(experiment == 6){
    main_title <- ''
  }
  
  # align 
  target <- patterns[[pattern_j]][(n_points[pattern_j]-n_learn[pattern_j]+1):n_points[pattern_j],]
  observed <- output_collector[n_learn[pattern_i]:nrow(output_collector),]
  pa <- phaseAlignPatterns(observed, target)
  
  ####################
  ### PLOT RESULTS ###
  ####################
  colfunc<-colorRampPalette(c("royalblue","springgreen"))
  par(mfrow=c(length(which_dims),3), mar=c(1,1,3,1), oma=c(0,0,3,0))
  for(i in 1:length(which_dims)){
    # plot the target pattern i in the left column
    plot(patterns[[pattern_i]][,which_dims[i]], type = 'l', lwd=2, 
         xaxt='n', yaxt = 'n', col = colfunc(3)[1],
         main = if(i==1) 'Boxing', cex.main=1.7)
    
    # plot morphed pattern in the middle column
    plot(output_collector[,which_dims[i]], type = 'l', 
         xlab = '', ylab = '', xaxt='n', yaxt = 'n', 
         col = colfunc(3)[2], lwd=3,
         main = if(i==1) 'Morph between boxing and cart wheel', cex.main=1.7)
    lines(tail(patterns[[pattern_i]][,which_dims[i]],n_learn[pattern_i]), 
          lwd=3, lty=2, col = colfunc(3)[1])
    lines((n_learn[pattern_i]+pa$shift+1):(n_learn[pattern_i]+pa$shift+n_learn[pattern_j]),
          tail(patterns[[pattern_j]][,which_dims[i]],n_learn[pattern_j]), 
          lwd=3, lty=2, col = colfunc(3)[3])
    
    # plot the target pattern j in the right column
    plot(patterns[[pattern_j]][,which_dims[i]], type = 'l', 
         lwd=3, col = colfunc(3)[3], 
         xaxt='n', yaxt = 'n',
         main = if(i==1) 'Cart wheel', cex.main=1.7)
  }
  title(main = 'Conceptors', outer = T, cex.main=3)
}
