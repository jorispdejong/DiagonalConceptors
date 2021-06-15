#############################
### Four Example Patterns ###
#############################
createExamplePatterns <- function(generate_new_data=F, L=1500,
                                  folder_path='data/periodic patterns/'){
  suppressPackageStartupMessages(library('numbers')) # mod() function
  
  paths <- paste0(folder_path, c("sin1.txt", 
                                 "sin2.txt",
                                 "5periodic1.txt",
                                 "5periodic2.txt"))
  
  if(generate_new_data){
    ### sin with period 2pi/8.83
    sin1 <- matrix(sin(2*pi*(1:L)/8.8342522),L,1)
    write.table(sin1, file = paths[1], row.names=FALSE, col.names=FALSE)
    
    ### sin with period 2pi/9.93
    sin2 <- matrix(sin(2*pi*(1:L)/9.8342522),L,1)
    write.table(sin2, file = paths[2], row.names=FALSE, col.names=FALSE)
    
    ### 5-periodic pattern random points
    random_points <- runif(5)
    random_points <- 1.8*(random_points-min(random_points))/(max(random_points)-min(random_points))-0.9;
    five_periodic1 <- matrix(random_points[mod((1:L)-1,5)+1],L,1)
    write.table(five_periodic1, file = paths[3], row.names=FALSE, col.names=FALSE)
    
    ### 5-periodic pattern random points + pertubation
    pertubation <- rnorm(5)
    random_pointsPert <- random_points + 0.2 * pertubation
    random_pointsPert <- 1.8*(random_pointsPert-min(random_pointsPert))/(max(random_pointsPert)-min(random_pointsPert))-0.9;
    five_periodic2 <- matrix(random_pointsPert[mod((1:L)-1,5)+1],L,1)
    write.table(five_periodic2, file = paths[4], row.names=FALSE, col.names=FALSE)
  }

  sin1 <- as.matrix(read.table(paths[1]))
  sin2 <- as.matrix(read.table(paths[2]))
  five_periodic1 <- as.matrix(read.table(paths[3]))
  five_periodic2 <- as.matrix(read.table(paths[4]))
  patterns <- list('sin 8.83'=sin1, 'sin 9.83'=sin2, 
                   '5-periodic'=five_periodic1, '5-periodic+pertubation'=five_periodic2)
  
  return(patterns)
}

##########################
### CHAOTIC ATTRACTORS ###
##########################
# returns the 4 chaotic attractors (Roessler, Lorenz, Mackey-Glass, Henon)
getChaoticAttractorsPatterns <- function(generate_new_data=F, L=1500,
                                         plot_result = T,
                                         folder_path='data/chaotic attractors/'){
  paths <- paste0(folder_path, c("RoesslerSeries.txt", 
                                 "LorenzSeries.txt",
                                 "MGSeries.txt",
                                 "HenonSeries.txt"))
  
  if(generate_new_data){
    ### Roessler
    Roessler_series <- generateRoesslerSeries(n_sample = L, plot_result = plot_result)
    write.table(Roessler_series, file = paths[1], row.names=FALSE, col.names=FALSE)
    
    ### Lorenz
    Lorenz_series <- generateLorenzSeries(n_sample = L, plot_result = plot_result)
    write.table(Lorenz_series, file = paths[2], row.names=FALSE, col.names=FALSE)
    
    ### Mackey-Glass
    MG_series <- generateMGSeries(n_sample = L, plot_result = plot_result)
    write.table(MG_series, file = paths[3], row.names=FALSE, col.names=FALSE)
    
    ### Lorenz
    Henon_series <- generateHenonSeries(n_sample = L, plot_result = plot_result)
    write.table(Henon_series, file = paths[4], row.names=FALSE, col.names=FALSE)
  }
  
  Roessler_series <- as.matrix(read.table(paths[1]))
  Lorenz_series <- as.matrix(read.table(paths[2]))
  MG_series <- as.matrix(read.table(paths[3]))
  Henon_series <- as.matrix(read.table(paths[4]))
  patterns <- list('Roessler'=Roessler_series, 'Lorenz'=Lorenz_series, 'MG'=MG_series, 'Henon'=Henon_series)
  
  return(patterns)
}
# generate chaotic attractor series
generateRoesslerSeries <- function(increments_per_unit=200, subsample_rate=150, n_sample, init_n_washout=5000, plot_result=T){
  rs <- c(0.5943, -2.2038, 0.0260) + 0.01 * rnorm(3)
  a <- 0.2 
  b <- 0.2 
  c <- 8.0 
  delta <- 1 / increments_per_unit # length of discrete approximation update interval
  Roessler_series <- matrix(0,2,n_sample)
  for(n in 1:init_n_washout){
    rs <- rs + delta * c(-(rs[2] + rs[3]), rs[1] + a * rs[2], b + rs[1]*rs[3] - c * rs[3])
  }
  for(n in 1:n_sample){
    for(k in 1:subsample_rate){
      rs <- rs + delta * c(-(rs[2] + rs[3]), rs[1] + a * rs[2], b + rs[1]*rs[3] - c * rs[3])
    }
    Roessler_series[,n] <- rs[1:2]
  }
  # normalize
  max_value <- apply(Roessler_series,1,max)
  min_value <- apply(Roessler_series,1,min)
  Roessler_series <- solve(diag(max_value - min_value)) %*% (Roessler_series - matrix(rep(min_value,n_sample),2,n_sample));
  
  # plot
  if(plot_result){
    par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
    plot(Roessler_series[1,], Roessler_series[2,], type = 'l')
  }
  
  return(t(Roessler_series))
}
generateLorenzSeries <- function(increments_per_unit=200, subsample_rate=15, n_sample, init_n_washout=5000, plot_result=T){
  # initialize Lorenz  state with a little random component
  ls <- c(10.036677794959058, 9.98674414052542, 29.024692318601613) + 0.01 * rnorm(3);
  sigma <- 10.0 
  b <- 8.0/3 
  r <- 28.0 
  delta <- 1 / increments_per_unit # length of discrete approximation update interval
  Lorenz_series <- matrix(0,2,n_sample)
  for(n in 1:init_n_washout){
    ls <- ls + delta * c(sigma * (ls[2]-ls[1]), r * ls[1] - ls[2] - ls[1]*ls[3], ls[1] * ls[2] - b * ls[3])
  }
  for(n in 1:n_sample){
    for(k in 1:subsample_rate){
      ls <- ls + delta * c(sigma * (ls[2]-ls[1]), r * ls[1] - ls[2] - ls[1]*ls[3], ls[1] * ls[2] - b * ls[3])
    }
    Lorenz_series[,n] <- ls[1:2]
  }
  # normalize
  max_value <- apply(Lorenz_series,1,max)
  min_value <- apply(Lorenz_series,1,min)
  Lorenz_series <- solve(diag(max_value - min_value)) %*% (Lorenz_series - matrix(rep(min_value,n_sample),2,n_sample));
  
  # plot
  if(plot_result){
    par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
    plot(Lorenz_series[1,], Lorenz_series[2,], type = 'l')
  }
  
  return(t(Lorenz_series))
}
generateMGSeries <- function(tau=17, increments_per_unit=10, subsample_rate=3, n_sample, init_n_washout=5000, plot_result=T){
  suppressPackageStartupMessages(library('numbers')) # mod() function
  
  gen_n_history <- tau * increments_per_unit
  seed <- 1.2 * rep(1,gen_n_history) + 0.2 * (runif(gen_n_history)-0.5)
  old_value <- 1.2
  gen_history <- seed
  
  MG_series <- matrix(0,2,n_sample)
  step <- 0
  for(n in 1:(n_sample+init_n_washout)){
    for(i in 1:(increments_per_unit*subsample_rate)){
      step <- step + 1
      tau_value <- gen_history[mod(step,gen_n_history)+1]
      new_value <- old_value + (0.2 * tau_value/(1.0 + tau_value^10) - 0.1 * old_value)/increments_per_unit
      gen_history[mod(step,gen_n_history)+1] <- old_value
      old_value <- new_value
    }
    if(n > init_n_washout) MG_series[,n - init_n_washout] <- c(new_value, tau_value)
  }
  
  # normalize
  max_value <- apply(MG_series,1,max)
  min_value <- apply(MG_series,1,min)
  MG_series <- solve(diag(max_value - min_value)) %*% (MG_series - matrix(rep(min_value,n_sample),2,n_sample));
  
  # plot
  if(plot_result){
    par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
    plot(MG_series[1,], MG_series[2,], type = 'l')
  }
  
  return(t(MG_series))
}
generateHenonSeries <- function(n_sample=L, init_n_washout=1000, plot_result=T){
  # initialize Henon  state with a little random component
  hs <- c(1.2677, -0.0278) + 0.01 * rnorm(2)
  a <- 1.4
  b <- 0.3
  
  Henon_series <- matrix(0,2,n_sample)
  n <- 1
  for(n in 1:init_n_washout){
    hs <- c(hs[2] + 1 - a * (hs[1]^2), b * hs[1])
    n <- n+1
  }
  for(n in 1:n_sample){
    hs <- c(hs[2] + 1 - a * hs[1]^2, b * hs[1])
    Henon_series[,n] <- hs
  }
  # normalize
  max_value <- apply(Henon_series,1,max)
  min_value <- apply(Henon_series,1,min)
  Henon_series <- solve(diag(max_value - min_value)) %*% (Henon_series - matrix(rep(min_value,n_sample),2,n_sample));
  
  # plot
  if(plot_result){
    par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
    plot(Henon_series[1,], Henon_series[2,])
  }
  return(t(Henon_series))
}

####################
### STICK PERSON ###
####################
# get raw data from files
getPatterns <- function(which_patterns=1:15, H=F){
  if(H){
    suppressPackageStartupMessages(library('R.matlab')) # readMat() function
    ### READ DATA ###
    import_path <- "data/human motion/raw/herbert/"
    box1_data <- readMat(paste0(import_path, "nnRawBox1.mat"))$nnRawDataBox1
    box2_data <- readMat(paste0(import_path, "nnRawBox2.mat"))$nnRawDataBox2
    box3_data <- readMat(paste0(import_path, "nnRawBox3.mat"))$nnRawDataBox3
    cart_wheel_data <- readMat(paste0(import_path, "nnRawCartWheel.mat"))$nnRawDataCartWheel
    crawl_data <- readMat(paste0(import_path, "nnRawCrawl.mat"))$nnRawDataCrawl
    exa_stride_data <- readMat(paste0(import_path, "nnRawExaStride.mat"))$nnRawDataExaStride
    get_down_data <- readMat(paste0(import_path, "nnRawGetDown.mat"))$nnRawDataGetdown
    get_seated_data <- readMat(paste0(import_path, "nnRawGetSeated.mat"))$nnRawDataGetSeated
    run_jog_data <- readMat(paste0(import_path, "nnRawRunJog.mat"))$nnRawDataRunJog
    sitting_data <- readMat(paste0(import_path, "nnRawSitting.mat"))$nnRawDataSitting
    slow_walk_data <- readMat(paste0(import_path, "nnRawSlowWalk.mat"))$nnRawDataSlowWalk
    standup_data <- readMat(paste0(import_path, "nnRawStandup.mat"))$nnRawDataStandup
    standup_from_stool_data <- readMat(paste0(import_path, "nnRawStandupFromStool.mat"))$nnRawDataStandupFromStool
    walk_data <- readMat(paste0(import_path, "nnRawWalk.mat"))$nnRawDataWalk
    waltz_data <- readMat(paste0(import_path, "nnRawWaltz.mat"))$nnRawDataWaltz
    
    ### PUT DATA IN LIST ###
    raw_patterns <- list('box1'=box1_data, 'box2'=box2_data,
                         'box3'=box3_data, 'cart_wheel'=cart_wheel_data,
                         'crawl'=crawl_data, 'exa_stride'=exa_stride_data,
                         'get_down'=get_down_data, 'get_seated'=get_seated_data,
                         'run_jog'=run_jog_data, 'sitting'=sitting_data,
                         'slow_walk'=slow_walk_data, 'standup'=standup_data,
                         'standup_from_stool'=standup_from_stool_data,
                         'walk'=walk_data, 'waltz'=waltz_data)
  }else{
    ### READ DATA ###
    import_path <- "data/human motion/raw/joris/"
    # subject 2
    bend_rise_data <- read.table(paste0(import_path, 'nnBendRise.txt'))
    punch_data <- read.table(paste0(import_path, 'nnPunch.txt'))
    jump_balance_data <- read.table(paste0(import_path, 'nnJumpBalance.txt'))
    walk_data <- read.table(paste0(import_path, 'nnWalk.txt'))
    # subject 54
    bear_data <- read.table(paste0(import_path, 'nnBear.txt'))
    chicken_data <- read.table(paste0(import_path, 'nnChicken.txt'))
    dragon_data <- read.table(paste0(import_path, 'nnDragon.txt'))
    ghost_data <- read.table(paste0(import_path, 'nnGhost.txt'))
    humming_bird_data <- read.table(paste0(import_path, 'nnHummingBird.txt'))
    monkey_data <- read.table(paste0(import_path, 'nnMonkey.txt'))
    penguin_data <- read.table(paste0(import_path, 'nnPenguin.txt'))
    praying_mantis_data <- read.table(paste0(import_path, 'nnPrayingMantis.txt'))
    roadrunner_data <- read.table(paste0(import_path, 'nnRoadrunner.txt'))
    squirrel_data <- read.table(paste0(import_path, 'nnSquirrel.txt'))
    superhero_data <- read.table(paste0(import_path, 'nnSuperhero.txt'))
    
    # put raw data into list
    raw_patterns <- list("bend_rise" = bend_rise_data, "punch" = punch_data, 
                         "jump_balance" = jump_balance_data, "walk" = walk_data, 
                         "bear" = bear_data, "chicken" = chicken_data, 
                         "dragon" = dragon_data, "ghost" = ghost_data, 
                         "humming_bird" = humming_bird_data, "monkey"=monkey_data, 
                         "penguin" = penguin_data, "praying_mantis" = praying_mantis_data,
                         "roadrunner" = roadrunner_data, "squirrel" = squirrel_data, 
                         "superhero"= superhero_data) 
  }
  return(raw_patterns[which_patterns])
}
# preprocess raw data (remove noise dimensions and normalize)
preprocessPatterns <- function(raw_patterns, save_data=T, folder_path='data/human motion/preprocessed/', H=T, smoothen=T){
  # H=TRUE means that we are using Herbert's data

  ### REMOVE DIMENSIONS THAT ARE CONSTANT OR ONLY NOISE
  # get patterns from folder
  column_indices <- if(H) c(5,17) else c(25, 26, 37, 38, 34, 46)
  removed_dimensions_patterns <- lapply(raw_patterns, FUN = function(x) x[,-column_indices])
  removed_dimensions <- lapply(raw_patterns, FUN = function(x) x[,column_indices])
  
  ### SMOOTHEN THE DATA USING SPLINES
  smooth_splines <- function(x, spar=0.5) smooth.spline(x,spar=spar)$y
  smooth_patterns <- lapply(removed_dimensions_patterns, FUN = function(x) apply(x, 2, FUN = smooth_splines))
  
  ### NORMALIZE THE DATA (COLUMN_WISE) 
  scalings <- getScalings(smooth_patterns) # get scaling factors of column of each pattern
  patterns <- normalizePatterns(smooth_patterns, scalings) # normalize each column of every pattern
  
  # set names
  names(patterns) <- names(raw_patterns)
  
  if(save_data) save(patterns, file = paste0(folder_path, ifelse(H,'herbert/','joris/'),'patterns.RData'))
  
  return(
    list(
      'patterns'=patterns, 
      'scalings'=scalings, 
      'removed_dimensions'=list(
        'column_indices'=column_indices, 
        'columns'=removed_dimensions)
    )
  )
}
