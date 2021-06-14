# source
source("scripts/data processing/MoCapPackage.R")

######################
### SCALE PATTERNS ###
######################
# loop over all patterns and scale patterns back to original scale, column-wise 
denormalizePatterns <- function(patterns, scalings){
  denorm_patterns <- lapply(matrix(1:length(scalings), nrow=1), 
                            FUN = function(j){
                              apply(matrix(1:ncol(scalings[[j]]), nrow=1), 2, 
                                    FUN = function(i){
                                      normalize(patterns[[j]][,i], 
                                                min = -1, 
                                                max = 1,
                                                a = scalings[[j]][1,i],
                                                b = scalings[[j]][2,i])
                                    }
                              )
                            }
  )
  return(denorm_patterns)
}

################################
### INSERT COLUMNS IN MATRIX ###
################################
# insert columns into a matrix
insertColumns <- function(M,cols,positions){
  suppressPackageStartupMessages(library('tibble'))
  
  for(i in 1:length(positions)){
    M <- add_column(as.data.frame(M), d=cols[,i], .after = positions[i]-1)
  }
  
  return(as.matrix(M))
}

postProcessData <- function(patterns, scalings, removed_dimensions, output_path){
  ### SCALE BACK TO ORIGINAL SCALE
  # denormalize each column of every pattern
  patterns_s <- denormalizePatterns(sg_outputs, scalings)
  
  # function that insert columns into the p-th pattern
  insertCol <- function(p){
    cols <- matrix(0,nrow(patterns_s[[p]]),length(removed_dimensions$column_indices))
    positions <- removed_dimensions$column_indices
    return(insertColumns(patterns_s[[p]],cols,positions))
  }
  
  # self-generated patterns, scaled appropriately
  patterns_o <- lapply(1:n_pattern, insertCol)
  
  ### CONVERT DATA TO XYZ COORDINATES
  asf_s2 <- readASF("data/stick person/mocap data/Subject2_VariousExpressions/2.asf")
  asf_s54 <- readASF("data/stick person/mocap data/Subject54_AnimalBehaviors/54.asf")
  # pattern 1:4 are of subject 2
  # pattern 5:15 are of subject 54
  print('Post processing patterns... Please wait')
  for(p in 1:length(patterns_o)){
    print(paste0('Pattern: ', names(patterns)[p]))
    
    asf <- if(p %in% 1:4) asf_s2 else asf_s54
    
    D <- patterns_o[[p]]
    amc <- readAMC(matrix2amc(D), asf)
    xyz <- AMC2XYZ(asf, amc)
    XYZ2Json(xyz, asf, fileName = paste0(output_path,names(patterns)[p]), frameRate = 120)
  }
}
