# libraries
suppressPackageStartupMessages(library('jsonlite')) # convert data frame to json


######################
### PRE PROCESSING ###
######################

### convert asf and amc to xyz
AMC2XYZ <- function(asf, amc, verbose = FALSE) {
  xyz <- rep(list(matrix(0, nrow = amc$nFrames, ncol = 3)), nrow(amc$skeleton))
  names(xyz) <- amc$skeleton$bone
  bonesVector <- split(amc$skeleton[, c("x", "y", "z")], seq(nrow(amc$skeleton[, c("x", "y", "z")])))
  names(bonesVector) <- amc$skeleton$bone
  xyz <- traverseSkeleton("root", asf, amc, xyz, 1, amc$nFrames,bonesVector, verbose)
  return(xyz)
}

### read an asf file into an object
readASF <- function(asfFilePath, in2cm = TRUE) {
  lines <- readLines(asfFilePath)

  skeleton <- data.frame(bone = character(0),
                         dirx = numeric(0), diry = numeric(0), dirz = numeric(0),
                         length = numeric(0),
                         angx = numeric(0), angy = numeric(0), angz = numeric(0),
                         dofx = numeric(0), dofy = numeric(0), dofz = numeric(0),
                         stringsAsFactors = FALSE)
  
  # read in length
  idx <- 1
  while (all(strsplit(lines[idx], "\\s+")[[1]] != "length")) idx <- idx + 1
  lengthLoc <- which(strsplit(lines[idx], "\\s+")[[1]] == "length") + 1
  len <- as.numeric(strsplit(lines[idx], "\\s+")[[1]][lengthLoc])
  
  # read in angle type
  idx <- idx + 1
  while (all(strsplit(lines[idx], "\\s+")[[1]] != "angle")) idx <- idx + 1
  angleTypeLoc <- which(strsplit(lines[idx], "\\s+")[[1]] == "angle") + 1
  angleType <- strsplit(lines[idx], "\\s+")[[1]][angleTypeLoc]
  
  # read in root
  idx <- idx + 1
  while (all(strsplit(lines[idx], "\\s+")[[1]] != "position")) idx <- idx + 1
  positionLoc <- which(strsplit(lines[idx], "\\s+")[[1]] == "position") + 1
  rootPosition <- as.numeric(
    strsplit(lines[idx], "\\s+")[[1]][positionLoc:(positionLoc + 2)])
  
  idx <- idx + 1
  while (all(strsplit(lines[idx], "\\s+")[[1]] != "orientation")) idx <- idx + 1
  orientLoc <- which(strsplit(lines[idx], "\\s+")[[1]] == "orientation") + 1
  rootOrientation <- as.numeric(
    strsplit(lines[idx], "\\s+")[[1]][orientLoc:(orientLoc + 2)])
  
  skeleton[nrow(skeleton) + 1, ] <- c("root",
                                      rootPosition[1], rootPosition[2],
                                      rootPosition[3],
                                      0,
                                      rootOrientation[1], rootOrientation[2],
                                      rootOrientation[3],
                                      1, 1, 1)
  
  # read in bones
  while (lines[idx] != ":hierarchy") {
    resBone <- readBone(idx, lines)
    idx <- resBone$idx + 1
    skeleton[nrow(skeleton) + 1, ] <- resBone$bone
  }
  
  # read in hierarchy
  childs <- rep(list(NULL), nrow(skeleton))
  names(childs) <- skeleton$bone
  skeleton$child <- 1:nrow(skeleton)
  skeleton$parent <- rep(0, nrow(skeleton))
  idx <- idx + 2
  while(trimws(lines[idx]) != "end") {
    parentChildren <- strsplit(trimws(lines[idx]), " ")[[1]]
    bone <- parentChildren[1]
    for (child in 2:length(parentChildren)) {
      childs[[bone]] <- c(childs[[bone]], parentChildren[child])
      skeleton[skeleton$bone == parentChildren[child], "parent"] <-
        skeleton[skeleton$bone == bone, "child"]
    }
    idx <- idx + 1
  }
  
  # make sure skeleton valus are numeric
  skeleton[, 2:ncol(skeleton)] <-
    apply(skeleton[, 2:ncol(skeleton)], 2, as.numeric)
  
  # get bones final vector: length x direction
  skeleton$x <- skeleton$length * skeleton$dirx * ifelse(in2cm, 2.54, 1) / len
  skeleton$y <- skeleton$length * skeleton$diry * ifelse(in2cm, 2.54, 1) / len
  skeleton$z <- skeleton$length * skeleton$dirz * ifelse(in2cm, 2.54, 1) / len
  
  # get C and CInv matrices
  CMatList <- apply(skeleton[, c("angx", "angy", "angz")], 1,
                    function(x) getCAxisAnglrMatrices(x, angleType))
  names(CMatList) <- skeleton$bone
  
  return(list(skeleton = skeleton,
              childs = childs,
              CMatList = CMatList,
              len = len,
              angleType = angleType))
}
### read an amc file into an object
readAMCFromFile <- function(amcFilePath, asf) {
  lines <- readLines(amcFilePath)
  res <- amc2Matrix(lines)

  D <- res$D
  D[1:3, ] <- D[1:3, ] * 2.54 / asf$len
  
  bones <- c(res$bones, setdiff(asf$skeleton$bone, res$bones))
  
  skeleton <- asf$skeleton[match(bones, asf$skeleton$bone),]
  
  degrees <- rep(list(matrix(0, nrow = ncol(D), ncol = 3)), length(bones))
  names(degrees) <- bones
  for (frame in 1:ncol(D)) {
    counter <- 4
    for (b in 1:length(bones)) {
      for (colNum in 1:3) {
        if (skeleton[b, 8 + colNum] == 1) {
          degrees[[bones[b]]][frame, colNum] <- D[counter, frame]
          counter <- counter + 1
        }
      }
    }
  }
  return(list(D = D, degrees = degrees,
              nFrames = ncol(D),
              skeleton = skeleton))
}
### read amc lines into an object
readAMC <- function(lines, asf) {
  res <- amc2Matrix(lines)
  
  D <- res$D
  D[1:3, ] <- D[1:3, ] * 2.54 / asf$len

  bones <- c(res$bones, setdiff(asf$skeleton$bone, res$bones))
  
  skeleton <- asf$skeleton[match(bones, asf$skeleton$bone),]
  
  degrees <- rep(list(matrix(0, nrow = ncol(D), ncol = 3)), length(bones))
  names(degrees) <- bones
  for (frame in 1:ncol(D)) {
    counter <- 4
    for (b in 1:length(bones)) {
      for (colNum in 1:3) {
        if (skeleton[b, 8 + colNum] == 1) {
          degrees[[bones[b]]][frame, colNum] <- D[counter, frame]
          counter <- counter + 1
        }
      }
    }
  }
  return(list(D = D, degrees = degrees,
              nFrames = ncol(D),
              skeleton = skeleton))
}

#######################
### POST PROCESSING ###
#######################
### function that exports a json that can be read by Unity
XYZ2Json <- function(xyzList, asf, fileName, frameRate){
  exportJson <- character(0)
  childs <- asf$childs
  
  # create progress bar
  pb <- txtProgressBar(min = 1, max = nrow(xyzList[[1]]), style = 3)  
  for(frame in 1:nrow(xyzList[[1]])){
    pos <- as.data.frame(t(sapply(xyzList, function(x) x[frame, ]))[, c(1, 3, 2)])
    colnames(pos) <- c("x", "z", "y")
    
    segmentsDataFrame <- data.frame(x1 = numeric(0), y1 = numeric(0),
                                    z1 = numeric(0), x2 = numeric(0),
                                    y2 = numeric(0), z2 = numeric(0))
    count <- 1
    for (parent in names(childs)) {
      for (child in childs[[parent]]) {
        pParent <- as.matrix(pos[rownames(pos) == parent,])
        pChild <- as.matrix(pos[rownames(pos) == child, ])
        segmentsListEntry <- c(pParent, pChild)
        segmentsDataFrame[count, ] <- segmentsListEntry
        count <- count + 1
      }
    }  
    
    frameJson = sprintf('{"t":%s, "segments": %s}', frame/frameRate, jsonlite::toJSON(segmentsDataFrame))
    if(frame == 1){ 
      exportJson <- paste0("[",frameJson)
    }else if(frame == nrow(xyzList[[1]])){
      exportJson <- paste0(exportJson,",",frameJson, "]")
      exportJson <- sprintf('{"segmentsOverTime":%s}', exportJson)
    }else{
      exportJson <- paste0(exportJson,",",frameJson) 
    }
    # update progress bar
    setTxtProgressBar(pb, frame)
  }
  write(exportJson, paste0(fileName, ".json"))
}
### function that converts the D matrix to bone segments and exports to json file
export2Json <- function(D, fileDestination, fileName, subject){
  # set subject variables
  if(subject == "s2"){
    import_path <- "./Data Processing/Data/Raw/Subject2_VariousExpressions/"
    asf <- readASF(paste0(import_path,'2.asf'))
  }else if(subject == "s54"){
    import_path <- "./Data Processing/Data/Raw/Subject54_AnimalBehaviors/"
    asf <- readASF(paste0(import_path,'54.asf'))
  }else{
    stop("subject must be 's2' or 's54'")
  }
  frame_rate <- as.numeric(strsplit(readLines(paste0(import_path, 'README.txt')), " ")[[1]][3])
  
  # scale first three columns
  D[,1:3] <- D[,1:3] / 2.54 * asf$len
  
  # convert D matrix to amc lines
  print("converting D to amc")
  amc_lines <- matrix2amc(D)
  # convert amc lines to amc object
  amc <- readAMC(amc_lines, asf)
  # convert amc object to xyz coordinates
  cat("converting amc to xyz")
  xyz <- AMC2XYZ(asf, amc)
  # convert xyz coordinates to json and export to chosen destination
  cat("converting xyz to json and exporting json to file destination")
  XYZ2Json(xyz, asf, paste0(fileDestination, fileName), frame_rate)  
}

########################
### HELPER FUNCTIONS ###
########################

### function to read a single bone
readBone <- function(idx, lines) {
  idxBegin <- idx
  while (all(strsplit(lines[idxBegin], "\\s+")[[1]] != "begin")) {
    idxBegin <- idxBegin + 1
  }
  idxEnd <- idxBegin
  while (all(strsplit(lines[idxEnd], "\\s+")[[1]] != "end")) {
    idxEnd <- idxEnd + 1
  }
  boneString <- paste(lines[idxBegin:idxEnd], collapse = " ")
  boneName <- strsplit(
    regmatches(boneString,
               regexpr("name [a-z]+", boneString)), "\\s+")[[1]][2]
  boneDirections <- as.numeric(
    strsplit(regmatches(boneString,
                        regexpr("direction (-?[0-9.]+e?[-+]?[0-9]*[ ]?){3}",
                                boneString)), "\\s+")[[1]][2:4])
  boneLength <- as.numeric(
    strsplit(regmatches(boneString,
                        regexpr("length [0-9\\.]+",
                                boneString)), "\\s+")[[1]][2])
  boneAngles <- as.numeric(
    strsplit(regmatches(boneString,
                        regexpr("axis (-?[0-9.]+e?[-+]?[0-9]*[ ]?){3}",
                                boneString)), "\\s+")[[1]][2:4])
  boneDOFX <- ifelse(grepl("rx", boneString), 1, 0)
  boneDOFY <- ifelse(grepl("ry", boneString), 1, 0)
  boneDOFZ <- ifelse(grepl("rz", boneString), 1, 0)
  
  bone <- list(boneName,
               boneDirections[1], boneDirections[2], boneDirections[3],
               boneLength,
               boneAngles[1], boneAngles[2], boneAngles[3],
               boneDOFX, boneDOFY, boneDOFZ)
  return(list(idx = idxEnd, bone = bone))
}
### function to convert amc file to matrix
amc2Matrix <- function(lines, verbose = FALSE) {

  # read-in header
  idx <- 2
  while (lines[idx - 1] != ":DEGREES") idx <- idx + 1
  
  D <- numeric()
  dims <- c(6, 3, 3, 3, 3, 3, 3, 2, 3, 1, 1, 2, 1, 2, 2, 3, 1, 1, 2, 1, 2, 3,
            1, 2, 1, 3, 1, 2, 1)
  locations <- c(1, 7, 10, 13, 16, 19, 22, 25, 27, 30, 31, 32, 34, 35, 37, 39,
                 42, 43, 44, 46, 47, 49, 52, 53, 55, 56, 59, 60, 62)
  segments <- c('root', 'lowerback', 'upperback', 'thorax', 'lowerneck',
                'upperneck', 'head', 'rclavicle', 'rhumerus', 'rradius',
                'rwrist', 'rhand', 'rfingers', 'rthumb', 'lclavicle',
                'lhumerus', 'lradius', 'lwrist', 'lhand', 'lfingers', 'lthumb',
                'rfemur', 'rtibia', 'rfoot', 'rtoes', 'lfemur', 'ltibia',
                'lfoot', 'ltoes')
  
  getSegemtnID <- function(idx) strsplit(lines[idx], " ")[[1]][1]
  
  # read-in data
  # labels can be in any order
  frame <- 1
  while (idx <= length(lines)) {
    if (verbose & frame %% 100 == 0) cat("Reading frame: ", frame, "\n")
    
    row <- rep(0, 62)
    
    # read angle label
    idx <- idx + 1
    frameSegments <- sapply(idx:(idx + 28), getSegemtnID)
    index <- match(segments, frameSegments)
    
    # where to put the data
    location <- locations[index]
    len <- dims[index]
    
    for (i in 1:29) {
      row[location[i]:(location[i] + len[i] - 1)] <-
        as.numeric(strsplit(lines[idx], " ")[[1]][2:(len[i] + 1)])
      idx <- idx + 1
    }
    
    D <- cbind(D, row)
    
    frame <- frame + 1
  }
  if (verbose) cat("Total number of frames read: ", frame - 1, "\n")
  return(list(D = D, bones = segments))
}
### function to convert matrix to amc file
matrix2amc <- function(D){
  # D is a matrix with 62 columns and ... rows
  
  header <- c('#!Matlab matrix to amc conversion', ':FULLY-SPECIFIED',':DEGREES')
  
  dims <- c(6, 3, 3, 3, 3, 3, 3, 2, 3, 1, 1, 2, 1, 2, 2, 3, 1, 1, 2, 1, 2, 3,
            1, 2, 1, 3, 1, 2, 1)
  locations <- c(1, 7, 10, 13, 16, 19, 22, 25, 27, 30, 31, 32, 34, 35, 37, 39,
                 42, 43, 44, 46, 47, 49, 52, 53, 55, 56, 59, 60, 62)
  segments <- c('root', 'lowerback', 'upperback', 'thorax', 'lowerneck',
                'upperneck', 'head', 'rclavicle', 'rhumerus', 'rradius',
                'rwrist', 'rhand', 'rfingers', 'rthumb', 'lclavicle',
                'lhumerus', 'lradius', 'lwrist', 'lhand', 'lfingers', 'lthumb',
                'rfemur', 'rtibia', 'rfoot', 'rtoes', 'lfemur', 'ltibia',
                'lfoot', 'ltoes')
  
  lines <- c()
  for(frame in 1:nrow(D)){
    details <- rep(NA, length(segments)+1)
    details[1] <- sprintf('%s',frame);
    for(i in 1:length(segments)){
      details[i+1] <- paste(segments[i],paste(D[frame,locations[i]:(locations[i]+dims[i]-1)], collapse = " "))
    }
    lines <- c(lines, details)
  }
  amc <- c(header, lines)
  return(amc)
}
### function to convert degrees to radians
deg2rad <- function(deg) deg * pi / 180
### function to get the position and initial root matrix
getRootMatrix <- function(frame, D, angleType, CMatList) {
  root_pos <- D[1:3, frame]
  root_angle <- D[4:6, frame]
  root_matrix <- getLRotationMatrix("root",
                                    root_angle[1], root_angle[2], root_angle[3],
                                    angleType, CMatList)
  return(list(root_pos = root_pos, root_matrix = root_matrix))
}
### function to get a rotation matrix around a single axis
rotationMatrix <- function(angle, axis) {
  if (axis == "x") {
    matrix(c(1, 0, 0,
             0, cos(angle), -sin(angle),
             0, sin(angle), cos(angle)),
           nrow = 3, byrow = TRUE)
  }else if (axis == "y") {
    matrix(c(cos(angle), 0, sin(angle),
             0, 1, 0,
             -sin(angle), 0, cos(angle)),
           nrow = 3, byrow = TRUE)
  } else if (axis == "z") {
    matrix(c(cos(angle), -sin(angle), 0,
             sin(angle), cos(angle), 0,
             0, 0, 1),
           nrow = 3, byrow = TRUE)
  } else {
    stop("unknown axis parameter, only x, y, z allowed.")
  }
}
### function to get L rotation matrix
getLRotationMatrix <- function(bone, ax, ay, az, angleType, CMatList) {
  ax <- ifelse(angleType == "deg", deg2rad(ax), ax)
  ay <- ifelse(angleType == "deg", deg2rad(ay), ay)
  az <- ifelse(angleType == "deg", deg2rad(az), az)
  
  Mx <- t(rotationMatrix(ax, "x"))
  My <- t(rotationMatrix(ay, "y"))
  Mz <- t(rotationMatrix(az, "z"))
  
  M <- Mx %*% My %*% Mz
  L <- CMatList[[bone]]$Cinv %*% M %*% CMatList[[bone]]$C
  
  return(L)
}
### function to get initial C and Cinv Angle-Axis matrices
getCAxisAnglrMatrices <- function(C_values, angleType) {
  ax <- ifelse(angleType == "deg", deg2rad(C_values[1]), C_values[1])
  ay <- ifelse(angleType == "deg", deg2rad(C_values[2]), C_values[2])
  az <- ifelse(angleType == "deg", deg2rad(C_values[3]), C_values[3])
  
  Cx <- t(rotationMatrix(ax, "x"))
  Cy <- t(rotationMatrix(ay, "y"))
  Cz <- t(rotationMatrix(az, "z"))
  
  C <- Cx %*% Cy %*% Cz
  Cinv <- solve(C)
  
  return(list(C = C, Cinv = Cinv))
}
### function to get L rotation matrix at a frame
getBoneMatrix <- function(bone, frame, degrees, angleType, CMatList) {
  bone_rotation <- degrees[[bone]][frame, ]
  L <- getLRotationMatrix(bone,
                          bone_rotation[1], bone_rotation[2], bone_rotation[3],
                          angleType, CMatList)
  return(L)
}
### traverse a skeleton calling DFS for every frame
traverseSkeleton <- function(bone, asf, amc, xyz, nStartFrame, nEndFrame, bonesVector, verbose = FALSE) {
  stack <- list()
  for (frame in nStartFrame:nEndFrame) {
    if (verbose & frame %% 100 == 0) cat("frame: ", frame, "\n")
    xyz <- DFS("root", asf, amc, xyz, bonesVector, frame, stack = stack)
  }
  return(xyz)
}
### Depth First Search with a stack at a frame
DFS <- function(bone, asf, amc, xyz, bonesVector, frame, parent = NULL,stack = NULL) {
  if (bone == "root") {
    root_mat_res <- getRootMatrix(frame, amc$D, asf$angleType, asf$CMatList)
    pos <- root_mat_res$root_pos
    mat <- root_mat_res$root_matrix
    xyz[[bone]][frame, ] <- pos
    stack <- c(stack, list(mat))
  } else {
    stack <- c(stack, list(
      getBoneMatrix(bone, frame, amc$degrees, asf$angleType, asf$CMatList) %*%
        stack[[length(stack)]]))
    xyz[[bone]][frame, ] <- xyz[[parent]][frame, ] +
      as.numeric(bonesVector[[bone]]) %*%
      stack[[length(stack)]]
  }
  
  for (b in asf$childs[[bone]]) {
    xyz <- DFS(b, asf, amc, xyz, bonesVector, frame, parent = bone,
               stack = stack)
  }
  stack <- stack[-length(stack)]
  return(xyz)
}

