# clear global environment
rm(list=ls(all.names = TRUE))

# source
source("scripts/data processing/MoCapPackage.R")

# global variables
export_path <- "data/stick person/raw/joris/" 

#################
### SUBJECT 2 ###
#################
# import and export relative paths to folders
import_path_s2 <- "data/stick person/mocap data/Subject2_VariousExpressions/" 
# frame rate give in the README.txt file
frame_rate_s2 <- as.numeric(strsplit(readLines(paste0(import_path_s2, 'README.txt')), " ")[[1]][3])
# asf file of subject 2
asf_s2 <- readASF(paste0(import_path_s2,'2.asf'))

### bend rise
amc_bend_rise <- readAMCFromFile(paste0(import_path_s2, 'bend_rise.txt'), asf_s2)
raw_data_bend_rise <- t(amc_bend_rise$D)[1:350,]
write.table(raw_data_bend_rise[1:350,], paste0(export_path, 'nnBendRise.txt'), 
            row.names = F, col.names = F)

### punch
amc_punch <- readAMCFromFile(paste0(import_path_s2, 'punch.txt'), asf_s2)
raw_data_punch <- t(amc_punch$D)
write.table(raw_data_punch[1:430,], paste0(export_path, 'nnPunch.txt'), 
            row.names = F, col.names = F)

### jump balance
amc_jump_balance <- readAMCFromFile(paste0(import_path_s2, 'jump_balance.txt'), asf_s2)
raw_data_jump_balance <- t(amc_jump_balance$D)
write.table(raw_data_jump_balance, paste0(export_path, 'nnJumpBalance.txt'), 
            row.names = F, col.names = F)

### walk
amc_walk <- readAMCFromFile(paste0(import_path_s2, 'walk.txt'), asf_s2)
raw_data_walk <- t(amc_walk$D)
write.table(raw_data_walk, paste0(export_path, 'nnWalk.txt'), 
            row.names = F, col.names = F)

#################
### SUBJECT 54 ##
#################
# import and export relative paths to folders
import_path_s54 <- "data/stick person/mocap data/Subject54_AnimalBehaviors/" 
# frame rate give in the README.txt file
frame_rate_s54 <- as.numeric(strsplit(readLines(paste0(import_path_s54, 'README.txt')), " ")[[1]][3])
# asf file of subject 54
asf_s54 <- readASF(paste0(import_path_s54,'54.asf'))

### bear
amc_bear <- readAMCFromFile(paste0(import_path_s54, 'bear.txt'), asf_s54)
raw_data_bear <- t(amc_bear$D)
write.table(raw_data_bear[800:1350,], paste0(export_path, 'nnBear.txt'), 
            row.names = F, col.names = F)

### chicken
amc_chicken <- readAMCFromFile(paste0(import_path_s54, 'chicken.txt'), asf_s54)
raw_data_chicken <- t(amc_chicken$D)
write.table(raw_data_chicken[1:700,], paste0(export_path, 'nnChicken.txt'), 
            row.names = F, col.names = F)

### dragon 
amc_dragon <- readAMCFromFile(paste0(import_path_s54, 'dragon.txt'), asf_s54)
raw_data_dragon <- t(amc_dragon$D)
write.table(raw_data_dragon[150:470,], paste0(export_path, 'nnDragon.txt'), 
            row.names = F, col.names = F)

### ghost
amc_ghost <- readAMCFromFile(paste0(import_path_s54, 'ghost.txt'), asf_s54)
raw_data_ghost <- t(amc_ghost$D)
write.table(raw_data_ghost[1:300,], paste0(export_path, 'nnGhost.txt'), 
            row.names = F, col.names = F)

### humming bird
amc_humming_bird <- readAMCFromFile(paste0(import_path_s54, 'hummingbird.txt'), asf_s54)
raw_data_humming_bird <- t(amc_humming_bird$D)
write.table(raw_data_humming_bird[120:500,], paste0(export_path, 'nnHummingBird.txt'), 
            row.names = F, col.names = F)

### monkey
amc_monkey <- readAMCFromFile(paste0(import_path_s54, 'monkey.txt'), asf_s54)
raw_data_monkey <- t(amc_monkey$D)
write.table(raw_data_monkey[900:1150,], paste0(export_path, 'nnMonkey.txt'), 
            row.names = F, col.names = F)


### penguin 
amc_penguin <- readAMCFromFile(paste0(import_path_s54, 'penguin.txt'), asf_s54)
raw_data_penguin <- t(amc_penguin$D)
write.table(raw_data_penguin[1100:1580,], paste0(export_path, 'nnPenguin.txt'), 
            row.names = F, col.names = F)


### praying mantis
amc_praying_mantis <- readAMCFromFile(paste0(import_path_s54, 'praying_mantis.txt'), asf_s54)
raw_data_praying_mantis <- t(amc_praying_mantis$D)
write.table(raw_data_praying_mantis[1270:1480,], paste0(export_path, 'nnPrayingMantis.txt'), 
            row.names = F, col.names = F)

### roadrunner
amc_roadrunner <- readAMCFromFile(paste0(import_path_s54, 'roadrunner.txt'), asf_s54)
raw_data_roadrunner <- t(amc_roadrunner$D)
write.table(raw_data_roadrunner[200:700,], paste0(export_path, 'nnRoadrunner.txt'), 
            row.names = F, col.names = F)

### squirrel
amc_squirrel <- readAMCFromFile(paste0(import_path_s54, 'squirrel.txt'), asf_s54)
raw_data_squirrel <- t(amc_squirrel$D)
write.table(raw_data_squirrel[250:700,], paste0(export_path, 'nnSquirrel.txt'), 
            row.names = F, col.names = F)


### superhero 
amc_superhero <- readAMCFromFile(paste0(import_path_s54, 'superhero.txt'), asf_s54)
raw_data_superhero <- t(amc_superhero$D)
write.table(raw_data_superhero[1:480,], paste0(export_path, 'nnSuperhero.txt'), 
            row.names = F, col.names = F)



