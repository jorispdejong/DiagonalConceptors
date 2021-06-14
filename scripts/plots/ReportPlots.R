### CLEAR GLOBAL ENVIRONMENT ###
rm(list=ls(all.names = TRUE))

### REPRODUCIBILITY ###
set.seed(123)

##################################
### RESERVOIR CONNECTIONS PLOT ###
##################################
# create evenly spaced points on disc
n_points <- 20
indices <- 1:n_points
r <- sqrt(indices/n_points)
theta <- pi * (1 + sqrt(5)) * indices
x <- r*cos(theta)
y <- r*sin(theta)

# plot points
par(mfrow=c(1,1))
plot(x,y, pch=19, cex=3)
arrows(x[1],y[1])
# create and plot connections
dens <- 0.2



