### CLEAR GLOBAL ENVIRONMENT ###
rm(list=ls(all.names = TRUE))

### REPRODUCIBILITY ###
set.seed(6)

### LIBRARIES ###
suppressPackageStartupMessages(library('phaseR'))

##########################
### STABILITY ANALYSIS ###
##########################
# set the parameters
aperture <- 3
Q <- 0.5
parameters <- c(aperture, Q)

# define the derivative function
derivative_c <- function(t, c, parameters){
  a <- parameters[1]
  Q <- parameters[2]
  dc <- -Q * c^3 + Q * c^2 - a^-2 * c
  list(dc)
}

# compute the solutions
if(1 - 4 * aperture^-2 * Q^-1 < 0){
  solutions <- c(0)  
}else{
  c_minus <- 0.5 - 0.5 * sqrt(1 - 4 * aperture^-2 * Q^-1)
  c_plus <- 0.5 + 0.5 * sqrt(1 - 4 * aperture^-2 * Q^-1)
  c_zero <- 0
  solutions <- c(c_plus, c_minus, c_zero)
}

# plot the direction field
par(mfrow=c(1,1), mar=c(5,5,2.5,3))
derivative_c_flowField <- flowField(derivative_c,
                                    xlim = c(0, 60), ylim = c(0,1),
                                    parameters = parameters, 
                                    system = "one.dim", 
                                    add = F,
                                    arrow.head = 0.1, frac = 80,
                                    points = 20,
                                    col = 1,
                                    xlab = '', ylab = '',
                                    yaxt='n',
                                    cex.lab=1.5, cex.axis=1, cex.main=1.5, cex.sub=1.5,
                                    main = expression('stability of c'['i']*', aperture=3, E[r'['i']*']=0.5'))
axis(side = 2, at = c(round(solutions,2),1))
mtext('t', side=1, line=2.5, cex=1.5)
mtext(expression('c'['i']), side=2, line=2.5, cex=1.5)

# plot the solutions
for(sol in solutions) abline(h=round(sol,2), col=2, lwd=2)


