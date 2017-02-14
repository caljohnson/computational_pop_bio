#assignment_4.R
#Carter Johnson
#ECL 233 Winter 2017

# Regardless of your choice, you'll need set up your ode function for multiple species:
# Initialize your population sizes with labels, e.g., n0=c(x=x0, y=y0, ...), then use these labels in your ode function
# odeFn = function(t, n, parms) {
# with(as.list(c(parms, n)), { # extract parameters from parms vector
# dx = ... # dx/dt
# dy = ... # dy/dt
# ...
# return(list(c(dx, dy, ...))) # return dn/dt as a list
# })
# }
# then running lsoda gets you each species as out$x, out$y, etc.

# Choice (a) competition
# fit paramecium data (from Leslie 1957 reading) to Lotka-Volterra competition model
# dn1/dt = r1*n1*(1-n1/K1-a12*n2/K1)
# dn2/dt = r2*n2*(1-n2/K2-a21*n1/K2)
# plot actual and simulated time series for both populations
# you get r's and K's - this is in the "parameciumrk.txt" file - don't forget to use colnames(your input matrix) = NULL to erase the default V1, V2, etc. column names
# you'll have to find the alpha's to the data in the "paramecium2.txt" file
# note that this means you have to reconstruct the parameters vector within optim - start with parmsrk=r and K, then do parms = c(parmsrk, a=a) to add a's within minimizing function

rm(list=ls()) # clear everything

#ODE RHS
competition = function(t, n, parms) {
  with(as.list(c(parms,n)), { # extract parameters from parms vector
    dn1 = r1*n1*(1-n1/K1 - a1*n2/K1)  #dn1/dt
    dn2 = r2*n2*(1-n2/K2 - a2*n1/K2)  # dn2/dt
    return(list(c(dn1,dn2))) # return dn/dt as a list
  })
}

#Optimize parameters in ODE
competitionMin = function(parms) {
  out = as.data.frame(lsoda(n0, times, competition, parms)) # run ode
  mse = (mean((out$n1-nTrue[,1])^2)+mean((out$n2-nTrue[,2])^2))/2 # calculate mean squared error between simulated and original data
  return(mse) # return mean squared error
}

# parameters and data
rK = as.matrix(read.table("parameciumrk.txt")) # actual data for r1,r2,K1,K2
colnames(rK) = NULL # take out default column names
parmsrk = c(r=c(rK[1,1],rK[1,2]), K=c(rK[2,1], rK[2,2])) #store r,K as labeled vector

nTrue = as.matrix(read.table("paramecium2.txt")) # actual data for a optimization
colnames(nTrue) = NULL # take out default column names
n0 = c(n=c(nTrue[1,1],nTrue[1,2])) # initial population sizes
tf = dim(nTrue)[1] # time
times = 1:tf # vector of times to run over

# need to start with guesses at the parameter values to give something for optim to start with
a1=0.5*rK[2,1]/rK[2,2]
a2=0.5*rK[2,2]/rK[2,1] #initial guess for a1, a2
parms0 = c(parmsrk,a1=a1, a2=a2) # set up parameter vector initial guess

# find fit and run simulations
optimOut = optim(parms0, competitionMin) # give optim initial guess and function
parms = optimOut$par # extract parameters
nSim = as.data.frame(lsoda(n0,times,competition,parms)) # run ode to get simulated populations

# plot results - actual and simulated data
pdf(file="competition.pdf")
matplot(times, nTrue, type="p", xlab="Time", ylab="Abundance", pch="o", col=1:2)
lines(nSim$time, nSim$n1, col=1)
lines(nSim$time, nSim$n2, col=2)
legend("topleft", c("N1", "N2"), lty=1, col=1:2)
dev.off()