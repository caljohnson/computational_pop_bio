# goals for today:
# - learn continuous-time integration
# - fit a model to data

rm(list=ls()) # clear everything

# we'll need the deSolve package: install from install packages, then load the package
library(deSolve)

# function to pass to lsoda:
# has to take in time, state variable(s), and parameter vector
# has to return state varaiable dynamics (function you're integrating)
logInt = function(t, n ,parms) {
  with(as.list(parms),{
    dn = r*n*(1-n/K) # logistic dn/dt
    return(list(dn)) # returns dn/dt as a list
  })
}

# function to minimize to get best fit to data: least squares
logIntMin = function(parms) {
  out = as.data.frame(lsoda(n0,times,logInt,parms)) # run logistic with current guess
  mse = mean((out$n-nTrue)^2)# calculate square differences
  return(mse) # return value we want minimized by optim
}

# parameters
r = 0.5 # intrinsic growth rate
K = 100 # carrying capacity
parms = c(r=r, K=K) # named vector of parameters
tf = 20 # end time
times=1:tf # vector of times to run over
n0 = c(n=2) # initial population size, have n= label means the output of lsoda will label population size as n

# run model
# use lsoda: inputs are:
# 1. initial conditions
# 2. time vector to run over
# 3. function to integrate
# 4. parameter vector
# as.data.frame puts the output in the format of a data frame
out = as.data.frame(lsoda(n0,times,logInt,parms))

# plot output
plot(out$time,out$n, type="l",xlab="time",ylab="population size")

# fit to data
nTrue = as.matrix(read.table("paramecium1.txt"))
colnames(nTrue)=NULL # takes out default column names from read.table
n0 = c(n=nTrue[1]) # pull first entry for initial population size
times = 1:length(nTrue)
# guesses for r and K
rguess = 1 # guess for growth rate
Kguess = 500 # guess for carrying capacity
parms0 = c(r=rguess,K=Kguess) # initial parameters vector

# find fit using optim
optimOut = optim(parms0,logIntMin)
# run function with best guess parameters
parms = optimOut$par
nSim = as.data.frame(lsoda(n0,times,logInt,parms))

# plot results
plot(times,nTrue,type="p",xlab="time",ylab="pop size")
lines(times,nSim$n)
