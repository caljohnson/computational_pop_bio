# Goals for the day:
# - learn continuous-time numerical integration
# - fit model to data

rm(list=ls()) # clear everything

# Numerically integrating continuous-time functions always
#  takes an extra "library" of R functions: deSolve
# To get this library, go to the menu for
#  Packages & Data -> Package Installer
#  then search for deSolve, select it,
#  and click "Install Selected"

#----------------------------------
# Once you've downloaded the deSolve libarary
#  via the Packages & Data menu,
#  you also have to load it each time
#  you start up R and plan to use it.
# You load it with this command:
library(deSolve) # loads the ODE numerical integration package

#----------------------------------
# functions

# First define the function you're going to numerically integrate
# In this case, the Lotka-Volterra competition model
#  all functions plugged into the numerical integration function 
#  have to have the same format:
#  fnName = function(t, n, parms) {...}
# See that inputs are first time t,
#  then state variable vector n,
#  then parameter values parms in a vector
#   where they're labeled by their names
# To create a vector of parameter values with corresponding names, 
#  you can use (at the prompt to provide an example)
v=c(x=1, y=4, z=0) 
v
#  this creates a vector v
#   with the first entry (1) labeled as x
#   the second entry (4) labeled as y
#   and the third entry (0) labeled as z
# Then you can use the with(as.list(parms), {..})
#  to run any commands using those parameters, 
# For example, with the v vector above:
#  first type in 
# print(x) #, you'll get an error, 
#  but if you do 
with(as.list(v), {print(x)}) #, it should work
# Note also that if you do 
c(a=1, b=c(2,3))
#  R automatically makes up names b1 and b2

# --------------
# back to the script...
# so let's make our function with our list of parameters
# logInt takes in time t, population size n, and the vector of parameters parms to numerically integrate the logistic function with lsoda
logInt = function(t, n, parms) {
	with(as.list(parms), { # extract parameters from parms vector
		dn = r*n*(K-n)/K # logistic dn/dt
		return(list(dn)) # return dn/dt as a list
		})
	}

# -----------------
# parameters
r = 0.5 # intrinsic growth rate
K = 100 # carrying capacity
parms = c(r=r, K=K) # named vector of parameters
tf = 20 # end time
times = 1:tf # vector of times to run over
n0 = c(n=2) # initial population size - having the n= label means that the output of lsoda will label the population size as n

# ----------------
# run model
# use the lsoda function - inputs are:
#  1. initial conditions
#  2. time vector to run over
#  3. function to integrate
#  4. parameter vector.
# Then convert the output to a "data frame" with as.data.frame
#   This lets you access pieces by labels as well as column/row numbers
out = as.data.frame(lsoda(n0, times, logInt, parms))

# at the prompt
out$time # check out the output
out$n 

# back to the script
plot(out$time, out$n, type="l", xlab="Time", ylab="Population size") # plot the output

# ----------------
# Model comparison to data

# We'll be using the R function optim to fit our model to data: find r and K given a time series for population size
# optim takes in a function and gives you the parameter values that minimize that function
# So what function would we want to minimize in order to get the best fit of our model to the data?
# Square difference between the data and model
# Assumption when using least square difference: errors are observation errors, normally distributed
# The function takes in the parameters, runs the model with those parameters, and finds the square difference from the data
# add this to the script under our functions
logIntMin = function(parms) {
  out = as.data.frame(lsoda(n0, times, logInt, parms)) # run ode
  mse = mean((out$n-nTrue)^2) # calculate mean squared error between simulated and original data
  return(mse) # return mean squared error
}

# parameters and data
nTrue = as.matrix(read.table("paramecium1.txt")) # actual data -- don't forget to make sure you're in the right directory
colnames(nTrue) = NULL # take out default column names -- otherwise this gets in the way in lsoda
n0 = c(n=nTrue[1]) # initial population size - us first data point, make sure to label
tf = length(nTrue) # time
times = 1:tf # vector of times to run over
# need to start with guesses at the parameter values to give something for optim to start with
rguess = 1 # guess for growth rate
Kguess = 500 # guess for carrying capacity
parms0 = c(r=rguess, K=Kguess) # set up parameter vector

# find fit and run simulations
optimOut = optim(parms0, logIntMin) # give optim initial guess and function
parms = optimOut$par # extract parameters
nSim = as.data.frame(lsoda(n0,times,logInt,parms)) # run ode to get simulated population

# plot results - actual and simulated data
plot(times, nTrue, type="p", xlab="Time", ylab="Abundance")
lines(nSim$time, nSim$n)

# task: take this to two interacting species
# For this assignment, you have two choices for simulating two interacting species: (a) competition or (b) predator-prey; pick one

# Regardless of your choice, you'll need set up your ode function for multiple species:
# Initialize your population sizes with labels, e.g., n0=c(x=x0, y=y0, ...), then use these labels in your ode function
# odeFn = function(t, n, parms) {
#	with(as.list(parms, n), { # extract parameters from parms vector
#		dx = ... # dx/dt
#		dy = ... # dy/dt
#		...
#		return(list(c(dx, dy, ...))) # return dn/dt as a list
#		})
#	}
# then running lsoda gets you each species as out$x, out$y, etc.

# Choice (a) competition
# fit paramecium data (from Leslie 1957 reading) to Lotka-Volterra competition model
# dn1/dt = r1*n1*(1-n1/K1-a12*n2/K1)
# dn2/dt = r2*n2*(1-n2/K2-a21*n1/K2)
# plot actual and simulated time series for both populations
# you get r's and K's - this is in the "parameciumrk.txt" file - don't forget to use colnames(your input matrix) = NULL to erase the default V1, V2, etc. column names
# you'll have to find the alpha's to the data in the "paramecium2.txt" file
# note that this means you have to reconstruct the parameters vector within optim - start with parmsrk=r and K, then do parms = c(parmsrk, a=a) to add a's within minimizing function

# Choice (b) predator-prey
# Create phase plane diagrams (predator density vs. prey density, where lines are trajectories over time) for the predator-prey mode:
# dN/dt = r*F(N)-b*G(N)*P # prey 
# dPdt = c*G(N)*P-m*P # predator 
# for four scenarios:
# (i) type I functional response (G(N)=N) and no density dependence in prey (F(N)=N)
# (ii) type I functional response (G(N)=N) and density dependence in prey (F(N)=N*(1-N/K))
# (iii) type II functional response (G(N)=N/(1+d*N)) and density dependence in prey (F(N)=N*(1-N/K))
# (iv) type III functional response (G(N)=N^2/(1+d*N^2)) and density dependence in prey (F(N)=N*(1-N/K))
# in addition to plotting trajectories, mark the starting and ending points.  Bonus if you also put a point at the internal equilibrium (equilibrium with nonzero values for both prey and predator) for each model.
# some parameter recommendations, you can alter these as you see fit:
# b = 0.01 # predator attack rate
# c = 0.1*b # predator conversion of preation into reproduction
# m = 0.2 # predator mortality
# a = 1 # rate of prey capture/unit prey and time
# d = 0.0015 # handling time
# r = 0.5 # prey growth rate
# K = 1000 # prey carrying capacity
# N0=300, P0=50 for initial population sizes, although you might want to try a few different initial population sizes per plot

# try together in class: set up Lotka-Volterra competition simulation (don't worry about data fitting just yet)
# some test parameters:
# r = c(0.8, 0.6) # population growth rates
# K = c(600, 200) # carrying capacities
# a = 0.5*c(K[1]/K[2], K[2]/K[1]) # competition coefficients (this will be in a parameter space with coexistence)
# run for 20 time units, starting with two individuals in each population 

# ----------------
