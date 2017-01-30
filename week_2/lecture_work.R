# today we're going to be doing loops lol
# we tested out if loops in the prompt
# now let's learn "while" loops lol

#While - performs an action while (XXX) is true, stop as soon as XXX is false

#clear workspace
rm(list = ls())

#Beverton-Holt function
bh = function(nt, alpha=1.2, beta=0.001){
  nt1 = alpha*nt/(1+beta*nt)
  return(nt1)
}

#Run until we hit equilibrium (or close to it) with a while loop

#minimum difference that signals reach equilibrium (reached margin of error)
minDiff=0.01
#maximum time to run
tmax=1000

#initial population size
n0 = 100

#start with two time points
n=rep(NA,2)
t=2
#set first time point
n[1] = n0
#second time point
n[2] = bh(n[1])

while(abs(n[t]-n[t-1])>minDiff & t<tmax){
  #update population size
  n[t+1] = bh(n[t])
  #update time
  t = t+1
  print(n)
  print(t)
}

#plot to show time series
plot(1:t,n, xlab="time", ylab="population size", type="l")


#for loops: perform an action for each value of x from xmin to xmax
#clear n to start over
rm(n)
#final time
tf= 100
#set aside a vector for the time series
n= rep(NA,tf)
#initial population size
n[1]=n0
#for loop
for(t in 1:(tf-1)){
  #step through bevertonholt
  n[t+1]=bh(n[t])
}
plot(1:tf, n, xlab="time", ylab="population size", type="l")

#run beverton-holt for a range of initial population sizes

#vector of initial population sizes
vn0 = seq(10,300, by=10)
#number of initial conditions
Nn0=length(vn0)
#create an empty matrix for time series for all initial conditions
tf=50
n=matrix(NA,nrow=tf,ncol=Nn0)
#fill in initial population sizes
n[1, ] = vn0
#loop through time
for(t in 1:(tf-1)){
  #step through beverton-holt
  n[t+1,]=bh(n[t,])
}

#plot
matplot(1:tf, n, type="l", xlab="time", ylab="population size")