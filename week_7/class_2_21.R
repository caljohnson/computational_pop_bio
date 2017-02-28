#Class - 2/21/17
#Environmental Stochasticity

#clear stuff
rm(list=ls())

#load checkerspot data
load("checkerspot.Rdata")

#plot data
N = N.data$N
N.yr = N.data$years
plot(N.yr, N, log="y", type="b")

rain=precip.data$precip
rain.yr=precip.data$years
plot(rain.yr, rain, type="b")

#model how population size changes as a function of precipitation
y=log(N)
n=length(N.yr)  #number of population densities
temp=N.yr[-1]-N.yr[-n]  #number of years between counts
index=which(temp==1) #data indices with single years between counts
r=y[index+1]-y[index] #realized per-capita growth rates
yrs=N.yr[index]
#want to get the precipitation data for the relevant years where we have pop data
rain.index = rain.yr %in% yrs
w = rain[rain.index]^(-2)
N2=N[index]

#build our little model
model=lm(r~N2+w)
#N[t+1]=N[t]*exp(a+b*N[t]+c*W[t])=f(N[t],W[t])
#take log
#logN[t+1]-logN[t]=a+b*N[t]+c*W[t]

f=function(N,W)N*exp(model$coefficients[1]+model$coefficients[2]*N+model$coefficients[3]*W)

#based on notes, we're going to calculate E[logR(t)] given we're randomly sampling from rain data
W.data=rain^(-2)
R = function(W)exp(model$coefficients[1]+model$coefficients[3]*W)

#pre 1971 data and compute this
index.pre=which(rain.yr<1971)
ER.pre=mean(log(R(W.data[index.pre])))
#post 1970
index.post=which(rain.yr>1970)
ER.post=mean(log(R(W.data[index.post])))

#simulate full model
simulator=function(N.init=1000, Tf=100, reps=1, W=W.data){
  Nt=matrix(N.init, Tf, reps)
  for(t in 1:(Tf-1)){
    temp = sample(W,size=reps,replace=TRUE)
    Nt[t+1,]=f(Nt[t,],temp)
  }
  return(Nt)
}

#run
Nt=simulator(Tf=200, reps=20, W=W.data[index.pre])
matplot(Nt, type="l", log="y")