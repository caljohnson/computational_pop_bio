#Class 2/28 - Week 8
#Modeling environmental fluctuations

#=========================================================================
#Finite State Markov chains
#=========================================================================
#Assume there is a finite # of environmental states
#Specify a transition matrix P[i,j], the prob. of going from state i to j
#Simulate using sample.int command

#Simplest is a two-state model (e.g., wet vs dry)
#Write a simulator for this environmental model

#Input: the transition matrix P
#       the length of the run Tf
#       the initial environmental state init.state
#Output: a single run of the model for length Tf

#1 - wet , 2 - dry
#Markov Chain simulator
MC.sim=function(P=rbind(c(0.9,0.1),c(0.5,0.5)), Tf=10, init.state=1){
  States=numeric(Tf)
  States[1]=init.state
  for (t in 2:Tf){
    #pick new state from set of states given transition probabilities of last state
    States[t]=sample.int(n=2, size=1, prob=P[States[t-1],])
  }
  return(States)
}

#run and plot simulation
out=MC.sim(Tf=100)
plot(out, type='b')

#the fraction of time spent in a dry year
Tf=10000
out=MC.sim(Tf=Tf)
counter=cumsum(out-1) #set dry year value to 1, wet to 0 then add up to get cumulative total # of dry years
fraction.dry=counter/1:Tf #normalize by current years count
plot(fraction.dry,type='l')

#things stabilize since this Markov chain converges to a unique stationary distribution
#given by the normalized left eigenvector of P
P=rbind(c(0.9,0.1),c(0.5,0.5))
v = eigen(t(P))$vector[,1]
v = v/sum(v)
abline(h=v[2], col='red')


#=========================================================================================
#Autoregressive process
#=========================================================================================
# Consider an AR(1)
# X[t]=a+b*X[t-1]+sigma*Z[t]
# where Z[t] are iid standard normals
AR1.sim=function(a=0, b=0.5, sigma=1, Tf=10, init.state=0){
  X=numeric(Tf)
  X[1]=init.state
  for (t in 2:Tf){
    #pick new state from set of states given transition probabilities of last state
    X[t]=a+b*X[t-1]+sigma*rnorm(1)
  }
  return(X)
}

#plot a simulation
out=AR1.sim(Tf=1000, b=-0.99)
plot(out,type='l')

#AR1 have a unique stationary distribution if |b|<1
#mean of the stationary process is a/(1-b)
#variance is sigma^2/(1-b^2)

#===========================================================================================
#Create a function that calculates r (growth rate) for the burning grass model
# N(t+1) = A(t+1)N(t)
# where A(t) are drawn from a finite # of matrices
# corresponding to a finite number of environmental states 1,...,k
# and where the environmental states vary stochastically 
# according to a Markov chain model with transition matrix P

# Input: A array of matrices, P transition matrix, Tf length of run to est. r
# Output: r (estimate)

Lyap = function(A, P, Tf){
  k = dim(P)[1] #number of environmental states
  n = dim(A)[1] #number of population states
  N=rep(1/n,n)
  State=1
  multipliers=numeric(Tf)
  for(t in 1:Tf){
    N = A[,,State]%*%N
    #rescale pop size to be 1, keep track of pop change factors
    multipliers[t]=sum(N)
    N=N/sum(N)
    State=sample.int(n=k,size=1,prob=P[State,])
  }
  r = mean(log(multipliers))
  return(r)
}

# Burnt years
Burnt=matrix(c(0.08,1.63,2.42,4.4,
               0.21,0.64,0.35,0.16,
               0,0.19,0.43,0.24,
               0,0.03,0.23,0.48),4,byrow=TRUE)
#Unburnt
Protected=matrix(c(0,0.706,0.391,3.59,
                   0.018,0.158,0.136,0.093,
                   0,0.08,0.07,0.21,
                   0,0,0.01,0.07),4,byrow=TRUE)

#combine into an array
A = array(0,dim=c(4,4,2))
A[,,1]=Burnt
A[,,2]=Protected

#run and print r
print(matrix(0.5,2,2))
print(Lyap(A=A,P=matrix(0.5,2,2),Tf=100000))

#lets end by plotting how r depends on frequency of burnt years 
#assuming that environmental state is uncorrelated in time

#we want to set P=rbind(c(p,1-p),c(p,1-p))
#where p is the probability of burnt year

the.burning=function(p){
  P = rbind(c(p,1-p),c(p,1-p))
  r = Lyap(A=A,P=P,Tf=10000)
  return(r)
}

please.stop=seq(0,1,length=20)
the.class=sapply(please.stop,the.burning)
plot(please.stop,the.class)
abline(h=0,lty=2)