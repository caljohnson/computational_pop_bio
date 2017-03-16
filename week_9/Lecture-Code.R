# This week, we examine discrete-time models accounting for
# demographic stochasticity. These models are all Markov chain (MC)
# models with a finite or countable number of states typically
# corresponds to the number of individuals in a given state or species.
# We will consider two models: a stochastic Ricker model and 
# Hubbell's neutral model of competing species. 

rm(list = ls())


#===========================================
# Stochastic ricker model
#------------------------------------------- 

# We will use the following data on Great tits (Parus major) from the GPDD to 
# parmaeterize a stochastic ricker model of the form 
# N[t+1] is poisson distributed with mean N[t]*exp(a+b*N[t])
# As the sum of independent Poissons are Poisson, this model
# has a simple interpretation. Namely, each individual
# independently replaces themselves with a Poisson number of 
# individuals with mean exp(a+b*N[t])

N.data=c(21,30,31,32,20,20,21,31,27,24,49,27,41,51,86,43,39,54,46,45,32,33,38,30,46,36,39,26,33,14,30,26,19,27,39,28,46,39,30,33,48,26,35,18)

# plot the data to see what we got 

plot(N.data,type="b")

# To estimate the parameters a and b, define
n=length(N.data)
x=N.data[-n]
y=N.data[-1]

# We can ask what values of a and b maximize the likelihood
# of observing y[i] if y[i] are given by a poisson random
# variable with mean exp(a+b*x[i])*x[i]? To this end, 
# we define the negative of the log likelihood function 
nLL=function(a=0.1,b=-0.1){
	return(-sum(dpois(y,lambda=x*exp(a+b*x),log=TRUE)))
}

# We use Ben Bolker's MLE package to 
# do model selection

require(bbmle)

model1=mle2(nLL,start=list(a=0.1,b=-0.1))
summary(model1)

# from this we can define the best fitting function
a=coef(model1)[1]
b=coef(model1)[2]
fit=function(x)x*exp(a+b*x)

# plotting the fit against the data in two ways; first as one step ahead predictions and then plotting n[t] vs n[t+1]

par(mfrow=c(1,2))
plot(N.data,type="b")
x=N.data[-length(N.data)]
for(i in 1:length(N.data))lines(x=c(i,i+1),y=c(x[i],fit(x[i])),col="red")

y=N.data[-1]
plot(x,y,type="p")
xs=seq(min(x),max(x),length=100)
lines(xs,fit(xs),col="red")

# seeing whether a model of the form N[t+1]=rpois(N[t]*exp(a)) does as well

nLL2=function(a=0.1){
	x=N.data[-length(N.data)]
	y=N.data[-1]
	return(-sum(dpois(y,lambda=x*exp(a),log=TRUE)))
}

model2=mle2(nLL2,start=list(a=0.1))

# computing AICc values for both models, the ricker model has a much lower AICc score, so it "wins"

AIC(model1)
AIC(model2)

# simulating the better model fit

Tf=5000
reps=1
N=numeric(Tf)
N[1]=N.data[1]

for(t in 2:Tf){
	N[t]=rpois(n=1,lambda=fit(N[t-1]))
	
}

plot(N)


# computing the transition matrix for the transient states of the better fitting model

k=100
Q=matrix(0,k,k)
for(i in 1:(k)){
	Q[i,]=dpois(1:k,lambda=fit(i))
}


# computing the quasi-stationary distribution

estuff=eigen(t(Q))

qsd=Re(estuff$vectors[,1])
qsd=qsd/sum(qsd)


# plotting the histogram of the simulation against the QSD
hist(N,freq=FALSE,40,col="red")
lines(1:k,qsd,lwd=2)


# Computing the Mean extinction times for this model

Means=solve(diag(k)-Q)%*%rep(1,k)
plot(Means[1:20])

# The intrinsic mean time to extinction
# corresponds to the mean time to extinction
# given that the population is following the
# quasi-stationary distribution. 

lambda=Re(estuff$values[1]) # the probability of persisting 
# to the next time step given the population is following 
# the QSD
1/(1-lambda) # mean time to extinction for QSD
# versus
Means[20] # from the simulation


# Next we consider Hubble's neutral model
# of competing species. In the "local" version
# of this model there S species and a total of N
# microsites which are all occupied by at least one 
# of the species. Species are demographically 
# equivalent. Hence, all individuals have an equal likelihood
# of dying (d) per year. These empty microsites 
# are colonized by immigrating individuals with
# probability I and by offspring of local individuals
# with the complementary probability 1-I. Consistent
# with neutrality, we assume that the colonizing immigrants
# are equal likely to come from any species. 
# When a non-immigrating individual colonizing a site, 
# the likelihood it is an offspring of an individual of species
# i is propotional to the local abundance X[t,i] of species i. 


d=0.1 # probability die
N=100 # number of microsites
S=10 # number of species
I=0.001 # immigration probability
Tf=100
hub=function(d=0.1,N=100,S=10,I=0.01,Tf=10){
	X=matrix(0,nrow=Tf,ncol=S)
	X[1,]=rep(N/S,S)
	
	for(t in 2:Tf){
		there=which(X[t-1,]>0)
		deaths=rbinom(n=length(there),size=X[t-1,there],prob=d)
		p=I*rep(1/S,S)+(1-I)*X[t-1,]/N
		colonizations=rmultinom(n=1,size=sum(deaths),prob=p)
		X.new=X[t-1,]
		X.new[there]=X.new[there]-deaths
		X.new=X.new+colonizations
		X[t,]=X.new
	}
	return(X)
}


# plotting time series
X=hub(Tf=10000,N=1000,S=100)
par(mfrow=c(1,2))
matplot(X,type="l")

# plotting species rank abundance curve at end of sim
species=which(X[10000,]>0)
rank=sort(X[10000,species],decreasing=TRUE)
plot(rank,log="y",pch=21,bg=rgb(0,0,1,0.8))

# estimating true diversity of community

entropy=function(p){k=length(p);
  temp=which(p>0)
  p.temp=p[temp]
  return(-sum(p.temp*log(p.temp)))}

true.diversity=function(X){
  Tf=dim(X)[1]
  k=dim(X)[2]
  entro.data=numeric(Tf)
  for(t in 1:Tf)entro.data[t]=entropy(X[t,]/sum(X[t,]))
  return(exp(mean(entro.data)))
}

true.diversity(X)

# what effect do you think I has on the true diversity? 

Is=seq(0.001,1,length=20)
temp=function(I)true.diversity(hub(Tf=10000,N=1000,S=100,I=I))
out=sapply(Is,temp)
plot(Is,out,type="b")



# HW recreate Mean Extinction Time (MTE) curves for the 
# P and NBd models in Fig 3 of the Melbourne and Hastings (2008) reading.  




