#Assignment 9 - Demographic Stochasticity
#Carter Johnson
#ECL233 Winter 2017 - Week 9

#recreate Mean Extinction Time (MTE) curves for the 
# P and NBd models in Fig 3 of the Melbourne and Hastings (2008) reading.  

rm(list=ls())

P_mean_time = function(R){
  #compute the QSD of the Poisson Ricker model and compute the mean time to extinction for given R
  #Poisson Ricker Model
  #N[t+1] = sum_i^N[t] S_i,t+1
  #where S_i,t+1 ~ Poisson(R exp(-aN[t]))
  #so that N[t+1] = Poisson(N[t]Rexp(-aN[t]))
  
  #adjust alpha so that the equilibrium population level is 30
  alpha = log(1/R)/(-30)
  
  # computing the transition matrix for the transient states of the model
  k=100
  Q=matrix(0,k,k)
  for(i in 1:(k)){
    Q[i,]=dpois(1:k,lambda=R*i*exp(-alpha*i))
  }
  # computing the quasi-stationary distribution
  estuff=eigen((Q))
  # Computing the Mean extinction times for this model from the simulation
  lambda=Re(estuff$values[1]) # the probability of persisting 
  # to the next time step given the population is following 
  # the QSD
  return(1/(1-lambda)) # mean time to extinction for QSD
}

NBd_mean_time = function(R, kD=0.5){
  #simulate the Negative Binomial dist. Ricker model and compute the mean time to extinction for given R
  #Neg. Binom. Dist. Ricker Model
  #N[t+1] = sum_i^N[t] S_i,t+1
  #where S_i,t+1 ~ NegBinom(R*exp(-a*N[t]),kD)
  #so that N[t+1] = Poisson(N[t]*R*exp(-a*N[t]),kD*N[t])

  #adjust alpha so that the equilibrium population level is 30
  alpha = log(1/R)/(-30)
  
  # computing the transition matrix for the transient states of the model
  k=100
  Q=matrix(0,k,k)
  for(i in 1:(k)){
    Q[i,]=dnbinom(x=1:k,mu=R*i*exp(-alpha*i),size=kD)
  }
  # computing the quasi-stationary distribution
  estuff=eigen((Q))
  # Computing the Mean extinction times for this model from the simulation
  lambda=Re(estuff$values[1]) # the probability of persisting 
  # to the next time step given the population is following 
  # the QSD
  return(1/(1-lambda)) # mean time to extinction for QSD
}

R = seq(from=1, to=20, by=0.1)
times = matrix(nrow=length(R),ncol=2)
#NBd_times = numeric(length(R))
for (i in 1:length(R)){
  times[i,1] = P_mean_time(R=R[i])
  times[i,2] = NBd_mean_time(R=R[i])
}
pdf(file="mean_time_extinctions.pdf")
matplot(R,log(times), type='l', ylab="Log(Mean Time to Extinction)")
legend("topright", c("P", "NBd"), lty=1:2, col=1:2)
dev.off()
