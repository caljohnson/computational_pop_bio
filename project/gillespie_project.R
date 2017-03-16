#Carter Johnson
#ECL 233 Term Project
#Gillespie Algorithm to introduce Demographic Stochasticity into
#Two-Population Variable-Allee-effect Model

#======================================================================================
#         Intro
#======================================================================================

#dx1 = x1(k1-(x1-1)^2) + D(x2-x1)
#dx2 = x2(k2-(x2-1)^2) + D(x1-x2)

# Events                   -    Scaled Rate     -          Outcome
# Allee-effected birth x1  - (k1+2x1-x1^2)x1    -      x1->x1+1, x2->x2
# Natural death of x1      - x1                 -      x1->x1-1, x2->x2
# Allee-effected birth x2  - (k2+2x2-x2^2)x2    -      x1->x1, x2->x2+1
# Natural death of x2      - x2                 -      x1->x1, x2->x2-1
# Migration of x1          - Dx1                -      x1->x1-1, x2->x2+1
# Migration of x2          - Dx2                -      x1->x1+1, x2->x2-1
#          sum -> total rate = R

rm(list=ls())
library(deSolve)

#======================================================================================
#             Functions
#======================================================================================

#==================================================================
#Deterministic Model (rs=1, gammas=1)
det_model = function(t,x,parms){
  with(as.list(c(parms,x)),{
    dx1 = x1*(k1-(x1-1)^2) + D*(x2-x1) #dx1/dt, pop 1
    dx2 = x2*(k2-(x2-1)^2) + D*(x1-x2) #dx2/dt, pop 2
    list(c(dx1,dx2))        #return both together
  })
}

#==================================================================
#Stochastic Model
#inputs: environmental parameters k1, k2, initial pop sizes of X1 and X2 
#        migration level D, max run time tf, extinction tolerance eps
stoch_model = function(k1,k2,D,X1_0, X2_0, tf=1000, eps=1.0e-1){
  #using r=1, update matrix scaled from 1 to 0.1
  
  pops = as.matrix(c(X1_0,X2_0)) #matrix for holding pops through time
  #start time 
  tstep=1 #time index
  vt = 0 #time values (will have whole vector after while loop)
  #matrix for updates - order: birth x1, death x1, birth x2, death x2, migration x1, migration x2
  dN = 0.1*cbind(c(1,0),c(-1,0),c(0,1),c(0,-1), c(-1,1), c(1,-1))
  nOut = ncol(dN) #number of outcomes
  
  
  #while loop through time
  while (vt[tstep]<tf){
    #determine scaled rates
    x1 = pops[1,tstep]  #current x1 pop size
    x2 = pops[2,tstep]  #current x2 pop size
    births_x1 = max(0, (k1+2*x1-x1^2)*x1) #use max to avoid negative birth rate
    births_x2 = max(0,(k2+2*x2-x2^2)*x2)
    #rates for each event - same order as dN matrix
    wt = c(births_x1, x1,births_x2, x2, D*x1, D*x2) 
    sum_wt = sum(wt) #total rate
    
    #calculate wait time to next event: exponential distribution based on total wait time
    dt = rexp(n=1, rate=sum_wt)
    
    #determine which event happens: use sample.int
    outIdx = sample.int(nOut,size=1, prob=wt/sum_wt) #normalize rates to get probs
    
    #update
    pops = cbind(pops,pops[,tstep] + dN[,outIdx]) #update pop sizes
    #check for extinction - break condition
    #print(pops[,tstep])
    if(pops[1,tstep+1]<=eps){
      #print("x1 went negative")
      #print(births_x1)
      pops[1,tstep+1]=0
    }
    if(pops[2,tstep+1]<=eps){
      #print("x2 went negative")
      #print(births_x2)
      pops[2,tstep+1]=0
      if (pops[1,tstep+1]==0){
        #print("extinction")
        vt[tstep+1]=vt[tstep]
        break}
    }
    vt[tstep+1]=vt[tstep]+dt #update time
    tstep = tstep+1          #update time step counter
  }  
  return(list(t=vt, N=pops)) #return both times and pop sizes
}

#=========================================================================================
# Comparing Stochastic and Deterministic Models
compare_models = function(k1, k2, D, X1_0, X2_0, plot_title="Title"){
  #compares stochastic and deterministic models
  #via plotting the two side by side for both X1 and X2
  
  #run stochastic
  stochOut = stoch_model(k1=k1,k2=k2,D=D, X1_0=X1_0, X2_0=X2_0)
  ext_time = tail(stochOut$t, n=1)
 
  #run deterministic - lsoda
  parms = c(k1=k1, k2=k2,D=D) #parameter vector
  x0 = c(x1=X1_0, x2=X2_0)     #initial pop sizes
  # times = seq(0, ext_time, by=0.01)
  detOut = as.data.frame(lsoda(x0, stochOut$t,det_model, parms))

  #print extinction time
  print("extinction time")
  print(ext_time)

  #plot
  pdf(file=paste(plot_title,".pdf"))
  #par(mar = rep(2, 4))
  mf=par(mfrow=c(2,1))
  plot(stochOut$t, stochOut$N[1,], type="l",xlab="time (generations)", ylab="X1")
  lines(stochOut$t, detOut$x1, col="blue")
  title(main=plot_title)
  plot(stochOut$t, stochOut$N[2,], type="l",xlab="time (generations)", ylab="X2")
  lines(stochOut$t, detOut$x2, col="blue")
  par(mfrow=mf)
  dev.off()
}
#=========================================================================================
# Mean Time to Extinction Computations
mean_extinction_times = function(k1, k2s, D=0.1, X1_0=1.5, X2_0=2.5, tf=1000, no.sims=500){
  #Computes the mean times to extinction for the 2-pop model
  #with a fixed K1 and variable K2
  
  #get number of experiments
  no.experiments=length(k2s)
  mean_times = numeric(no.experiments)
  
  #simulate each expiriment no.sims times
  ext_time = numeric(no.sims) #vector to hold extinction times
  
  #loop over experiments (finding mean extinction time given a K2 value)
  for (j in 1:no.experiments){
    #loop over simulations to find mean extinction time
    for (i in 1:no.sims){
      #run stochastic model until extinction with K2=k2s[j]
      stochOut = stoch_model(k1=k1,k2=k2s[j],D=D, X1_0=X1_0, X2_0=X2_0, tf)
      #get extinction time
      ext_time[i] = tail(stochOut$t,n=1)
    }
    #get mean extinction time for given K2
    mean_times[j]=mean(ext_time)
  }
  #plot how mean times to extinction depends on K2
  pdf(file="mean_time_extinctions.pdf")
  plot(k2s,mean_times, xlab="K2", ylab="Mean Time to Extinction", type="l")
  dev.off()
}

#======================================================================================
#         Running Simulations
#======================================================================================

# Run Stochastic vs Deterministic Model Comparisons
# for various Allee-effect pairings

# Weak-Fatal Pairing Parameters
k1 = 1.5 #environmental parameters of site 1
k2 = -0.25 #environmental parameters of site 2
X1_0 = 1.5 #initial pop 1 size
X2_0 = 2.5  #initial pop 2 size
D = 0.1 #migration level
#compare_models(k1=k1, k2=k2, X1_0 = X1_0, X2_0=X2_0, D=D, plot_title="Weak-Fatal Pairing")

# Strong-Strong Pairing Parameters
k1 = 0.5 #environmental parameters of site 1
k2 = 0.5 #environmental parameters of site 2
X1_0 = 1.5 #initial pop 1 size
X2_0 = 1.5  #initial pop 2 size
D = 0.1 #migration level
#compare_models(k1=k1, k2=k2, X1_0 = X1_0, X2_0=X2_0, D=D, plot_title="Strong-Strong Pairing")

# Weak-Strong Pairing Parameters
k1 = 1.5 #environmental parameters of site 1
k2 = 0.25 #environmental parameters of site 2
X1_0 = 1.5 #initial pop 1 size
X2_0 = 2.5  #initial pop 2 size
D = 0.1 #migration level
#compare_models(k1=k1, k2=k2, X1_0 = X1_0, X2_0=X2_0, D=D, plot_title="Weak-Strong Pairing")

# Strong-Fatal Pairing Parameters
k1 = 0.5 #environmental parameters of site 1
k2 = -0.25 #environmental parameters of site 2
X1_0 = 1.5 #initial pop 1 size
X2_0 = 2.5  #initial pop 2 size
D = 0.1 #migration level
#compare_models(k1=k1, k2=k2, X1_0 = X1_0, X2_0=X2_0, D=D, plot_title="Strong-Fatal Pairing")

# Loss of Resilience in Strong-(Fatal/Strong) System
#    K1 Strong, vary K2 from Fatal to Strong
#    use mean time to extinction to examine effect on resilience
k1 = 0.25 #strong allee effect
k2s = seq(-0.5,0.25, by=0.025) #vary K2 from Fatal to Strong
#run mean extinction time computation across these K2 values and plot
mean_extinction_times(k1=0.25, k2s=k2s, no.sims=1500, tf=2000)



