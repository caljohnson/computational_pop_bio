#Lecture 10 - ECL 233 
#Carter Johnson

# SIR Model and Demographic Stochasticity in Continuous Time Models

# dS/dt = µ(S+I+R) - ßSI/N -µS
# dI/dt = ßSI/N - µI - ∂I
# dR/dt = ∂I - µR

# dN/dt = 0, N=S+I+R
# Nondim'liz †=(µ+∂)t
# π = µ/(µ+∂)
# å = ∂/(µ+∂)= 1-π
# r0 = ß/(µ+∂)

# dS/d† = πN - πS - r0(SI/N)
# dI/d† = r0SI/N - (å+π)S
# dR/d† = åI - πR

# Events     - Scaled Rate - Outcome
# birth of S -  πN         - S->S+1, I->I, R->R
# death of S -  πS         - S->S-1, I->I, R->R
# infection  - r0SI/N      - S->S-1, I->I+1, R->R
# recovery   - åI          - S->S, I->I-1, R->R+1
# death of I - πI          - S->S, I->I-1, R->R
# death of R - πR          - S->S, I->I, R->R-1 
#          sum -> total rate = r

# (1) Calculate wait time to next event 
      #exponentially distributed with mean 1/r
# (2) Determine which event happens
      #random uniform draw weighted by relative rates

#=================================================================================================
#GOAL: Implement Demographic Stochasticity in SIR model with births and deaths, compare to deterministic

rm(list=ls())
library(deSolve)

#deterministic SIR model with births and deaths, time non-dimlzd
SIRdet = function(t,x,parms){
  with(as.list(c(parms,x)),{
    dS = p*(N-S)-R0*S*I/N #dS/dt, susceptibles
    dI = R0*S*I/N - I     #dI/dt, infecteds
    #R = N-S-I, dont need diff eq
    list(c(dS,dI))        #return both together
  })
}

#stochastic SIR model
#inputs: disease reproductive number R0, birth/death rate p, total pop N, 
#        initial infecteds I0, total run time tf
SIRstoch = function(R0,p, N,I0,tf=100){
  #error-checking  
  if(p>1) stop("Cannot have birth/death rate more than 1")
  #could also check for negative parm values and variables
  #force all population sizes to be integers
  N = as.integer(N)
  I0 = as.integer(I0)
  if(I0>N) stop("Cannot start with more infecteds than the total pop size")
  
  #initialize: 
  #population sizes of S/I/R
  S0 = floor(0.5*(N-I0))
  Rec0 = N-S0-I0 #no. recovereds to start
  alpha = 1-p #recovery rate
  SIR = as.matrix(c(S0,I0,Rec0)) #matrix for holding SIR pops through time
  #start time 
  tstep=1 #time index
  vt = 0 #time values (will have whole vector after while loop)
  #matrix for updates - order: birth S, death S, infection, recovery, death I, death R 
  dN = cbind(c(1,0,0),c(-1,0,0), c(-1,1,0), c(0,-1,1), c(0,-1,0), c(0,0,-1))
  nOut = ncol(dN) #number of outcomes
  
  #while loop through time
  while (vt[tstep]<tf){
    #determine scaled rates
    Nt = sum(SIR[,tstep]) #total pop size
    #rates for each event - same order as dN matrix
    wt = c(p*Nt,p*SIR[1,tstep], R0*SIR[1,tstep]*SIR[2,tstep]/Nt, alpha*SIR[2,tstep], p*SIR[2,tstep], p*SIR[3,tstep])
    sum_wt = sum(wt) #total rate
    #calculate wait time to next event: exponential distribution based on total wait time
    dt = rexp(n=1, rate=sum_wt)
    #determine which event happens: use sample.int
    outIdx = sample.int(nOut,size=1, prob=wt/sum_wt) #normalize rates to get probs
    #update
    vt[tstep+1]=vt[tstep]+dt #update time
    SIR = cbind(SIR,SIR[,tstep] + dN[,outIdx]) #update pop sizes
    tstep = tstep+1          #update time step counter
  }  
  return(list(t=vt, N=SIR)) #return both times and pop sizes
}

#==================================================================
# Parameters
R0 = 8   #infection rate
p = 0.04 #scaled birth/death rate
N = 1000 #pop size
I0 = floor(0.01*N) #initial infecteds

#==================================================================
#run stochastic
stochOut = SIRstoch(R0=R0, p=p, N=N, I0=I0)

#run deterministic - lsoda
parms = c(R0=R0,p=p,N=N) #parameter vector
x0 = c(S=N-I0, I=I0)     #initial pop sizes
detOut = as.data.frame(lsoda(x0, stochOut$t,SIRdet, parms))

#plot
par(mar = rep(2, 4))
mf=par(mfrow=c(2,1))
plot(stochOut$t, stochOut$N[1,], type="l",xlab="time (disease generations)", ylab="Susceptibles")
lines(stochOut$t, detOut$S, col="blue")
plot(stochOut$t, stochOut$N[2,], type="l",xlab="time (disease generations)", ylab="Infecteds")
lines(stochOut$t, detOut$I, col="blue")
par(mfrow=mf)
