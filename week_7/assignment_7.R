#Assignment 7 - The Storage Effect
#Carter Johnson
#ECL233 Winter 2017

# To understand when all species coexist, it is useful to look 
# at the stochastic growth rate of each species when rare. 
# This growth rate for, say species 1, is given by running the 
# simulation without species 1 and computing the average of 
# the following quantity 
# log(1-d+d*kids[1]/sum(kids*N[t-1,]))
# during the simulation. For the assignment, compute 
# stochastic growth rates when rare for each species in the 
# community when (i) Sigma=diag(10) and 
# (ii) Sigma=diag(10)*0.5+0.5. 
# Create a pair of barplots in a single figure showing 
# these stochastic growth rates for all the species. 
# HINT: To make the simulations run a bit faster, 
# it isn't necessary to store N at every point in time, 
# just keep track of the current vector N of population frequencies

#clear workspace
rm(list=ls())
#load MASS library for multivariate normals
library("MASS")

#function to compute stochastic growth rates by doing lottery model simulation 
#with leaving out species left_out
lottery_leave_one_out=function(k=10,N0=rep(0.1,10),Sigma=diag(10),mu=seq(1,1.25,length=10),Tf=100,d=0.1,left_out){
  #implements the lottery function leaving out species k
  #to calculate the stochastic growth rate of species k
  N0[left_out]=0
  N_old=N0 # initial values
  quantity=rep(0,Tf-1)
  for(t in 2:Tf){ # starting the loop
    kids=exp(mvrnorm(n=1,mu=mu,Sigma=Sigma)) # random draw creating the lognormal number of kids
    N_new=(1-d)*N_old+d*kids*N_old/sum(kids*N_old) # determining the next community state
    quantity[t-1] = log(1-d+d*kids[left_out]/sum(kids*N_old))
    N_old <- N_new
  }
  gr=mean(quantity)
  return(gr) # returning the stochastic growth rate
}

#(i) Sigma=diag(10)
gr1 = rep(0,10)
for (i in 1:10){
  gr1[i]=lottery_leave_one_out(Sigma=diag(10),Tf=10000,left_out=i)
}

#(ii) Sigma=diag(10)*0.5+0.5.
gr2 = rep(0,10)
for (i in 1:10){
  gr2[i]=lottery_leave_one_out(Sigma=(diag(10)*0.5+0.5),Tf=10000,left_out=i)
}

#plot pair of barplots
pdf(file="stoch_growth_rates.pdf")
par(mfrow=c(2,1))
par(mar = rep(4, 4))
barplot(gr1, xlab="species", ylab="stochastic growth rate", main="Sigma=diag(10)")
barplot(gr2, xlab="species", ylab="stochastic growth rate", main="Sigma=diag(10)*0.5+0.5")
dev.off()