#Assignment 8 - Temporal Correlations
#Carter Johnson
#ECL233 Winter 2017 - Week 8

# Consider a stage structured population with juveniles and adults
# Each year, adults have F(t) kids which enter the juvenile class
# immediately prior to the next year, 50% of juveniles survive and
# mature to adults, and 50% of adults surive. Each adult
# produces a random number, F(t), of offspring. To model
# F(t), use an autoregressive process X(t) given by 
# X(t)=rho*X(t-1)+sigma*sqrt(1-rho^2)*Z(t)
# where Z(t) are iid std normals.
# This form of the autoregressive process has 
# mean 0, variance sigma^2 and temporal autocorrelation rho.
# To ensure that we can separetely the mean of F(t)
# and its variance (determined by sigma), define
# F(t)=exp(a-sigma^2/2+X(t))
# Then the mean of F(t) is given by exp(a)

rm(list=ls())

# X[t] = rho*X[t-1]+sigma*sqrt(1-rho^2)*Z[t]
# where Z[t] are iid standard normals
autoregressive_process=function(rho, sigma=1, Tf=10, init.state=0){
  X=numeric(Tf)
  X[1]=init.state
  for (t in 2:Tf){
    #pick new state from set of states given transition probabilities of last state
    X[t]=rho*X[t-1]+sigma*sqrt(1-rho^2)*rnorm(1)
  }
  return(X)
}

# F(t)=exp(a-sigma^2/2+X(t))
random_kids=function(a=0.25, rho, sigma, Tf=10){
  F = numeric(Tf)
  X = autoregressive_process(rho, sigma, Tf)
  for (t in 1:Tf){
    F[t] = exp(a-sigma^2/2+X[t])
  }
  return(F)
}

#Find growth rate r
calculate_r=function(rho, sigma, surv_rate=0.5, Tf){
  #initialize population, half in stage 1- juvi, half in stage 2-adult
  N=c(0.5,0.5)
  #set up multipliers vector
  multipliers=numeric(Tf)
  #get time series for random offspring function
  F = random_kids(rho=rho,sigma=sigma, Tf=Tf)

  for(t in 1:Tf){
    #get transition matrix for the model for this time:
    #top row- 0 juvis result from juvis, F[t] juvis result from adults
    #bottom row - 0.5 juvis age up to adult, 0.5 adults survive as adults
    P=rbind(c(0,F[t]), c(surv_rate,surv_rate))
    #grow juvinilles, calculate survivorship
    N=P%*%N
    #compute how much bigger total pop size is from 1
    multipliers[t]=sum(N)
    #normalize total pop size to 1
    N=N/sum(N)
  }
  r=mean(log(multipliers))
  return(r)
}

# Create a contour plot of r with a=0.25
# sigma varying between 0 and 1
# rho varying between -0.95 and 0.95
sigma = seq(0,1, by= 0.1)
rho = seq(-0.95, 0.95, by=0.05)
r = matrix(0,length(sigma),length(rho))
for (i in 1:length(sigma)){
  for (j in 1:length(rho)){
  r[i,j] = calculate_r(rho=rho[j],sigma=sigma[i], Tf=10000)  
  }
}

pdf(file="growth_rate_contour.pdf")
filled.contour(sigma,rho, r, xlab="sigma", ylab="rho", main="r contours")
dev.off()

