#Class - 2/14/17
#Beginning of Probability Half of Course

#clear
rm(list=ls())

#binomial distribution
rbinom(n=1, size=240000,prob=0.5) #random number of heads in 240000 coin flips, experiment outcome
dbinom(12000,size=24000,prob=0.5) #likelihood of obtaining value 12000 in above experiment
sum(dbinom(12102:24000, size=24000, prob=0.5)) #prob of getting at least 12101 heads
#OR use 1-cumulative prob.
1-pbinom(12101, size=24000, prob=0.5) #pbinom gives cumulative distribution, prob getting 12101 or less

#what does the distribution look like?
no.trials=24000
barplot(dbinom(0:no.trials, size=no.trials, prob=0.5))
#it looks normal because of CENTRAL LIMIT THEOREM

#CLT exploration
#choose 5 numbers
my.data=c(9,11,42,911,230)
#10,000 sampling experiments with 100 trials
no.trials=100
k=10000
#sample no.trials*k w/ replacement from my.data
X = sample(my.data, size=no.trials*k,replace=TRUE)
#rewrite X as matrix with 100 rows and 10000 columns - each column is one experiment
XX = matrix(X,nrow=100, ncol=10000)
#sums of columns
sum.XX = colSums(XX)
#histogram of sums - show density (amount of area in each box)
hist(sum.XX,freq=FALSE)

#get the right Gaussian to plot over this
meanX = mean(my.data)*no.trials
varX = mean((my.data-mean(my.data))^2)*no.trials
sdX = sqrt(varX)
#plot on top of histogram
xs=seq(min(sum.XX),max(sum.XX),length=100)
lines(xs, dnorm(xs,mean=meanX,sd=sdX))

#Poisson distribution deals in discrete integers - can be approx by a binom distribution
lambda=5
n=120
plot(dbinom(0:(3*lambda),size=n,prob=lambda/n),type="b",col="red")
lines(dpois(0:(3*lambda),lambda=lambda),type="b")

#read in some SARS data, goal:fit some distributions to it
load("SARS.Rdata") #data file lists number of infections from each person in data set
require(fitdistrplus)

model1=fitdist(SARS,distr="pois",method="mle")
model2=fitdist(SARS,distr="nbinom",method="mle")

cdfcomp(list(model1,model2),legendtext = c("Poisson", "Negative Binomial"))

#simulate a branching process (disease outbreak)

#sample from the SARS data to determine the number of folks infected by each inidividual in
#a given generation of the outbreak

reps=500 #number of outbreaks
Tf=15 #number of generations for simulation

N=matrix(NA,Tf,reps)
N[1,] = 1

#run a double loop to go through reps and time
for (i in 1:reps){
  for (j in 1:(Tf-1)){
    if (N[j,i]>0){
    N[j+1,i]=sum(sample(SARS,size=N[j,i],replace=TRUE))
    }else{
      N[j+1,i]=0
    }
  }
}

#plot 
matplot(N,type="b", pch=".", log="y")