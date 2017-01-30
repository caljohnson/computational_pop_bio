#assignment2.R
#Carter Johnson

# two plots using Ricker model (N_t+1=N_t*exp(r*(1-N_t/K))):
# 1. plot N vs. t for four different values of r: 1.5, 2.3, 2.6, 3.0 (show different behaviors - one-cycle, two-cycle, four-cycle, chaos)
# 2. bifurcation plot for r from 1.5 to 3.6

# other parameter suggestions (can use different values, your choice):
# K=100
# N0=50
# tf=100 in 1st plot, 500 in 2nd

#clear workspace
rm(list = ls())

#Ricker function - takes in population size N_t, growth rate r, carrying capacity K
ricker = function(N, r, K=100){
  N2 = N*exp(r*(1-N/K))
  return(N2)
}

#Ricker loop
ricker_loop=function(n,r,tf){
  for(t in 1:(tf-1)){
    #step through ricker model
    n[t+1,]=ricker(n[t,],r)
  }
  return(n)
}

#plot N vs. t for four different values of r: 1.5, 2.3, 2.6, 3.0 (show different behaviors - one-cycle, two-cycle, four-cycle, chaos)

#set vector of r values
r = c(1.5,2.3,2.6,3.0)
Nr = length(r)

#create an empty matrix for time series for all initial conditions
tf=100
n=matrix(NA,nrow=tf,ncol=Nr)

#set initial pop sizes
N0=50
n[1,] = N0

#step through ricker model
n=ricker_loop(n,r,tf)

# create file for saving plot as a pdf
pdf(file="vary_r_plot_ricker.pdf")
#plot
matplot(1:tf, n, type="l", xlab="time", ylab="population size")
#add legend
legend("topleft", c("r=1.5", "r=2.3", "r=2.6", "r=3.0"), lty=1:4, col=1:4)
#shut down current plot (paired with "pdf" command above)
dev.off()


# 2. bifurcation plot for r from 1.5 to 3.6

#set vector of r values
r = seq(1.5,3.6, by=.0001)
Nr = length(r)

#create an empty matrix for time series for all initial conditions
tf=500
n=matrix(NA,nrow=tf,ncol=Nr)

#set initial pop sizes
N0=50
n[1,] = N0

#step through ricker model
n=ricker_loop(n,r,tf)

#capture steady states/cycles of system
#set empty matrix to capture last 20 states for each simulation with value r
ss = matrix(NA,nrow=Nr,ncol=20)
for(i in 0:20){
  ss[,i]=n[tf-i,]
}

# create file for saving plot as a pdf
pdf(file="ricker_bif_diag_r.pdf")
#plot
matplot(r, ss, type="p", xlab="r", ylab="N*", pch=".", col=1:1)
#shut down current plot (paired with "pdf" command above)
dev.off()

