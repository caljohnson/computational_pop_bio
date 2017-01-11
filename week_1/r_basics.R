#Carter Johnson - ECL 233 Winter 2017
#Assignment 1: R Basics

# plot N_t+1 vs. N_t for Ricker and Beverton-Holt on the same plot for the same K and r values
# Ricker: N_t+1 = N_t*exp(r*(1-N_t/K))
# Beverton-Holt: N_t+1 = (exp(r)*N_t)/(1+(exp(r)-1)*N_t/K)

#clear workspace
rm(list = ls())

#Ricker function - takes in population size N_t, growth rate r, carrying capacity K
ricker = function(N, r, K){
  N2 = N*exp(r*(1-N/K))
  return(N2)
}

#Beverton-Holt function - takes in population size N_t, growth rate r, carrying capacity K
bevertonHolt = function(N, r, K){
  N2 = (exp(r)*N)/(1+(exp(r)-1)*N/K)
  return(N2)
}

#create vector for N_t
N = seq(0,10,0.1)

#set parameters
K = 6
r = 1.2

#PLOT N_t+1 vs N_t
# create file for saving plot as a pdf
pdf(file="r_basics.pdf")

#compute ricker and beverton-holt curves as a matrix and plot as two vectors against N_t
matplot(N, cbind(ricker(N, r, K), bevertonHolt(N,r,K)), type="l", xlab=expression(N[t]), ylab=expression(N[t+1]))

#add legend
legend("topleft", c("Ricker", "Beverton-Holt"), lty=1:2, col=1:2)

#shut down current plot (paired with "pdf" command above)
dev.off()

