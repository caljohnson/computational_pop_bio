#Assignment 4 - Critical Domain Size
#Carter Johnson

# Consider a population living along a one-dimensional 
#transect for which only the interval [-L,L] is habitable. 
#At the beginning of each year, individuals living in this interval, 
#reproduce and die with each individual producing R offspring. 
#Each of the offspring disperse following a Laplacian dispersal 
#kernel i.e. the probability of dispersing from x to y is
# k(y,x)=exp(-abs(y-m-x)/a)/(2*a)
# 
# where m is the mean displacement and 2a^2 is the variance in the displacement.
#Offspring that land in the interval [-L,L] live to reproduce in the next year.
#Offspring landing outside of this interval die.
# 
# The integral difference equation for this model is given by
# 
# N[t+1,y]=integral from -L to L of R*k(y,x)N[t,x]dx
# 
# For the assignment, assume that R=4, a=2, and m=1.

#clear workspace
rm(list = ls())

#Create a Laplacian dispersal kernel
k=function(y,x, a=2, m=1){ exp(-abs(y-m-x)/a)/(2*a)}

#Create full integral kernel
K=function(y,x, R=4){ R*k(y,x)}

# Create a finite number of locations, store as a vector
locations=function(L, n.locs){
  # Create size bins using midpoints
  dx=2*L/n.locs
  locs_vec=seq(-L+dx/2,L-dx/2,length=n.locs)
  return(locs_vec)
}

# Approximate the infinite dimensional IPM with a finite matrix
IPM=function(L, n.locs){
  # Create size bins using midpoints
  dx=2*L/n.locs
  locs_vec=seq(-L+dx/2,L-dx/2,length=n.locs)
  # Discretize the IPM using the outer command. 
  IPM=outer(locs_vec,locs_vec,K)*dx
  # Visualize kernel
  #filled.contour(locs_vec,locs_vec,IPM)
  return(IPM)
}

# (1) For L=10, create a plot of the population densities at times 
# T=1,2,..,5 and N[1,x]=1 for -2<x<2
tf = 5
L=10
n.locs=100
#create locations vector
locs_vec = locations(L, n.locs)

#create IPM matrix
IPM.matrix = IPM(L, n.locs)

#create population matrix, rows are times, cols are locations
N = matrix(0,nrow=tf,ncol=n.locs)

#set initial population density at time t=1 to N[1,x]=1 for -2<x<2, 0 outside
for(x in (n.locs/2-n.locs/10):(n.locs/2+n.locs/10)){N[1,x]=1}

#loop over time
for(t in 1:(tf-1)){
  N[t+1,] = IPM.matrix%*%N[t,]
}

# create file for saving plot as a pdf
pdf(file="pop_densities_vs_time.pdf")
#plot population densities, turn locations vector into a column vector, transpose N
matplot(cbind(locs_vec), t(N), type="l", xlab=expression(x), ylab="N[t,x]")
#add legend
legend("topleft", c("t=1", "t=2", "t=3", "t=4","t=5"), lty=1:5, col=1:5)
#shut down current plot (paired with "pdf" command above)
dev.off()



# (2) Create a plot of the dominant eigenvalue as a function of L 
#for L between 1 and 5. Add a dashed line at height 1.
dL = 0.01
L = seq(1-dL, 5+dL, by=dL)

#number of dominant eigenvalues to calculate
n.eigens = length(L)
eigvals = rep(0,n.eigens)

#loop over L values ells
for(i in 1:n.eigens){
  n.locs=100
  #create IPM matrix
  IPM.matrix = IPM(L[i], n.locs)
  #calculate dominant eigenvalues of the IPM matrix
  ev = eigen(IPM.matrix)
  eigvals[i] = Re(ev$values[1])
}

# create file for saving plot as a pdf
pdf(file="IPM_eigvals.pdf")
#plot dominant eigenvalues of IPM matrix as function of L
plot.window(c(1,5), c(0,3))
plot(L, eigvals, type="l", xlab=expression(L), ylab="Dominant Eigenvalue of IPM matrix")
abline(1,0, col=1, lty=2)
#shut down current plot (paired with "pdf" command above)
dev.off()