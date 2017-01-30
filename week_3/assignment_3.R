#assignment_3.R
#Carter Johnson

#Problem statement: 
#For this HW, you will add density-dependence to the teasel matrix model 
#and plot how the dominant eigenvalue at low densities and 
#the total abundance at equilibria vary with the maximum fecundity 
#of the reproductive individuals.

#To this end, consider the nonlinear matrix model given by 
#N(t+1)=A(N(t))%*%N(t) where N is the population state vector 
#of length 7, and A(N) is given by the teasel matrix TS with 
#the 1,7 entry of 431 replaced by f/(1+a*sum(N)). For this nonlinear model,
#the linearization at the origin (i.e. what happens at low densities) is given
#by the TS matrix with the 1,7 entry replaced by f.

#For these plots use a=0.01 and let f range from 0 to 75. To get estimates of the
#stable equilibrium, use the initial condition N=(2,0,0,0,0,0,0) and iterate 
#the nonlinear model for 500 time steps. Use the final vector of densities 
#as your estimate for the stable equilibrium (Note: one can prove that this 
#model has a globally stable equilibrium). When running this nonlinear model,
#keep in mind that the matrix A(N) changes with each update as it depends on 
#the densities from the previous time step.

#function to alter TS data by adding density dependence
A=function(N,f,a=0.01){
  A.temp=TS;
  A.temp[1,7]=f/(1+a*sum(N));
  return(A.temp) 
  
}

#linearized density dependent TS matrix
linA = function(f,a=0.01){
  A.temp=TS;
  A.temp[1,7]=f;
  return(A.temp) 
  
}

# upload Teasel data
TS = as.matrix(read.table("teasel_stage.txt"))

#Create matrix of population data
#final time step count
tf = 500
#number of stages in Teasel model
nstg = ncol(TS)
#initial population vector
n0 = c(2, rep(0,nstg-1))
#set up matrix for pop data
Nt = matrix(NA,nrow=nstg,ncol=tf)
#fill in inital conditions
Nt[,1] = n0

#create vector for eigenvalues of linearization
#and for the total stable pop sizes
TS_eigvals = rep(0,75)
stable_pops = rep(0,75)

#Vary f from 0 to 75 to obtain a function for the eigenvalue of the linearization
#and a function for the total population size at the stable equilibrium 
for(f in 0:75){
  
  #loop through time to track population growth
  for(t in 1:(tf-1)){
    #update population with density-dependent Leslie matrix (altered TS)
    Nt[,t+1] = A(Nt[,t],f)%*%Nt[,t]
  }

  #calculate stable population as a sum of last pop vector
  stable_pops[f] = sum(Nt[,tf])
  
  #calculate eigenvalues of linearized density dependent TS matrix
  ev = eigen(linA(f))
  TS_eigvals[f] = Re(ev$values[1])
  
}

#plot everything
pdf(file="lin_eigs.pdf")
plot(0:74,TS_eigvals,type="l",col="black",xlab="f",ylab="Linearized Eigenvalues")
abline(1,0)
dev.off()
pdf(file="stable_pop.pdf")
plot(0:74,stable_pops,type="l",col="black",xlab="f",ylab="Stable population size")
dev.off()
