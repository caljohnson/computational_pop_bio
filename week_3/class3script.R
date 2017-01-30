# 1. upload Teasel data
# 2. run for 50 time steps staring with 2 individuals in stage 1
# 3. plot actual growth factors over time: total N[t+1]/N[t]
# 4. plot expected growth factor over time: leading eigenvalue
# how much deviation do you get from not being in stable age distribution

# upload Teasel data
TS = as.matrix(read.table("teasel_stage.txt"))

# run through time
nstg = ncol(TS) # number of stages in the Teasel model
n0 = c(2, rep(0,nstg-1)) # initial population vector
tf = 50 # final run time
Nt = matrix(NA,nrow=nstg,ncol=tf)
Nt[,1] = n0 # fill in initial conditions# for loop through time
for(t in 1:(tf-1)) {
  Nt[,t+1] = TS%*%Nt[,t]
}

# calculate annual population growth factors
Ntot_t = colSums(Nt) # total population sizes through time
gt = Ntot_t[2:tf]/Ntot_t[1:(tf-1)] # growth factor
# alternate: Ntot_t[-1]/Ntot_t[-tf] # take out 1st entry in numerator, take out last entry in demoninator
lambda = Re(eigen(TS)$values[1]) # leading eigenvalue
# plot
plot(1:(tf-1),gt,type="l",col="black",xlab="Time",ylab="Growth factor")
lines(1:(tf-1), rep(lambda,tf-1),col="red")

#------------ now Sebastian: elasticities
# Goal #5. test approximating at time tf-1
# compute approximation
v = Re(eigen(TS)$vectors[,1]) # right eigenvector: stable age distribution
v = v/sum(v)# normalize
w = Re(eigen(t(TS))$vectors[,1]) # left eigenvector: reproductive values
w = w/sum(w*v) # normalize w*v
Napprox = sum(w*n0)*v*lambda^(tf-1) # approximate population size
# plot approximation against actual simulation
barplot(cbind(Napprox,Nt[,tf]),beside=TRUE)

# goal 6: sensitivity analysis
# create matrix of sensitivities i.e. ((d[lambda]/d[L_ij])) = S
S = w%o%v # outer productive of w and v: matrix of sensitivities
barplot(S,beside=TRUE)  # visualize with a bar plot
# create matrix of elasticities
E = S*TS/lambda
# can check: sum(E) should =1
barplot(E,beside=TRUE) # visualize elasticities
