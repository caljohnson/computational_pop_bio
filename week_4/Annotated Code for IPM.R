# Integral Projection Models
# ==========================

# This code creates an Integral Projection Model (IPM) for a data set of the species Dracocephalum austriacum 
# (Nicole et al. 2011). Size x is defined as as the logarithm of the product of plant height and number of stems
# The IPM is of the form N[y,t+1]=integral of (s[x]*G[y,x]+F[y,x])*N[x,t]dx
# where 
# N[x,t] is the density of individuals of size x in year t
# s[x] is the survivorship of individuals of size x
# G[y,x] is the infinitesmal probability of surviving individuals of size x growing to a size y
# F[y,x] is the density of offspring of size y produced by individuals of size x

# Goals of today are:
# 1. Read in and plot the data
# 2. Use the GLM package to estimate s[x], G[y,x], and F[y,x]
# 3. Plot the data against the fitted functions
# 4. Create a discretized version of the model using midpoint rule
# 5. Do an elasticity analysis of the discretized model

# clear variables
rm(list = ls())

#####
# Goal 1: Read and plot some data
####
# read the csv file of data
D=read.csv("IPM-data.csv")
# we can use names(D) to see the names of all of the columns of data
# e.g. we can access the "size" data with the command D$size
# Plotting the survival data 
plot(D$size,D$surv,pch=21,bg=rgb(1,0,0,0.5))
# Plotting the growth of surviving individuals
winners=which(D$surv==1)
plot(D$size[winners],D$sizeNext[winners],pch=21,bg=rgb(1,0,0,0.5))
# Plotting the probability of flowering
plot(D$size,D$fec.flower,pch=21,bg=rgb(1,0,0,0.5))
# Plotting the fecundity data
parents=which(D$fec.flower==1)
plot(D$size[parents],D$fec.seed[parents],pch=21,bg=rgb(1,0,0,0.5))

#####
# Goals 2 and 3: Build the components of the IPM and plot
####
# To build these components, we are going to use the GLM package
# GLM stands for Generalized Linear Models
# The generalization over linear regression corresponds to 
# (i) can have non-normally distributed data e.g. Poisson, Binomial
# (ii) nonlinear transformations of the data i.e. "link functions"

# For example, for the survival data, the data consists of 0s and 1s
# as each individual either survived or didn't. Hence, this data
# is best modeled as a binomial random variable where the probability 
# of success (i.e. surviving) is a function of size
# As we expect this probability to increase and saturate with size,
# the default link function is the logit function i.e. the 
# probability of surviving increases in a logistic manner with size

# The GLM package will find parameters that maximize the likelihood
# of the model producing the observed data. We can do this as follows
# for survivorship

model.survive=glm(surv~size,data=D,family=binomial)

# The first argument is the "formula" (i.e. we want to model survival as a function 
# of size). The second argument is the name of our data frame, the final argument
# specifies that we are modeling survivorship using binomials

# To see how well our fit does, we can use the following command

summary(model.survive)

# Next we might want to plot the model against the data
# We need the range of sizes (i.e. our x)
size.range=range(cbind(D$size,D$sizeNext),na.rm=TRUE)
a=size.range[1] # the minimum size
b=size.range[2] # the maximum size
# We need to create a vector of x values for our function
n=100 # number of x values
xs=seq(a,b,length=n) # x values
# Next, we need to define the survival kernel
# input is size and output is probability of surviving
s=function(x)predict(model.survive,data.frame('size'=x),type='response')
ys=s(xs)
# Plotting yields
plot(D$size,jitter(D$surv,factor=0.1),pch=21,bg=rgb(1,0,0,0.5))
lines(xs,ys,lwd=3)


# what about growth???
# The simplest model is to view the data as 
# normally distributed around a line.
# Hence, we can do a regular linear regression
model.growth=glm(sizeNext~size,data=D)
# Checking out the stats shows an excellent fit
summary(model.growth)
# Lets define the mean trend as follows
G.mean=function(x)predict(model.growth,data.frame('size'=x),type='response')
# Plotting things
plot(D$size,D$sizeNext,pch=21,bg=rgb(1,0,0,0.5))
lines(xs,G.mean(xs),lwd=3)

# plot the residuals to see if things are roughly normal
hist(residuals(model.growth))

# create the growth kernel where x is size and y is sizeNext
G=function(y,x){dnorm(y,mean=G.mean(x),sd=sd(residuals(model.growth)))}

# The fecundity kernel
# first, did you have kids or not
model.flower=glm(fec.flower~size,data=D,family=binomial)
# create the flowering function
F.flower=function(x)predict(model.flower,data.frame('size'=x),type='response')
# second, how many kids given that I flowered?
# create a new version of the data that eliminates the non-flowering
data2=data.frame(size=D$size[parents],fec.seed=D$fec.seed[parents])
model.kids=glm(fec.seed~size,data=data2,family=poisson)
F.kids=function(x)predict(model.kids,data.frame('size'=x),type='response')
# plot these
par(mfrow=c(1,2))
plot(D$size,D$fec.flower,pch=21,bg=rgb(1,0,0,0.5))
lines(xs,F.flower(xs))
plot(D$size[parents],D$fec.seed[parents],pch=21,bg=rgb(1,0,0,0.5))
lines(xs,F.kids(xs))
#estimate the establishment probability
est.prob=sum(is.na(D$size))/sum(D$fec.seed,na.rm=TRUE)
mean.size=mean(D$sizeNext[is.na(D$size)])
# estimate the distribution of sizes in first year
sd.size=sd(D$sizeNext[is.na(D$size)])
# Putting together the fecundity function
F.all=function(y,x)F.flower(x)*F.kids(x)*est.prob*dnorm(y,mean=mean.size,sd=sd.size)

# PUT The full kernel TOGETHER!!!!
K=function(y,x)s(x)*G(y,x)+F.all(y,x)

#####
# Goal 4: Create the discretized model
####
# We can approximate the infinite dimensional model
# by creating a finite number of size classes, say n.sizes
n.sizes=100
# Create size bins using midpoints
dx=(b-a)/n.sizes
sizes=seq(a+dx/2,b-dx/2,length=n.sizes)
# Nex we discretize the IPM using the outer command. 
IPM=outer(sizes,sizes,K)*dx
# lets see what the kernel looks like
filled.contour(sizes,sizes,IPM)



#####
# Goal 5: Apply our matrix model methods
####

# Compute the dominant eigenvalue and eigenvector
estuff=eigen(IPM)
lambda=Re(estuff$values[1])
v=Re(estuff$vectors[,1])
v=v/sum(v)
# Compute the reproductive values
w=Re(eigen(t(IPM))$vectors[,1])
w=w/sum(v*w)
# Plot stable size distribution and reproductive values
par(mfrow=c(1,2)) # create two subplots
plot(sizes,v)
plot(sizes,w)

# compute sensitivity and elasticity matrices
S=w%o%v
E=S*IPM/lambda
# plot elasticity matrix
filled.contour(sizes,sizes,E)


