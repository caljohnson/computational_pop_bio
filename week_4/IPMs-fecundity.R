
#######################
# Estimating the fecundity kernel
#######################

#plot the seed data
plot(data$size,data$fec.seed,pch=21,bg=rgb(1,0,0,0.5))
# use a general linear model assuming poisson distributed
model3=glm(fec.seed~size,data=data,family=poisson())
a3=model3$coefficients
# estimate probability of seed establishment
establishment.prob=sum(is.na(data$size))/sum(data$fec.seed,na.rm=TRUE)
# estimate mean recruit size and variation around that mean
recruit.size.mean=mean(data$sizeNext[is.na(data$size)])
recruit.size.sd=sd(data$sizeNext[is.na(data$size)])
# put together the fecundity function
Fecundity=function(x,y) { 		
	establishment.prob*
	dnorm(y,mean=recruit.size.mean,sd=recruit.size.sd)*
	exp(a3[1]+a3[2]*x)
  }



# # create simpson integrating vector

# simpson=(beta-alpha)*c(1,rep(c(4,2),(n-2)/2),4,1)/n

# simpson=outer(rep(1,n),simpson)

# # create discretized IPM

# # IPM=simpson*K.discrete
