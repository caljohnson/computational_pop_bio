
###################
# Estimating the survival kernel
###################

#plot the survivorship data using the jitter option to make the values clearer
plot(data$size,jitter(data$surv,factor=0.1),pch=21,bg=rgb(1,0,0,0.5))

# model survivorship with a logistic function 

model2=glm(surv~size,data=data,family=binomial)
a2=model2$coefficients
Survival=function(x){1/(exp(-a2[1]-a2[2]*x)+1)}

# add the plot of the function to the data

lines(xs,Survival(xs),lwd=4,col="blue")