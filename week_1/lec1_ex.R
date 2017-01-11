#In Class Work 1/10/17

#plotting exponential growth for different values of r and n0

#clear workspace
rm(list=ls())

#define a function for exponential growth
#takes in initial population size n0, population growth rate r, 
#and time t (either a final scalar time or vector of times)
expGrowth = function(n0, r, t){
	#exponential growth function
	n = n0*exp(r*t)
	#function output 
	return(n)
}

#plot exponential growth for different values of r and n0, one t series

#parameters
#time series
t = seq(0,10,0.1)
#baseline growth rate
r = 1.2
#baseline initial population size
n0 = 2

#call expGrowth on each combination on the same t series
#baseline exponential growth
n = expGrowth(n0,r, t)
#second r value exp growth
n2 = expGrowth(n0, r*2, t)
#second n0 value exp growth
n3 = expGrowth(n0*2, r, t)

#plot different ivp solutions
# create file for saving plot as a pdf
pdf(file="lec1_ex.pdf")

#create a matrix with these data
y = cbind(n,n2,n3)
# the matrix plot command that plots each column
# options used are type="l" to connect the points (lines), and log="y" to get 
# the y axis on a log scale
matplot(t,y,log="y",type="l", xlab="Time", ylab="Pop size") 
# #R sets the plot window size based on first plot, so use biggest first
# plot(t, n.r2, type="l", col="blue", xlab="Time", ylab="Pop size", log="y")

# #add the other two solutions
# lines(t, n.base, col="black")
# lines(t, n.n02, col="red") 

#add legend
legend("topleft", c("Base", "r doubled", "n0 doubled"), lty=1:3, col=1:3)

# shut down current plot (paired with "pdf" command above)
dev.off() 