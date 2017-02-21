#Assignment 6
#Carter Johnson
#ECL233 Winter 2017

# Consider two approaches to preventing a disease outbreak
# Partial absolute control: control a fraction q of the
# population so that these individuals can infect no others.
# Homogenous partial control: reduce the probability of
# every transmission by a factor q.

#clear workspace
rm(list=ls())

#read in some SARS data, goal:fit some distributions to it
load("SARS.Rdata") #data file lists number of infections from each person in data set

#function to do all of the simulation and computation
run_models=function(q=0.4, ploton=1){ 
  reps=10000   #number of outbreaks
  Tf=10        #number of generations for simulation
  N_init = 1 #initial infected population size

  N1=matrix(NA,Tf,reps)
  N2=matrix(NA,Tf,reps)
  N1[1,] = N_init
  N2[1,] = N_init

  #run a double loop to go through reps and time
  for (i in 1:reps){
    for (j in 1:(Tf-1)){
      if (N1[j,i]>q){ #Model 1 - Partial Absolute Control
        #get random infection counts from SARS data, 
        #number of "trials" should be number of infected individuals in last generation
        #minus the fraction q which have been quarantined, rounding up (since we start w/ 1 unquarantined individual)
        N1[j+1,i]=sum(sample(SARS,size=(ceiling((1-q)*N1[j,i])),replace=TRUE))
      }else{
        N1[j+1,i]=0
      }
      if (N2[j,i]>0){ #Model 2 - Homogenous Partial Control
        #get random infection counts from SARS data, 
        #number of "trials" should be number of infected individuals in last generation
        #probability of infection has been reduced by q for each infection incident, so can
        #include this effect as a factor of 1-q to the total number of infections
        #rounding up the infectious pool for consistency w/ model 1
        N2[j+1,i]=(1-q)*sum(sample(SARS,size=ceiling(N2[j,i]),replace=TRUE))
      }else{
        N2[j+1,i]=0
      }
    }
  }

  #plot histogram of number of infections at gen 10 across experiments
  tenth_gen_N_model1 = N1[10,]
  tenth_gen_N_model2 = N2[10,]
  if (ploton==1){
    pdf(file="histogram_model1.pdf")
    hist(tenth_gen_N_model1, freq=TRUE, main="Histogram of Gen 10, Model 1 (Partial Abs. Control)", xlab="N1[10]")
    title()
    dev.off()
    pdf(file="histogram_model2.pdf")
    hist(tenth_gen_N_model2, freq=TRUE,main="Histogram of Gen 10, Model 2 (Homog. Partial Control)", xlab="N2[10]")
    title()
    dev.off()
  }

  #calculate probability of extinction by gen 10
  nonextinct1 = 0 #number of non-extinctions in gen 10 for model 1
  nonextinct2 = 0 #number of non-extinctions in gen 10 for model 2
  for (i in 1:reps){
    if (tenth_gen_N_model1[i]>0){
      nonextinct1 = nonextinct1+1
    }
    if (tenth_gen_N_model2[i]>0){
      nonextinct2 = nonextinct2+1
    }
  }
  ext_prob_1 = 1-nonextinct1/reps #probability of extinction by gen 10 for model 1
  ext_prob_2 = 1-nonextinct2/reps #probability of extinction by gen 10 for model 2

  return(c(ext_prob_1, ext_prob_2))
}  

# (A) Simulate the model for 10 generations separately
# for each form of control for 10,000 reps and q=0.4. Plot
# histograms for the total number infected 
# (conditioned on non-extinction) at the 10th generation
run_models()


# (B) Create a plot of the probability of extinction by
# generation 10 as function of q for both control strategies.
#run models for various q
dq=0.01                              #q seq step size
q=seq(from=0, to=1, by=dq)          #sequence of q values
ext_probs = matrix(NA, (1/dq+1),2)  #matrix to store ext probs
for (i in 1:(1/dq+1)){
  #run models to get extinction probabilities at gen 10
  ext_probs[i,] = run_models(q=q[i],ploton=0)
}
#plot probability of extinction as a function of q
pdf(file="ext_probs.pdf")
matplot(q, ext_probs, type="l", xlab=expression(q), ylab="Extinction Probability")
title("Extinction Probability at Gen 10 as a function of q")
legend("bottomright", c("Model 1", "Model 2"), lty=1, col=1:2)
dev.off()