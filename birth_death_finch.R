library(ape)
library(BioGeoBEARS)

#R code for a birth-death model
#This model is a model which allows species to speciate with a rate lamba, as well as to die with a rate mu
#In this model, clades are expected growing exponentially at the net diversification rate.


#load tree data into R
tree <- read.nexus("http://www.r-phylo.org/w/images/0/02/Geospiza.nex")
birthdeath(tree)


#Birth_rate = lambda
#Death_rate = mu
#Net Diversification Rate (r) = Birth_rate - Death_rate
#Relative Extinction Rate (a)= Death_rate / Birth_rate
lambda <- 5
mu <- 2
a  = mu / lambda
r = lambda - mu

#log likelihood function

log_likelihood <- function(lambda, mu, phy) {
  
  #creat an empty vector for Lnl
 
  a  = mu / lambda
  r = lambda - mu
  #The number of speices in the tree
  N = length(tree$tip.label)
  #branching time of each species
  x <- c(NA, branching.times(tree))
 
 #Set up the condition for the function 
 #If the speciation rate is smaller than 0 or the extinction rate is larger than 1, then return to a very small number   
    
    if (r  < 0 || a  > 1){
      Lnl = exp(0.001)
    }   else{
 #If the speciation rate is larger than 0 or the extinction rate is smaller than 1, then do the log likelihood calculation
 #
  
      
      Lnl = lfactorial(N - 1) + (N - 2) * log(r) + r * sum(x[3:N]) + N * log(1 - a) - 2 * sum(log(exp(r * x[2:N]) - a))
    }
    
print(paste("lambda=",lambda,"mu=", mu, "Lnl=",Lnl)) 
  
 

    
 
  return(Lnl)
  
}

#optimise the birth and death rate
#Find value of lambda and mu by maximizing the log_likelihood given the tree.
optim(c(lambda,mu), fn = function(p)log_likelihood(lambda = p[1],mu=p[2],phy = tree), method=c("L-BFGS-B"),lower = 0, control = list(fnscale = -1))



