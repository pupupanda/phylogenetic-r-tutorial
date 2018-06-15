library(ape)
library(BioGeoBEARS)


#R code for the ML for a birth-death model
#This model is a model which allows species to speciate with a rate lamba, as well as to die with a rate mu
#In this model, clades are expected growing exponentially at the net diversification rate.


#load tree data into R
tree <- read.tree("elopomorph.tre")



#Birth_rate = lambda
#Death_rate = mu
#Net Diversification Rate (r) = Birth_rate - Death_rate
#Relative Extinction Rate (a)= Death_rate / Birth_rate
#lambda <- 5
#mu <- 2
#a  = mu / lambda
#r = lambda - mu

#log likelihood function

log_likelihood <- function(lambda, mu, tree) {
  
  #creat an empty vector for Lnl
 
  a  = mu / lambda
  r = lambda - mu
  #The number of speices in the tree
  N = length(tree$tip.label)
  #branching time of each species
  x <- c(NA, branching.times(tree))
 plot(tree)
 #Set up the condition for the function 
  #If the speciation rate is smaller than 0 or the extinction rate is larger than 1
  #then there won't be a tip to be observed. 
  #return the log likelehood to a very small number   
    Lnl = 0
    if (r  < 0 || a  > 1){
      Lnl = exp(0.001)
    }   else{
  #If the speciation rate is larger than 0 or the extinction rate is smaller than 1, then do the log likelihood calculation
  
      ##The log_likelihood calculation break down
      
      # There are (n-1)! possible topologies for any set of n-1 waiting times
      # Log (factorial(numtips-1))
      
      # the number of internal nodes is the number of tips -2 (n-2)
      # the expected waiting time is net diversification rate 
      # log expected waiting time  for internal branches
      # (N - 2) * log(r) 
       
      
      
      # total length  of internal nodes not including root: sum(x[3:N])
      # total amount of diversification on these 2 branches
      # r * sum(x[3:N]) 
      
      # N is the number of tip which means there is no speciation nor extinction.
      # log(Waiting time of nothing happen)
      # N * log(1 - a) 
      
      # there are two branches associatied with each node
      # - 2 * sum(log(exp(r * x[2:N]) - a))
  
      # ages of each internal branches
      # x[2:N]
      
    
      # net diversification from each node
      # r * x[2:N]
      
      # expected amount of diversification from each node
      # exp(r * x[2:N])
      
      # expected amount of diversification from each node, 
      # minus number of lineages expected to go extinct from each node
      # exp(r * x[2:N]) - a
      
      
      # amount of waiting time from each internal node
      # sum(log(exp(r * x[2:N]) - a))
      
      
  
      
      Lnl = lfactorial(N - 1) + (N - 2) * log(r) + r * sum(x[3:N]) + N * log(1 - a) - 2 * sum(log(exp(r * x[2:N]) - a))
    }
    
print(paste("lambda=",lambda,"mu=", mu, "Lnl=",Lnl)) 
  
 

    
 
  return(Lnl)
  
}

#optimise the birth and death rate
#Find the maximum value of lambda and mu by using optim on  the log_likelihood function
optim(c(0.01,0.001), fn = function(p)log_likelihood(lambda = p[1],mu=p[2],tree = tree), method="L-BFGS-B",lower = 0, control = list(fnscale = -1))

#output
#$par
#[1] 0.02973697 0.00000000

#$value
#[1] -77.77948


#the birth-death function from ape package, we can compare the output
birthdeath(tree)

#Estimation of Speciation and Extinction Rates
#with Birth-Death Models

#Phylogenetic tree: tree 
#Number of tips: 61 
#Deviance: 155.559 
#Log-likelihood: -77.77948 
#Parameter estimates:
#  d / b = 0   StdErr = 0 
#b - d = 0.02973589   StdErr = 0.002746608 
#(b: speciation rate, d: extinction rate)
#Profile likelihood 95% confidence intervals:
#  d / b: [-0.9543733, 0.2203064]
#b - d: [0.02278106, 0.03798103]

