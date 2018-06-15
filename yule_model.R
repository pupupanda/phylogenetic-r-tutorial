library(ape)
library(BioGeoBEARS)

##Pure Birth model
#This model assum that each species has a constant probability to speciation and no extinction

#Read in data
phy <- read.tree("elopomorph.tre")
#the yule function from ape package
yule(phy)
#$lambda
#[1] 0.02973639

#$se
#[1] 0.003871348

#$loglik
#[1] -77.77948

# birth rate (speciation rate)
lambda =1

log_likelihood <-function (lambda, phy) 
{



  
  #Total length of the tree
  X <- sum(phy$edge.length)

  
  #Number of internal nodes in the tree
  nb.node <- phy$Nnode 
  
  

  
  #log_likelihood
  
  loglik <- -lambda * X + lfactorial(phy$Nnode) + (nb.node - 1) * log(lambda)
 
  print(paste("lambda=",lambda," loglik=",loglik))
  return(loglik)
}

optim(par=c(lambda), fn = log_likelihood, lower = 0.0000001, control = list(fnscale = -1),method=c("L-BFGS-B"), phy=phy)

#$par
[1] 0.02974368

#$value
[1] -77.77948




