library(ape)
library(BioGeoBEARS)

##Pure Birth model
#This model assum that each species has a constant probability to speciation and no extinction

#Read in data
phy <- read.nexus("http://www.r-phylo.org/w/images/0/02/Geospiza.nex")
#the yule function from ape package
yule(phy)


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

yule = optim(par=c(lambda), fn = log_likelihood, lower = 0, control = list(fnscale = -1),method=c("L-BFGS-B"), phy=phy)






