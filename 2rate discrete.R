library(phytools)

#2 rate Markov model for a discrete binary charater

wd <- ("D:/2018 semster 1/Bioinf 702/assignment 4")
setwd(wd)
#read in data
X<-read.csv("elopomorph.csv",row.names=1)
feed.mode<-setNames(X[,1],rownames(X))
feed.mode
#read in tree
eel.tree<-read.tree("elopomorph.tre")
eel.tree
eel.tree <- reorder(eel.tree, "postorder")
#number of species
nb.tip <- length(eel.tree$tip.label)
nb.tip
#number of nodes
nb.node <- eel.tree$Nnode
nb.node
plot(eel.tree)


#number of states of the character
nl <- nlevels(feed.mode)
nl
#the level of characters
lvls <- levels(feed.mode)
lvls

#transform the charaters into integer
feed.mode.in <- as.integer(feed.mode)
feed.mode.in

#rate matrix
Q <- matrix(0, nl, nl)
Q

#recorder the tree in a decendent form
eel.tree <- reorder(eel.tree, "postorder")

#ancestor edge
e1 <- eel.tree$edge [,1]
e1
#decendent edge
e2 <- eel.tree$edge [,2]
e2

#length of each edges
EL <- eel.tree$edge.length
EL

#Tips of the tree
TIPS <- 1:nb.tip
TIPS

#likelihood matrix
liks <- matrix(0, nb.tip + nb.node , nl)
liks[cbind(TIPS, feed.mode.in)] <- 1


#intitial value for transition rate from suction to bite and from bite to suction
r1 = 0.5
r2 = 0.5

#set parameters for the log_likelihood function 
params = c(r1, r2)

log_likelihood <- function (params){
  
  comp <- numeric(nb.tip+nb.node-1)
  r1 = params[1]
  r2 = params[2]
  Q[1,2] <- r1
  Q[2,1] <- r2
  
  diag(Q) <- -rowSums(Q)
  
  for (i in seq(from = 1, by = 2, length.out = nb.node)) {
    j <- i + 1L
    #ancestor
    anc <- e1[i]
    #because this is a binary tree, each parent can have two possible children
    #decendent 1 
    des1 <- e2[i]
    #decendent 2
    des2 <- e2[j]
    
    #multiply the rate matrix and each possible descendants
    v.l <- expm(Q * EL[i]) %*% liks[des1, ]
    v.r <- expm(Q * EL[j]) %*% liks[des2, ]
    v <- v.l * v.r
   
    
    comp[anc] <- sum(v)
    liks[anc, ] <- v/comp[anc]
  }
  
  # Lnl <- -liks[anc, ]
  node_lnLs = log(comp[-TIPS]) 
  lnL = sum(node_lnLs)
  
  cat(params, lnL, "\n")
  return(lnL)
}

optim(par=c(r1,r2), fn=log_likelihood, method=c("L-BFGS-B"),lower = 0.001, control = list(fnscale = -1))



fitARD<-ace(feed.mode.in,eel.tree,model="ARD",type="discrete")
fitARD
