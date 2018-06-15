##This code is for a 3 state model which included 
#1. all equal rate model
#2. all different rate model
#3. Symmetric model
#4. some different rates



library(ape)
library(phytools)


#read in data
X<-read.csv("elopomorph.csv",row.names=1)
feed.mode<-setNames(X[,1],rownames(X))
feed.mode

feed.mode = as.factor(feed.mode) # NJM (may not be necessary but just in case)
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
nl <- nlevels(feed.mode)+1
nl
#
lvls <- levels(feed.mode)
lvls
#
feed.mode.in <- as.integer(feed.mode)
feed.mode.in

# NJM keep the names
names(feed.mode.in) = names(feed.mode)
feed.mode.in

#rate matrix
Q <- matrix(0, nl, nl,nl)
Q
eel.tree <- reorder(eel.tree, "postorder")
#ancestor
e1 <- eel.tree$edge [,1]
e1
#decendent
e2 <- eel.tree$edge [,2]
e2
EL <- eel.tree$edge.length
EL
TIPS <- 1:nb.tip
TIPS

# set up a 3 by n matrix
liks <- matrix(0, nb.tip + nb.node , nl)
head(liks)

#> head(liks)
#[,1] [,2] [,3]
#[1,]    0    0    0
#[2,]    0    0    0
#[3,]    0    0    0
#[4,]    0    0    0
#[5,]    0    0    0
#[6,]    0    0    0

# sort the tip data, so the names are in the same order as 
# the tip labels in the tree
order_of_tips_in_tree = match(x=eel.tree$tip.label, table=names(feed.mode.in))
order_of_tips_in_tree

# reorder accordingly
feed.mode.in = feed.mode.in[order_of_tips_in_tree]


# input states on the tip into the matrix
liks[cbind(TIPS, feed.mode.in)] <- 1
head(liks)

# set an arbitrary initial value for each rate
r1 = 0.1
r2 = 0.12
r3 = 0.15
r4 = 0.18
r5 = 0.09
r6 = 0.16


##Same rate model

#There is only 1 single rate
params_Same = c(r1)

log_likelihood_Same <- function (params_Same){
  
  comp <- numeric(nb.tip+nb.node-1)
  r1 = params_Same[1]
 
  
  
  # 3x3 matrix rate matrix
  Q = matrix(data=0, nrow=3, ncol=3)
  
  # assign different initial value for the rate matrix
  Q[1,2] = r1
  Q[2,1] = r1
  Q[1,3] = r1
  Q[3,1] = r1
  Q[2,3] = r1
  Q[3,2] = r1
  
  #the rate matrix, the rate from some state to any state is the same 
  Q
  #     [,1] [,2] [,3]
  #[1,]  0.0  0.1  0.1
  #[2,]  0.1  0.0  0.1
  #[3,]  0.1  0.1  0.0
  
  diag(Q) <- -rowSums(Q)
  
  #passing down from each tip
  for (i in seq(from = 1, by = 2, length.out = nb.node)) {
    
    j <- i + 1L
    #ancestor
    anc <- e1[i]
    #because this is a binary tree, each parent can have two possible children
    #decendent 1 
    des1 <- e2[i]
    #decendent 2
    des2 <- e2[j]
    
    #multiply the rate matrix with the transition matrix
    v.l <- expm(Q * EL[i]) %*% liks[des1, ]
    v.r <- expm(Q * EL[j]) %*% liks[des2, ]
    v <- v.l * v.r
    #print(v)
    comp[anc] <- sum(v)
    liks[anc, ] <- v/comp[anc]
  }
  # Lnl <- -liks[anc, ]
  node_lnLs = log(comp[-TIPS]) 
  lnL = sum(node_lnLs)
  #print (lnL)
  cat(params, lnL, "\n")
  return(lnL)
  
}


# All different rate


params_Diff = c(r1,r2,r3,r4,r5,r6)

log_likelihood_Diff <- function (params_Diff){
  
  comp <- numeric(nb.tip+nb.node-1)
  r1 = params_Diff[1]
  r2 = params_Diff[2]
  r3 = params_Diff[3]
  r4 = params_Diff[4]
  r5 = params_Diff[5]
  r6 = params_Diff[6]
  
  # 3x3 matrix rate matrix
  Q = matrix(data=0, nrow=3, ncol=3)
  
  # assign different initial value for the rate matrix
  Q[1,2] = r1
  Q[2,1] = r2
  Q[1,3] = r3
  Q[3,1] = r4
  Q[2,3] = r5
  Q[3,2] = r6
  
  #the rate matrix. rate from any state to any state is different
  Q
  #     [,1] [,2] [,3]
  #[1,]  0.00 0.10 0.30
  #[2,]  0.20 0.00 0.25
  #[3,]  0.15 0.35 0.00
  
  diag(Q) <- -rowSums(Q)
  
  #passing down from each tip
  for (i in seq(from = 1, by = 2, length.out = nb.node)) {
    
    j <- i + 1L
    #ancestor
    anc <- e1[i]
    #because this is a binary tree, each parent can have two possible children
    #decendent 1 
    des1 <- e2[i]
    #decendent 2
    des2 <- e2[j]
    
    #multiply the rate matrix with the transition matrix
    v.l <- expm(Q * EL[i]) %*% liks[des1, ]
    v.r <- expm(Q * EL[j]) %*% liks[des2, ]
    v <- v.l * v.r
    #print(v)
    comp[anc] <- sum(v)
    liks[anc, ] <- v/comp[anc]
  }
  # Lnl <- -liks[anc, ]
  node_lnLs = log(comp[-TIPS]) 
  lnL = sum(node_lnLs)
  #print (lnL)
  cat(params, lnL, "\n")
  return(lnL)
  
}



# Symmetric model:
params_Symmetric = c(r1, r2, r3)

log_likelihood_symmetric <- function (params_Symmetric){
  
  comp <- numeric(nb.tip+nb.node-1)
  r1 = params_Symmetric[1]
  r2 = params_Symmetric[2]
  r3 = params_Symmetric[3] 
  
  
  # 3x3 matrix rate matrix
  Q = matrix(data=0, nrow=3, ncol=3)

  # assign different initial value for the rate matrix
  Q[1,2] = r1
  Q[2,1] = r1
  Q[1,3] = r2
  Q[3,1] = r2
  Q[2,3] = r3
  Q[3,2] = r3
  
  #the rate matrix is symmetrical, i.e the rate from state a any state b is the same from any state b to a
  Q
  #     [,1] [,2] [,3]
  #[1,]  0.0  0.1  0.2
  #[2,]  0.1  0.0  0.3
  #[3,]  0.2  0.3  0.0
  
  diag(Q) <- -rowSums(Q)
  
  #passing down from each tip
  for (i in seq(from = 1, by = 2, length.out = nb.node)) {
    
    j <- i + 1L
    #ancestor
    anc <- e1[i]
    #because this is a binary tree, each parent can have two possible children
    #decendent 1 
    des1 <- e2[i]
    #decendent 2
    des2 <- e2[j]
    
    #multiply the rate matrix with the transition matrix
    v.l <- expm(Q * EL[i]) %*% liks[des1, ]
    v.r <- expm(Q * EL[j]) %*% liks[des2, ]
    v <- v.l * v.r
    #print(v)
    comp[anc] <- sum(v)
    liks[anc, ] <- v/comp[anc]
  }
  # Lnl <- -liks[anc, ]
  node_lnLs = log(comp[-TIPS]) 
  lnL = sum(node_lnLs)
  #print (lnL)
  cat(params, lnL, "\n")
  return(lnL)
  
}



# some different rates

params_partial = c(r1, r2, r3, r4)

log_likelihood_partial <- function (params_partial){
  
  comp <- numeric(nb.tip+nb.node-1)
  r1 = params_partial[1]
  r2 = params_partial[2]
  r3 = params_partial[3] 
  r4 = params_partial[4]  
  
  
  # 3x3 matrix rate matrix
  Q = matrix(data=0, nrow=3, ncol=3)
  
  # assign different initial value for the rate matrix
  Q[1,2] = r1
  Q[2,1] = r4
  Q[1,3] = r2
  Q[3,1] = r4
  Q[2,3] = r3
  Q[3,2] = r4
  
  #the rate matrix is symmetrical, i.e the rate from state a any state b is the same from any state b to a
  Q
  #     [,1] [,2] [,3]
  #[1,]  0.0  0.1  0.2
  #[2,]  0.1  0.0  0.3
  #[3,]  0.2  0.3  0.0
  
  diag(Q) <- -rowSums(Q)
  
  #passing down from each tip
  for (i in seq(from = 1, by = 2, length.out = nb.node)) {
    
    j <- i + 1L
    #ancestor
    anc <- e1[i]
    #because this is a binary tree, each parent can have two possible children
    #decendent 1 
    des1 <- e2[i]
    #decendent 2
    des2 <- e2[j]
    
    #multiply the rate matrix with the transition matrix
    v.l <- expm(Q * EL[i]) %*% liks[des1, ]
    v.r <- expm(Q * EL[j]) %*% liks[des2, ]
    v <- v.l * v.r
    #print(v)
    comp[anc] <- sum(v)
    liks[anc, ] <- v/comp[anc]
  }
  # Lnl <- -liks[anc, ]
  node_lnLs = log(comp[-TIPS]) 
  lnL = sum(node_lnLs)
  #print (lnL)
  cat(params, lnL, "\n")
  return(lnL)
  
}

### Maximum likelihood for each model

##Same rate

params_Symmetric = c(r1)
log_likelihood_Same(params_Same=params_Same)

optim(par=c(r1), fn=log_likelihood_Same, method=c("L-BFGS-B"),lower = 0.00000000001, control = list(fnscale = -1))
#output of the maximum log_likelihood

#$par
#[1] 0.00543165
#$value
#[1] -48.73838



##All different rate
params_Diff = c(r1,r2,r3,r4,r5,r6)
log_likelihood_Diff(params_Diff=params_Diff)

optim(par=c(r1,r2,r3,r4,r5,r6), fn=log_likelihood_Diff, method=c("L-BFGS-B"),lower = 0.00000000000000001, control = list(fnscale = -1))

#output of the maximum log_likelihood

#$par
#[1] 1.556471e-02 1.734934e-02 1.000000e-17 9.601901e-01 1.000000e-17 9.697112e-01
#$value
#[1] -35.90021


#Symmetric
params_Symmetric = c(r1, r2,r3)
log_likelihood_symmetric(params_Symmetric=params_Symmetric)

optim(par=c(r1,r2,r3), fn=log_likelihood_symmetric, method=c("L-BFGS-B"),lower = 0.00000000001, control = list(fnscale = -1))

#output of the maximum log_likelihood

#$par
#[1] 1.585683e-02 1.000000e-11 1.000000e-11
#$value
#[1] -36.33993




# some different
params_partial = c(r1, r2,r3 ,r4)
log_likelihood_partial(params_partial=params_partial)

optim(par=c(r1,r2,r3, r4), fn=log_likelihood_partial, method=c("L-BFGS-B"),lower = 0.00000000001, control = list(fnscale = -1))

#output of the maximum log_likelihood
#$par
#[1] 1.603536e-02 1.000000e-11 1.000000e-11 1.859640e-02

#$value
#[1] -36.0235
