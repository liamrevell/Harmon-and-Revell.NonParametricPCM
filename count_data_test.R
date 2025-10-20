# install.packages(c("ape","phytools"))  # if needed
library(ape)
library(phytools)
library(parallel)
library(doParallel)
# Simulate count data on a tree under a stepwise Mk model with two rates
simulate_mk_counts <- function(tree,
                               q_up, q_down,           # rates for +1 and -1
                               max_state = 50,         # upper cap for state-space
                               root_state = 0,         # starting count at the root
                               nsim = 1) {
  
  stopifnot(inherits(tree, "phylo"))
  stopifnot(root_state >= 0, root_state <= max_state)
  
  # Build ordered-state Q matrix: states 0..max_state
  states <- 0:max_state
  lab    <- as.character(states)
  Q <- matrix(0, length(states), length(states),
              dimnames = list(lab, lab))
  
  for (i in seq_along(states)) {
    s <- states[i]
    if (s < max_state) Q[i, i + 1] <- q_up
    if (s > 0)         Q[i, i - 1] <- q_down
  }
  diag(Q) <- -rowSums(Q)
  Q<-t(Q)
  
  one_sim <- function() {
    sm <- phytools::sim.history(tree, Q, anc = as.character(root_state))
    tip_states  <- phytools::getStates(sm, "tips")

    # Convert to integer counts
    tip_counts  <- setNames(as.integer(tip_states),  names(tip_states))

    tip_counts
  }
  
  if (nsim == 1) return(one_sim())
  sapply(seq_len(nsim), function(i) one_sim())
}



picRegression <- function (tree, x, y, method="standard", sigTest="permutation",
                           ...) 
{
  ## original version was missing nperm argument so I added it (LR)
  if(hasArg(nperm)) nperm<-list(...)$nperm
  else nperm<-100
  if (!inherits(tree, "phylo")) 
    stop("tree should be object of class \"phylo\".")
  
  if (method %in% c("standard", "sign", "rank") == FALSE) {
    cat("  Invalid model. Setting model=\"standard\"\n\n")
    model <- "standard"
  }
  
  if (sigTest %in% c("analytic", "permutation") == FALSE) {
    cat("  Invalid model. Setting model=\"permutation\"\n\n")
    sigTest <- "permutation"
  }
  
  if((method=="sign" | method=="rank")&sigTest=="analytic") {
    cat("  No analytic p-value method exists for method ", method, ". Setting model=\"permutation\"\n\n")
    sigTest <- "permutation"
  }
  
  if(method=="standard") {
    icx<-ape::pic(x, tree)
    icy<-ape::pic(y, tree)
    
    if(sigTest=="analytic") {
      res<-summary(lm(icy~icx+0))
      testStat<-res$coefficients[1,3]
      pVal<-res$coefficients[1,4]
    }
    if(sigTest=="permutation") {
      dicx<-c(icx, -icx)	
      dicy<-c(icy, -icy)
      res<-summary(lm(dicy~dicx))
      testStat<-res$coefficients[2,3]
      nullDist<-numeric(nperm)
      for(i in 1:nperm) {
        picx<-sample(icx)
        dpicx<-c(picx, -picx)
        res<-summary(lm(dicy~dpicx))
        nullDist[i]<-res$coefficients[2,3]
      }
      
      pValueHigh<-2*(sum(nullDist >= testStat)+1)/(nperm+1)
      pValueLow<-2*(sum(nullDist <= testStat)+1)/(nperm+1)
      
      pVal<-min(pValueHigh, pValueLow)
    }
    
  }
  
  if(method=="sign") {
    icx<-ape::pic(x, tree)
    icy<-ape::pic(y, tree)
    xPos<-icx>0
    yPos<-icy>0
    testStat<-sum(xPos==yPos & icx!=icy)
    pValueLow<-2*pbinom(testStat, sum(icx!=icy), prob=0.5)
    pValueHigh<-2*pbinom(length(icx)-testStat, sum(icx!=icy), prob=0.5)
    pVal<-min(pValueHigh, pValueLow)
  }
  
  if(method=="rank") {
    icx<-ape::pic(x, tree)
    icy<-ape::pic(y, tree)
    
    dicx<-c(icx, -icx)	
    dicy<-c(icy, -icy)
    
    rx<-rank(dicx, ties.method="average")
    ry<-rank(dicy, ties.method="average")
    res<-summary(lm(ry~rx))
    testStat<-res$coefficients[2,3]
    
    nullDist<-numeric(nperm)
    for(i in 1:nperm) {
      picx<-sample(icx)
      dpicx<-c(picx, -picx)
      prx<-rank(dpicx, ties.method="average")
      res<-summary(lm(ry~prx))
      nullDist[i]<-res$coefficients[2,3]
    }
    
    pValueHigh<-2*(sum(nullDist >= testStat)+1)/(nperm+1)
    pValueLow<-2*(sum(nullDist <= testStat)+1)/(nperm+1)
    
    pVal<-min(pValueHigh, pValueLow)
    
  }
  
  
  obj <- list(testStat=testStat, pVal=pVal)
  class(obj) <- "picRegression"
  obj
}


# ---------- EXAMPLE ----------
set.seed(123)
tr <- ape::rtree(80)

# Up-steps twice as likely as down-steps; lower bound at 0
sim1 <- simulate_mk_counts(tr, q_up = 0.6, q_down = 0.3,
                           max_state = 60, root_state = 0, nsim=2)

picRegression(tr, sim1[,1], sim1[,2], method="sign")

# repeat fig 2


ntaxa<-seq(10,100,by=10)
## number of simulations
ns<-1000

## simulate trees (we'll re-utilize these)
sim_trees<-foreach(i=1:length(ntaxa))%dopar%{
  phytools::pbtree(n=ntaxa[i],scale=1,nsim=ns)
}
r0<-r1<-r2<-r3<-r4<-
  vector(mode="list",length=ns)


## create matrices for results
pStandard<-pContrasts<-pContrastsPerm<-pRank<-
  pSign<-tContrasts<-tContrastsPerm<-tStandard<-
  tRank<-tSign<-matrix(nrow=ns, ncol=length(ntaxa),
                       dimnames=list(1:ns,ntaxa))

for(i in 1:length(sim_trees)) {
  for(j in 1:length(sim_trees[[i]])) {
    t<-sim_trees[[i]][[j]]
    while(1) {
      ss<-simulate_mk_counts(t, q_up=1, q_down=1, nsim=2)
      if(length(unique(ss[,1]))>1 && length(unique(ss[,2]))>1) break;
    }
    x<-ss[,1]
    y<-ss[,2]
    r0[[j]]<-summary(lm(y~x))
    r1[[j]]<-picRegression(t, x, y, method="standard", 
                           sigTest="analytic")
    r2[[j]]<-picRegression(t, x, y, method="standard", 
                           sigTest="permutation")
    r3[[j]]<-picRegression(t, x, y, method="sign", 
                           sigTest="permutation")
    r4[[j]]<-picRegression(t, x, y, method="rank", 
                           sigTest="permutation")
  }
  tStandard[,i]<-sapply(r0,function(x) x$coefficients[2,4])
  tContrasts[,i]<-sapply(r1,function(x) x$testStat)
  tContrastsPerm[,i]<-sapply(r2,function(x) x$testStat)
  tSign[,i]<-sapply(r3,function(x) x$testStat)
  tRank[,i]<-sapply(r4,function(x) x$testStat)
  ## populate matrices with P-values
  pStandard[,i]<-sapply(r0,function(x) x$coefficients[2,4])
  pContrasts[,i]<-sapply(r1,function(x) x$pVal)
  pContrastsPerm[,i]<-sapply(r2,function(x) x$pVal)
  pSign[,i]<-sapply(r3,function(x) x$pVal)
  pRank[,i]<-sapply(r4,function(x) x$pVal)
  ## print feedback
  cat(paste("Done with N = ",ntaxa[i]," ...\n",sep=""))
}
countSignif<-function(p) {
  sum(p<0.05)/length(p)
}
t1Standard<-apply(pStandard, 2, countSignif)
t1Contrasts<-apply(pContrasts, 2, countSignif)
t1ContrastsPerm<-apply(pContrastsPerm, 2, countSignif)
t1Sign<-apply(pSign, 2, countSignif)
t1Rank<-apply(pRank, 2, countSignif)

x<-ntaxa
plot(x, t1Standard,  type="l", ylim=c(0, 1), lwd=2, 
     xlab="number of taxa", ylab="type I error rate",
     las=1,bty="n")
lines(c(1, 110), c(0.05, 0.05), col="red")
lines(x, t1Contrasts, lwd=2, col="black", lty=2)
lines(x, t1ContrastsPerm, lwd=2, col="black", lty=3)
lines(x, t1Rank, lwd=2, col="grey")
lines(x, t1Sign, lwd=2, col="grey", lty=2)
grid()
legend("topleft", lwd=2,, cex=0.5, col=c("black", "black", "black", "grey", "grey"), legend=c("Standard", "IC", "IC Perm", 
                                                         "PRC" , "CST"), lty=c(1, 2,3, 1, 2), bg="white")
mtext("Type I error rate under count data",
      line=2,adj=0)
