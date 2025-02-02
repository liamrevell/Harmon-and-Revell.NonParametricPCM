## load phylogenetics libraries
library(geiger)
library(phytools)

## to parallelize
library(parallel)
library(doParallel)
library(foreach)

## type I error analysis

## set seed
set.seed(99)

## set up cluster
ncores<-max(10,detectCores())
mc<-makeCluster(ncores,type="PSOCK")
registerDoParallel(cl=mc)

## first simulate under uncorrelated Brownian motion
## for two traits

## this is a test of type I error under the null 
## for standard contrasts

## number of taxa
ntaxa<-seq(10,100,by=10)

## number of simulations
ns<-1000

## simulate trees (we'll re-utilize these)
sim_trees<-foreach(i=1:length(ntaxa))%dopar%{
  phytools::pbtree(n=ntaxa[i],scale=1,nsim=ns)
}

## simulate data
sim_bm<-foreach(i=1:length(ntaxa))%dopar%{
  lapply(sim_trees[[i]],phytools::fastBM,nsim=2)
}

## create matrices for results
pStandard<-pContrasts<-pContrastsPerm<-pRank<-
  pSign<-tContrasts<-tContrastsPerm<-tStandard<-
  tRank<-tSign<-matrix(nrow=ns, ncol=length(ntaxa),
    dimnames=list(1:ns,ntaxa))

## iterate over ntaxa
for(i in 1:length(ntaxa)){
  r0<-r1<-r2<-r3<-r4<-
    vector(mode="list",length=ns)
  for(j in 1:ns){
    t<-sim_trees[[i]][[j]]
    x<-sim_bm[[i]][[j]][,1]
    y<-sim_bm[[i]][[j]][,2]
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
  ## populate matrices with test statistics
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

## second simulate under uncorrelated Ornstein-Uhlenbeck
## for two traits with varied alpha levels

## this is a test of type I error under the null 
## for standard contrasts

## levels of alpha
alpha<-c(0.1,0.2,0.5,1,2,5,10,20,50,100)

## simulate data
foo<-function(phy,alpha){
  ou<-geiger:::rescale.phylo(phy,model="OU",alpha=alpha)
  phytools::fastBM(ou,nsim=2)
}

nn<-which(ntaxa==50)

sim_ou<-foreach(i=1:length(alpha))%dopar%{
  lapply(sim_trees[[nn]],foo,alpha=alpha[i])
}

## create matrices for results
pStandard<-pContrasts<-pContrastsPerm<-pRank<-
  pSign<-tContrasts<-tContrastsPerm<-tStandard<-
  tRank<-tSign<-matrix(nrow=ns, ncol=length(alpha),
    dimnames=list(1:ns,alpha))

## iterate over alpha
for(i in 1:length(alpha)){
  r0<-r1<-r2<-r3<-r4<-
    vector(mode="list",length=ns)
  for(j in 1:ns){
    t<-sim_trees[[nn]][[j]]
    x<-sim_ou[[i]][[j]][,1]
    y<-sim_ou[[i]][[j]][,2]
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
  ## populate matrices with test statistics
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
  cat(paste("Done with alpha = ",alpha[i]," ...\n",sep=""))
}

countSignif<-function(p) {
  sum(p<0.05)/length(p)
}
t2Standard<-apply(pStandard, 2, countSignif)
t2Contrasts<-apply(pContrasts, 2, countSignif)
t2ContrastsPerm<-apply(pContrastsPerm, 2, countSignif)
t2Sign<-apply(pSign, 2, countSignif)
t2Rank<-apply(pRank, 2, countSignif)

## subdivide plotting area
par(mfrow=c(1,2))

## create plot panel (a)
x<-ntaxa
plot(x, t1Standard,  type="l", ylim=c(0, 0.6), lwd=2, 
  xlab="number of taxa", ylab="type I error rate",
  las=1,bty="n")
lines(c(1, 100), c(0.05, 0.05), col="red")
lines(x, t1Contrasts, lwd=2, col="black", lty=2)
lines(x, t1ContrastsPerm, lwd=2, col="black", lty=3)
lines(x, t1Rank, lwd=2, col="grey")
lines(x, t1Sign, lwd=2, col="grey", lty=2)
grid()
legend("topleft", lwd=2, col=c("black", "black", "black", 
  "grey", "grey"), legend=c("Standard", "IC", "IC Perm", 
    "IC rank test", "IC sign test"), lty=c(1, 2,3, 1, 2))
mtext("a) type I error rate under Brownian motion",
  line=2,adj=0)

## create plot panel (b)
x<-alpha
plot(x, t2Standard,  type="l", ylim=c(0, 0.6), lwd=2, 
  xlab=expression(paste(alpha," of OU model")), 
  ylab="type I error rate",
  las=1,bty="n",log="x")
lines(c(0.001, 101), c(0.05, 0.05), col="red")
lines(x, t2Contrasts, lwd=2, col="black", lty=2)
lines(x, t2ContrastsPerm, lwd=2, col="black", lty=3)
lines(x, t2Rank, lwd=2, col="grey")
lines(x, t2Sign, lwd=2, col="grey", lty=2)
grid()
legend("topleft", lwd=2, col=c("black", "black", "black", 
  "grey", "grey"), legend=c("Standard", "IC", "IC Perm", 
    "IC rank test", "IC sign test"), lty=c(1, 2,3, 1, 2))
mtext("b) type I error rate under Ornstein-Uhlenbeck",
  line=2,adj=0)

stopCluster(mc)

## save workspace & clean up (except trees)
save.image(file="typeI-error-analysis.Rdata")

rm(list=ls()[-which(ls()%in%c("sim_trees","picRegression","ns"))])

########################################################

## power analysis

## set seed
set.seed(99)

## set up cluster
ncores<-max(10,detectCores())
mc<-makeCluster(ncores,type="PSOCK")
registerDoParallel(cl=mc)

## third simulate correlated Brownian motion
## for two traits

## this is a test of power

## number of taxa
ntaxa<-seq(10,100,by=10)

## simulate data
cor_mat<-matrix(c(1,0.5,0.5,1),2,2)
sim_corbm<-foreach(i=1:length(ntaxa))%dopar%{
  lapply(sim_trees[[i]],phytools::sim.corrs,vcv=cor_mat)
}

## create matrices for results
pStandard<-pContrasts<-pContrastsPerm<-pRank<-
  pSign<-tContrasts<-tContrastsPerm<-tStandard<-
  tRank<-tSign<-matrix(nrow=ns, ncol=length(ntaxa),
    dimnames=list(1:ns,ntaxa))

## iterate over ntaxa
for(i in 1:length(ntaxa)){
  r0<-r1<-r2<-r3<-r4<-
    vector(mode="list",length=ns)
  for(j in 1:ns){
    t<-sim_trees[[i]][[j]]
    x<-sim_corbm[[i]][[j]][,1]
    y<-sim_corbm[[i]][[j]][,2]
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
  ## populate matrices with test statistics
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
  ## give feedback
  cat(paste("Done with N = ",ntaxa[i]," ...\n",sep=""))
}

countSignif<-function(p) {
  sum(p<0.05)/length(p)
}
t3Standard<-apply(pStandard, 2, countSignif)
t3Contrasts<-apply(pContrasts, 2, countSignif)
t3ContrastsPerm<-apply(pContrastsPerm, 2, countSignif)
t3Sign<-apply(pSign, 2, countSignif)
t3Rank<-apply(pRank, 2, countSignif)

## fourth simulate correlated evolution on OU tree
## for two traits

## this is a test of power

## levels of alpha
alpha<-c(0.1,0.2,0.5,1,2,5,10,20,50,100)

## simulate data
cor_mat<-matrix(c(1,0.5,0.5,1),2,2)

foo<-function(phy,alpha){
  ou<-geiger:::rescale.phylo(phy,model="OU",alpha=alpha)
  phytools::sim.corrs(ou,vcv=cor_mat)
}

nn<-which(ntaxa==50)

sim_corou<-foreach(i=1:length(alpha))%dopar%{
  lapply(sim_trees[[nn]],foo,alpha=alpha[i])
}

## create matrices for results
pStandard<-pContrasts<-pContrastsPerm<-pRank<-
  pSign<-tContrasts<-tContrastsPerm<-tStandard<-
  tRank<-tSign<-matrix(nrow=ns, ncol=length(alpha),
    dimnames=list(1:ns,alpha))

## iterate over alpha
for(i in 1:length(alpha)){
  r0<-r1<-r2<-r3<-r4<-
    vector(mode="list",length=ns)
  for(j in 1:ns){
    t<-sim_trees[[nn]][[j]]
    x<-sim_corou[[i]][[j]][,1]
    y<-sim_corou[[i]][[j]][,2]
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
  ## populate matrices with test statistics
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
  ## give feedback
  cat(paste("Done with alpha = ",alpha[i]," ...\n",sep=""))
}

countSignif<-function(p) {
  sum(p<0.05)/length(p)
}
t4Standard<-apply(pStandard, 2, countSignif)
t4Contrasts<-apply(pContrasts, 2, countSignif)
t4ContrastsPerm<-apply(pContrastsPerm, 2, countSignif)
t4Sign<-apply(pSign, 2, countSignif)
t4Rank<-apply(pRank, 2, countSignif)

## subdivide plotting area
par(mfrow=c(1,2))

## create panel (a) of plot
x<-ntaxa
plot(x, t3Standard,  type="l", ylim=c(0, 1), lwd=2, 
  xlab="number of taxa", ylab="power to reject null",
  las=1,bty="n")
lines(x, t3Contrasts, lwd=2, col="black", lty=2)
lines(x, t3ContrastsPerm, lwd=2, col="black", lty=3)
lines(x, t3Rank, lwd=2, col="grey")
lines(x, t3Sign, lwd=2, col="grey", lty=2)
grid()
legend("bottomright", lwd=2, col=c("black", "black", "black", 
  "grey", "grey"), legend=c("Standard", "IC", "IC Perm", 
    "IC rank test", "IC sign test"), lty=c(1, 2,3, 1, 2))
mtext("a) power under Brownian motion with r = 0.5",
  line=2,adj=0)

## create panel (b) of plot
x<-alpha
plot(x, t4Standard,  type="l", ylim=c(0, 1), lwd=2, 
  xlab=expression(paste(alpha," of OU model")), 
  ylab="power to reject null",
  las=1,bty="n",log="x")
lines(x, t4Contrasts, lwd=2, col="black", lty=2)
lines(x, t4ContrastsPerm, lwd=2, col="black", lty=3)
lines(x, t4Rank, lwd=2, col="grey")
lines(x, t4Sign, lwd=2, col="grey", lty=2)
grid()

legend("bottomright", lwd=2, col=c("black", "black", "black", 
  "grey", "grey"), legend=c("Standard", "IC", "IC Perm", 
    "IC rank test", "IC sign test"), lty=c(1, 2,3, 1, 2))

mtext("b) power under Ornstein-Uhlenbeck with r = 0.5",
  line=2,adj=0)

## stop cluster
stopCluster(mc)

## save workspace
save.image(file="power-analysis.Rdata")

rm(list=ls()[-which(ls()%in%c("sim_trees","picRegression","ns"))])

########################################################

## effect size

## load workspace
load(file="typeI-error-analysis.Rdata")
rm(list=ls()[-which(ls()%in%c("sim_trees","picRegression","ns"))])

ns<-length(sim_trees[[1]])

## set seed
set.seed(99)

## set up cluster
ncores<-max(10,detectCores())
mc<-makeCluster(ncores,type="PSOCK")
registerDoParallel(cl=mc)

## number of taxa
ntaxa<-seq(10,100,by=10)

## correlations
r_sim<-seq(0,1,by=0.1)

## set trees to use
nn<-which(ntaxa==50)

## BM function
foo<-function(phy,r){
  cor_mat<-matrix(c(1,r,r,1),2,2)
  if(r<1) return(phytools::sim.corrs(phy,vcv=cor_mat))
  else {
    tmp<-phytools::fastBM(phy)
    return(cbind(tmp,tmp))
  }
}

## simulate under BM
sim_corbm<-foreach(i=1:length(r_sim))%dopar%{
  lapply(sim_trees[[nn]],foo,r=r_sim[i])
}

## OU function
foo<-function(phy,r,alpha=10){
  cor_mat<-matrix(c(1,r,r,1),2,2)
  ou<-geiger:::rescale.phylo(phy,model="OU",alpha=alpha)
  if(r<1) return(phytools::sim.corrs(ou,vcv=cor_mat))
  else {
    tmp<-phytools::fastBM(ou)
    return(cbind(tmp,tmp))
  }
}

## simulate under OU
sim_corou<-foreach(i=1:length(r_sim))%dopar%{
  lapply(sim_trees[[nn]],foo,r=r_sim[i])
}

## create matrices for results
pStandard<-pContrasts<-pContrastsPerm<-pRank<-
  pSign<-tContrasts<-tContrastsPerm<-tStandard<-
  tRank<-tSign<-matrix(nrow=ns, ncol=length(r_sim),
    dimnames=list(1:ns,r_sim))

## first BM

## iterate over r_sim
for(i in 1:length(r_sim)){
  r0<-r1<-r2<-r3<-r4<-
    vector(mode="list",length=ns)
  for(j in 1:ns){
    t<-sim_trees[[nn]][[j]]
    x<-sim_corbm[[i]][[j]][,1]
    y<-sim_corbm[[i]][[j]][,2]
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
  ## populate matrices with test statistics
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
  ## give feedback
  cat(paste("Done with r = ",r_sim[i]," ...\n",sep=""))
}

countSignif<-function(p) {
  sum(p<0.05)/length(p)
}
t5Standard<-apply(pStandard, 2, countSignif)
t5Contrasts<-apply(pContrasts, 2, countSignif)
t5ContrastsPerm<-apply(pContrastsPerm, 2, countSignif)
t5Sign<-apply(pSign, 2, countSignif)
t5Rank<-apply(pRank, 2, countSignif)

## now OU

## create matrices for results
pStandard<-pContrasts<-pContrastsPerm<-pRank<-
  pSign<-tContrasts<-tContrastsPerm<-tStandard<-
  tRank<-tSign<-matrix(nrow=ns, ncol=length(r_sim),
    dimnames=list(1:ns,r_sim))

## iterate over r_sim
for(i in 1:length(r_sim)){
  r0<-r1<-r2<-r3<-r4<-
    vector(mode="list",length=ns)
  for(j in 1:ns){
    t<-sim_trees[[nn]][[j]]
    x<-sim_corou[[i]][[j]][,1]
    y<-sim_corou[[i]][[j]][,2]
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
  ## populate matrices with test statistics
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
  ## give feedback
  cat(paste("Done with r = ",r_sim[i]," ...\n",sep=""))
}

countSignif<-function(p) {
  sum(p<0.05)/length(p)
}
t6Standard<-apply(pStandard, 2, countSignif)
t6Contrasts<-apply(pContrasts, 2, countSignif)
t6ContrastsPerm<-apply(pContrastsPerm, 2, countSignif)
t6Sign<-apply(pSign, 2, countSignif)
t6Rank<-apply(pRank, 2, countSignif)


## subdivide plotting area
par(mfrow=c(1,2))

## create panel (a) of plot
x<-r_sim
plot(x, t5Standard,  type="l", xlim=c(0,1), 
  ylim=c(0, 1), lwd=2, 
  xlab="r", ylab="power to reject null",
  las=1,bty="n")
lines(x, t5Contrasts, lwd=2, col="black", lty=2)
lines(x, t5ContrastsPerm, lwd=2, col="black", lty=3)
lines(x, t5Rank, lwd=2, col="grey")
lines(x, t5Sign, lwd=2, col="grey", lty=2)
grid()
legend("bottomright", lwd=2, col=c("black", "black", "black", 
  "grey", "grey"), legend=c("Standard", "IC", "IC Perm", 
    "IC rank test", "IC sign test"), lty=c(1, 2,3, 1, 2))
mtext("a) power under Brownian motion with various r",
  line=2,adj=0)

## create panel (b) of plot
x<-r_sim
plot(x, t6Standard,  type="l", xlim=c(0,1), 
  ylim=c(0, 1), lwd=2, 
  xlab="r", ylab="power to reject null",
  las=1,bty="n")
lines(x, t6Contrasts, lwd=2, col="black", lty=2)
lines(x, t6ContrastsPerm, lwd=2, col="black", lty=3)
lines(x, t6Rank, lwd=2, col="grey")
lines(x, t6Sign, lwd=2, col="grey", lty=2)
grid()
legend("bottomright", lwd=2, col=c("black", "black", "black", 
  "grey", "grey"), legend=c("Standard", "IC", "IC Perm", 
    "IC rank test", "IC sign test"), lty=c(1, 2,3, 1, 2))
mtext("b) power under OU with various r",
  line=2,adj=0)

## stop cluster
stopCluster(mc)

## save workspace
save.image(file="effect-size.Rdata")


