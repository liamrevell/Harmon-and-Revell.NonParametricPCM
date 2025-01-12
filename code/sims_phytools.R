require(geiger)
require(phytools)
require(TreeSim)

ns<-100


pStandard<-matrix(nrow=ns, ncol=10)
pContrasts<-matrix(nrow=ns, ncol=10)
pContrastsPerm<-matrix(nrow=ns, ncol=10)
pRank<-matrix(nrow=ns, ncol=10)
pSign<-matrix(nrow=ns, ncol=10)

tContrasts<-matrix(nrow=ns, ncol=10)
tContrastsPerm<-matrix(nrow=ns, ncol=10)
tStandard<-matrix(nrow=ns, ncol=10)
tRank<-matrix(nrow=ns, ncol=10)
tSign<-matrix(nrow=ns, ncol=10)

par(mfrow=c(2,2))

## this is type I error when x & y are simulated under Brownian motion

mm<-cbind(c(1,0),c(0,1))
for(i in 1:10) {
  treeSet<-sim.bd.taxa(n=i*5, numbsim=ns, lambda=1, mu=0)
  simSet<-lapply(treeSet, function(x) sim.char(x, mm))
  for(j in 1:ns) {
    t<-treeSet[[j]]
    x<-simSet[[j]][,1,1]
    y<-simSet[[j]][,2,1]
    r0<-summary(lm(y~x))
    r1<-picRegression(t, x, y, method="standard", sigTest="analytic")
    r2<-picRegression(t, x, y, method="standard", sigTest="permutation")
    r3<-picRegression(t, x, y, method="sign", sigTest="permutation")
    r4<-picRegression(t, x, y, method="rank", sigTest="permutation")
    
    tStandard[j,i]<-r0$coefficients[2,3]
    tContrasts[j,i]<-r1$testStat
    tContrastsPerm[j,i]<-r2$testStat
    tSign[j,i]<-r3$testStat
    tRank[j,i]<-r4$testStat
    
    pStandard[j,i]<-r0$coefficients[2,4]
    pContrasts[j,i]<-r1$pVal
    pContrastsPerm[j,i]<-r2$pVal
    pSign[j,i]<-r3$pVal
    pRank[j,i]<-r4$pVal
    
  }
}

countSignif<-function(p) {
  sum(p<0.05)/length(p)
}
t1Standard<-apply(pStandard, 2, countSignif)
t1Contrasts<-apply(pContrasts, 2, countSignif)
t1ContrastsPerm<-apply(pContrastsPerm, 2, countSignif)
t1Sign<-apply(pSign, 2, countSignif)
t1Rank<-apply(pRank, 2, countSignif)

x<-1:10*5
plot(x, t1Standard,  type="l", ylim=c(0, 1), lwd=2, 
  xlab="Tree size", ylab="Type I error rate",
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

## now how about type I error when non-Brownian

mm<-cbind(c(1,0),c(0,1))
for(i in 1:10) {
  treeSet<-sim.bd.taxa(n=i*5, numbsim=ns, lambda=1, mu=0)
  simSet<-lapply(treeSet, function(x) sim.char(x, mm))
  for(j in 1:ns) {
    t<-treeSet[[j]]
    x<-simSet[[j]][,1,1]
    y<-simSet[[j]][,2,1]
    r0<-summary(lm(y~x))
    r1<-picRegression(t, x, y, method="standard", sigTest="analytic")
    r2<-picRegression(t, x, y, method="standard", sigTest="permutation")
    r3<-picRegression(t, x, y, method="sign", sigTest="permutation")
    r4<-picRegression(t, x, y, method="rank", sigTest="permutation")
    
    tStandard[j,i]<-r0$coefficients[2,3]
    tContrasts[j,i]<-r1$testStat
    tContrastsPerm[j,i]<-r2$testStat
    tSign[j,i]<-r3$testStat
    tRank[j,i]<-r4$testStat
    
    pStandard[j,i]<-r0$coefficients[2,4]
    pContrasts[j,i]<-r1$pVal
    pContrastsPerm[j,i]<-r2$pVal
    pSign[j,i]<-r3$pVal
    pRank[j,i]<-r4$pVal
    
  }
}

countSignif<-function(p) {
  sum(p<0.05)/length(p)
}
t1Standard<-apply(pStandard, 2, countSignif)
t1Contrasts<-apply(pContrasts, 2, countSignif)
t1ContrastsPerm<-apply(pContrastsPerm, 2, countSignif)
t1Sign<-apply(pSign, 2, countSignif)
t1Rank<-apply(pRank, 2, countSignif)

x<-1:10*5
plot(x, t1Standard,  type="l", ylim=c(0, 1), lwd=2, 
  xlab="Tree size", ylab="Type I error rate",
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








pStandard<-matrix(nrow=ns, ncol=10)
pContrasts<-matrix(nrow=ns, ncol=10)
pContrastsPerm<-matrix(nrow=ns, ncol=10)
pRank<-matrix(nrow=ns, ncol=10)
pSign<-matrix(nrow=ns, ncol=10)

tContrasts<-matrix(nrow=ns, ncol=10)
tContrastsPerm<-matrix(nrow=ns, ncol=10)
tStandard<-matrix(nrow=ns, ncol=10)
tRank<-matrix(nrow=ns, ncol=10)
tSign<-matrix(nrow=ns, ncol=10)


mm<-cbind(c(1,0.5),c(0.5,1))
for(i in 1:10) {
  treeSet<-sim.bd.taxa(n=i*5, numbsim=ns, lambda=1, mu=0)
  simSet<-lapply(treeSet, function(x) sim.char(x, mm))
  for(j in 1:ns) {
    t<-treeSet[[j]]
    x<-simSet[[j]][,1,1]
    y<-simSet[[j]][,2,1]
    r0<-summary(lm(y~x))
    r1<-picRegression(t, x, y, method="standard", sigTest="analytic")
    r2<-picRegression(t, x, y, method="standard", sigTest="permutation")
    r3<-picRegression(t, x, y, method="sign", sigTest="permutation")
    r4<-picRegression(t, x, y, method="rank", sigTest="permutation")
    
    tStandard[j,i]<-r0$coefficients[2,3]
    tContrasts[j,i]<-r1$testStat
    tContrastsPerm[j,i]<-r2$testStat
    tSign[j,i]<-r3$testStat
    tRank[j,i]<-r4$testStat
    
    pStandard[j,i]<-r0$coefficients[2,4]
    pContrasts[j,i]<-r1$pVal
    pContrastsPerm[j,i]<-r2$pVal
    pSign[j,i]<-r3$pVal
    pRank[j,i]<-r4$pVal
    
  }
}

t1Standard<-apply(pStandard, 2, countSignif)
t1Contrasts<-apply(pContrasts, 2, countSignif)
t1ContrastsPerm<-apply(pContrastsPerm, 2, countSignif)
t1Sign<-apply(pSign, 2, countSignif)
t1Rank<-apply(pRank, 2, countSignif)

x<-1:10*5
plot(x, t1Standard,  type="l", ylim=c(0, 1), lwd=2, 
  xlab="Tree size", ylab="Type 1 error rate")
lines(c(1, 100), c(0.05, 0.05), col="red")
lines(x, t1Contrasts, lwd=2, col="black", lty=2)
lines(x, t1ContrastsPerm, lwd=2, col="black", lty=3)
lines(x, t1Rank, lwd=2, col="grey")
lines(x, t1Sign, lwd=2, col="grey", lty=2)

legend("topleft", lwd=2, col=c("black", "black", "black", 
  "grey", "grey"), legend=c("Standard", "IC", "IC Perm", 
    "IC rank test", "IC sign test"), lty=c(1, 2,3, 1, 2))


