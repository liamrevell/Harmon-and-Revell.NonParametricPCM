# nonpar phytools style

picRegression <- function (tree, x, y, method="standard", sigTest="permutation") 
{
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
    icx<-pic(x, tree)
    icy<-pic(y, tree)
    
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
    icx<-pic(x, tree)
    icy<-pic(y, tree)
    xPos<-icx>0
    yPos<-icy>0
    testStat<-sum(xPos==yPos)
    pValueLow<-2*pbinom(testStat, length(icx), prob=0.5)
    pValueHigh<-2*pbinom(length(icx)-testStat, length(icx), prob=0.5)
    pVal<-min(pValueHigh, pValueLow)
  }
  
  if(method=="rank") {
    icx<-pic(x, tree)
    icy<-pic(y, tree)
    
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