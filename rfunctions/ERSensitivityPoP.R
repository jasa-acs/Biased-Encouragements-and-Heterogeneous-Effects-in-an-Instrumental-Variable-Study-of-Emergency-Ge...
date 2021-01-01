########################
#ERSensitivityPoP
########################
#Function to perform a Studentized Sensitivity analysis on the effect ratio in a paired observational study
#Incorporates Covariates to improve Standard Errors, using Pairs of Pairs Standard Error

#######
#Input
######

#Y: vector of encouraged-minus-not paired differences
#D: Vector of encouraged-minus-not treatment levels
#null: the value of the attributable effect
#alpha: desired Type I error rate
#alternative: can be "less", "greater", or "two.sided"
#Gamma: vector of values for Gamma at which to perform the sensitivity analysis
#Changepoint: if true, function returns the maximal Gamma for which the test rejects at level alpha
#SensitivityInterval: if true, function returns (100-alpha) sensitivity intervals. They will be one-sided if the alternative is less than or greater than, and two-sided if the alternative is two.sided



#######
#Output
#######
#Gamma: vector of Gammas for which the sensitivity analysis was performed
#pval: p-values for each value of Gamma
#GammaPval: matrix combining Gamma and pval
#Changepoint: maximal Gamma for which the test rejected at level alpha
#SensitivityInterval: upper and lower bounds for 100(1-alpha) sensitivity intervals for each value of Gamma

ERSensitivityPoP = function(Y, D, index = NULL, null = 0, alpha = 0.05, alternative = "greater", Gamma=1, nperm = 50000, Changepoint = T, SensitivityInterval = T)
{
  if(any(Gamma < 1))
  {
    stop("Values for Gamma must be >= 1")
  }
  if(alternative!="less" & alternative!= "greater" & alternative != "two.sided")
  {
    stop("Values for alternative are `less', `greater', or `two.sided'")
  }
  if(length(null) > 1)
  {
    stop("Value under the null must be a scalar")
  }
  if(alpha < 0 | alpha > 0.5)
  {
    stop("alpha must be between 0 and 0.5")
  }
  Ytrue = Y
  Dtrue = D
  alphatrue = alpha
  I = length(Y)
  Adjust = Y - null*D
  if(alternative == "less")
  {
    Adjust = -Adjust
  }
  if(alternative == "two.sided")
  {
    alpha = alphatrue/2
    if(mean(Adjust) < 0) 
    {
      Adjust = -Adjust
    }
  }
  pval = rep(0, length(Gamma))
  if(is.null(index)){
    for(i in 1:length(Gamma))
    {
      D2 = (Adjust) - (Gamma[i]-1)/(1+Gamma[i])*abs(Adjust)
      obs = mean(D2)/(sd(D2)/sqrt(I))
      Adjmat = matrix(abs(Adjust), I, nperm)
      Zmat = matrix(runif(I*nperm) <Gamma[i]/(1+Gamma[i]), I, nperm)
      Dmat = (2*Zmat-1)*(Adjmat) - (Gamma[i]-1)/(1+Gamma[i])*Adjmat
      perm = colMeans(Dmat)/(sqrt(colVars(Dmat)/I)) 
      pval[i] = (1+sum(perm>=obs))/(nperm + 1)
    }
  }else{
    for(i in 1:length(Gamma))
    {
      D2 = (Adjust) - (Gamma[i]-1)/(1+Gamma[i])*abs(Adjust)
      seH = SEPoP(D2, index)
      obs = mean(D2)/(seH)
      Adjmat = matrix(abs(Adjust), I, nperm)
      Zmat = matrix(runif(I*nperm) <Gamma[i]/(1+Gamma[i]), I, nperm)
      Dmat = (2*Zmat-1)*(Adjmat) - (Gamma[i]-1)/(1+Gamma[i])*Adjmat
      sevec = SEPoPMat(Dmat, index)
      perm = colMeans(Dmat)/sevec 
      pval[i] = (1+sum(perm>=obs))/(nperm + 1)
    }
  }
  Pmatrix = cbind(Gamma, pval)
  colnames(Pmatrix) = c("Gamma", "P-value")
  if(Changepoint == T)
  {
    proceed = ERSensitivityPoP(Ytrue, Dtrue, index, null, alphatrue, alternative, Gamma=1, nperm, Changepoint = F, SensitivityInterval = F)$pval <= alpha
    change = 1
    if(proceed)
    {
      change = uniroot(ERChangepointPoP, interval = c(1, 30), Y = Ytrue, D = Dtrue, index = index, null = null, alpha = alphatrue, alternative = alternative, nperm = nperm, extendInt = "upX")$root
    } 	
  }	
  if(SensitivityInterval == T)
  {
    lb = rep(-Inf, length(Gamma))
    ub = rep(Inf, length(Gamma))
    for(i in 1:length(Gamma))
    {		
      
      UB = mean(Ytrue)/mean(Dtrue) + 4*sd(Ytrue)/sqrt(I)
      LB = mean(Ytrue)/mean(Dtrue) - 4*sd(Ytrue)/sqrt(I)
      SUB = Inf
      SLB = -Inf
      if(alternative == "greater")
      {
        SLB = uniroot(ERSIPoP, interval = c(UB-4*sd(Ytrue)/sqrt(I), UB), extendInt = "upX", Gamma = Gamma[i], Y=Ytrue, D = Dtrue, index = index, alternative = "greater", alpha = alpha, nperm = nperm)$root
      }
      if(alternative == "less")
      {
        SUB = uniroot(ERSIPoP, interval = c(LB, LB + 4*sd(Ytrue)/sqrt(I)), extendInt = "downX", Gamma = Gamma[i], Y=Ytrue, D = Dtrue, index = index, alternative = "less", alpha = alpha, nperm = nperm)$root
      }
      if(alternative == "two.sided")
      {
        SLB = uniroot(ERSIPoP, interval = c(UB-4*sd(Ytrue)/sqrt(I), UB), extendInt = "upX", Gamma = Gamma[i], Y=Ytrue, D = Dtrue, index = index, alternative = "greater", alpha = alpha, nperm = nperm)$root
        SUB = uniroot(ERSIPoP, interval = c(LB, LB+4*sd(Ytrue)/sqrt(I)), extendInt = "downX", Gamma = Gamma[i], Y=Ytrue, D = Dtrue, index = index, alternative = "less", alpha = alpha, nperm = nperm)$root
      }
      lb[i] = SLB
      ub[i] = SUB
    }
    SImat = cbind(Gamma, lb, ub)
    colnames(SImat) = c("Gamma", "Lower Bound", "Upper Bound")
  }
  if(Changepoint == F & SensitivityInterval == F)
  {
    return(list(Gamma=Gamma, pval = pval, GammaPval = Pmatrix))
  }
  if(Changepoint == F & SensitivityInterval == T)
  {
    return(list(Gamma = Gamma, pval = pval, GammaPval = Pmatrix, SensitivityInterval = SImat))
  }
  if(Changepoint == T & SensitivityInterval == F)
  {
    return(list(Gamma = Gamma, pval = pval, GammaPval = Pmatrix, Changepoint = change))
  }
  if(Changepoint == T & SensitivityInterval == T)
  {
    return(list(Gamma = Gamma, pval = pval, GammaPval = Pmatrix, Changepoint = change, SensitivityInterval = SImat))
  }
  
}




####These are auxiliary functions used for root finding and for calculating columnwise variances in StudentizedSensitivity
ERChangepointPoP = function(Gamma, Y, D, index, null, alternative, alpha, nperm)
{
  alphachange = alpha
  if(alternative=="two.sided")
  {
    alphachange = alpha/2
  }
  ERSensitivityPoP(Y, D, index, null, alpha, alternative, Gamma, nperm, Changepoint = F, SensitivityInterval = F)$pval - alphachange 
}

ERSIPoP = function(null,  Gamma, Y, D, index, alternative, alpha, nperm)
{
  ERSensitivityPoP(Y, D, index, null, alpha, alternative, Gamma, nperm, Changepoint = F, SensitivityInterval = F)$pval - alpha  
  
}

colVars <- function(x) {
  N = nrow(x)
  (colSums(x^2) - colSums(x)^2/N) / (N-1)
}



indexPoP= function(X.t, X.c)
{
  require(nbpMatching)
  X.avg = (X.t + X.c)/2
  X.avg = as.matrix(X.avg)
  Sigmahat = cov(X.avg)
  Sinv = solve(Sigmahat)
  npair = nrow(X.avg)
  dist.pair = matrix(0, npair, npair)
  for(i in 1:npair)
  {
    dist.pair[i,] = mahalanobis(X.avg, X.avg[i,], Sinv, inverted = T)
    
  }
  dist.pair = dist.pair*1000
  dPoP = suppressWarnings(distancematrix(dist.pair))
  nbm = nonbimatch(dPoP)
  Pairs = as.matrix(nbm$halves[,c(2,4)])
  Pairsfin = Pairs
  if(npair%%2==1)
  {
    ghostind = which(Pairs == max(Pairs), arr.ind = T)
    notmatched = Pairs[ghostind[1], -ghostind[2]]
    Pairs2 = Pairs[-ghostind[1],]
    PoPavg = (X.avg[Pairs2[,1],,drop=F] + X.avg[Pairs2[,2],,drop=F])/2
    dist.trip = mahalanobis(PoPavg, X.avg[notmatched,,drop = F], Sinv, inverted = T)
    tripind = which.min(dist.trip)
    Pairsfin = Pairs2[-tripind,]
    triple = c(Pairs2[tripind,], notmatched)
  }
  index = rep(0, npair)
  for(i in 1:nrow(Pairsfin))
  {
    index[Pairsfin[i,]] = i
  }
  if(npair%%2==1)
  {
    index[triple] = (nrow(Pairsfin) + 1)
  }
  
  index
}

PoPMatch = function(X.t, X.c, nbm = NULL)
{
  X.avg = (X.t + X.c)/2
  X.avg = as.matrix(X.avg)
  Sigmahat = cov(X.avg)
  Sinv = solve(Sigmahat)
  npair = nrow(X.avg)
  dist.pair = matrix(0, npair, npair)
  for(i in 1:npair)
  {
    dist.pair[i,] = mahalanobis(X.avg, X.avg[i,], Sinv, inverted = T)
    
  }
  dist.pair = dist.pair*1000
  dPoP = suppressWarnings(distancematrix(dist.pair))
  nbm = nonbimatch(dPoP)
}


QPoP= function(X.t, X.c, nbm = NULL)
{
  X.avg = (X.t + X.c)/2
  X.avg = as.matrix(X.avg)
  Sigmahat = cov(X.avg)
  Sinv = solve(Sigmahat)
  npair = nrow(X.avg)
  if(is.null(nbm))
  {
    dist.pair = matrix(0, npair, npair)
    for(i in 1:npair)
    {
      dist.pair[i,] = mahalanobis(X.avg, X.avg[i,], Sinv, inverted = T)
      
    }
    dist.pair = dist.pair*1000
    dPoP = suppressWarnings(distancematrix(dist.pair))
    nbm = nonbimatch(dPoP)
  }
  Pairs = as.matrix(nbm$halves[,c(2,4)])
  Pairsfin = Pairs
  if(npair%%2==1)
  {
    ghostind = which(Pairs == max(Pairs), arr.ind = T)
    notmatched = Pairs[ghostind[1], -ghostind[2]]
    Pairs2 = Pairs[-ghostind[1],]
    PoPavg = (X.avg[Pairs2[,1],,drop=F] + X.avg[Pairs2[,2],,drop=F])/2
    dist.trip = mahalanobis(PoPavg, X.avg[notmatched,,drop = F], Sinv, inverted = T)
    tripind = which.min(dist.trip)
    Pairsfin = Pairs2[-tripind,]
    triple = c(Pairs2[tripind,], notmatched)
  }
  Q = matrix(0, npair, floor(npair/2))
  for(i in 1:nrow(Pairsfin))
  {
    Q[Pairsfin[i,],i] = 1
  }
  if(npair%%2==1)
  {
    Q[triple,ncol(Q)] = 1
  }
  
  Q
}




SEPoP = function(L, index)
{
  npair = length(L)
  ni = table(index)
  vv = tapply(L, index, FUN=var2)*ni
  se = (1/npair)*sqrt(sum(vv))
  se
}
SEPoPMat = function(Lmat, index)
{
  Lmat = as.matrix(Lmat)
  npair = nrow(Lmat)
  npop = floor(npair/2)
  ni = rep(2, floor(npair/2))
  if(npair%%2==1)
  {
    ni[length(ni)] = 3
  }
  VV = matrix(0, npop, ncol(Lmat))
  for(i in 1:(npop))
  {
    VV[i,] = (ni[i]^2/(ni[i]-1))*(colMeans(Lmat[index==i,,drop=F]^2) - colMeans(Lmat[index==i,,drop=F])^2)
  }
  sevec = (1/npair)*sqrt(colSums(VV))
  sevec
}

SEPoP2 = function(L, index)
{
  npair = length(L)
  ni = rep(2, floor(npair/2))
  if(npair%%2==1)
  {
    ni[length(ni)] = 3
  }
  mult = (ni/(ni-1))[index]
  mm = tapply(L, index, mean)
  res = as.vector(L) - mm[index]
  vv = sum(mult*res^2)
  se = (1/npair)*sqrt(sum(vv))
  se
}

var2 = function(x)
{
  (mean(x^2) - mean(x)^2)*length(x)/(length(x)-1)
}

SEPoP2= function(L, X.t, X.c, nbm = NULL)
{
  X.avg = (X.t + X.c)/2
  X.avg = as.matrix(X.avg)
  Sigmahat = cov(X.avg)
  Sinv = solve(Sigmahat)
  npair = nrow(X.avg)
  if(is.null(nbm))
  {
    dist.pair = matrix(0, npair, npair)
    for(i in 1:npair)
    {
      dist.pair[i,] = mahalanobis(X.avg, X.avg[i,], Sinv, inverted = T)
      
    }
    dist.pair = dist.pair*1000
    dPoP = suppressWarnings(distancematrix(dist.pair))
    nbm = nonbimatch(dPoP)
  }
  Pairs = as.matrix(nbm$halves[,c(2,4)])
  Pairsfin = Pairs
  if(npair%%2==1)
  {
    ghostind = which(Pairs == max(Pairs), arr.ind = T)
    notmatched = Pairs[ghostind[1], -ghostind[2]]
    Pairs2 = Pairs[-ghostind[1],]
    PoPavg = (X.avg[Pairs2[,1],,drop=F] + X.avg[Pairs2[,2],,drop=F])/2
    dist.trip = mahalanobis(PoPavg, X.avg[notmatched,,drop = F], Sinv, inverted = T)
    tripind = which.min(dist.trip)
    Pairsfin = Pairs2[-tripind,]
    triple = c(Pairs2[tripind,], notmatched)
  }
  vv = 0
  for(i in 1:nrow(Pairsfin))
  {
    vv = vv + 2*sum((L[Pairsfin[i,]] - mean(L[Pairsfin[i,]]))^2)
  }
  if(npair%%2==1)
  {
    vv = vv + (3/2)*sum((L[triple] - mean(L[triple]))^2)
  }
  
  se = 1/npair*sqrt(vv)
  se
}




