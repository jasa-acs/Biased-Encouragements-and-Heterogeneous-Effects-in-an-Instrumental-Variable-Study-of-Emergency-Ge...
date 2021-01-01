########################
#PermTIV
########################
#Function to perform a Sensitivity analysis on the effect ratio in a paired observational study
#under the assumption of proportional doses

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

PermTIV = function(Y, D, null = 0, alpha = 0.05, alternative = "greater", Gamma=1, nperm = 50000, Changepoint = T, SensitivityInterval = T)
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
  for(i in 1:length(Gamma))
  {
    D2 = (Adjust) - (Gamma[i]-1)/(1+Gamma[i])*abs(Adjust)
    obs = mean(D2)
    Adjmat = matrix(abs(Adjust), I, nperm)
    Zmat = matrix(runif(I*nperm) <Gamma[i]/(1+Gamma[i]), I, nperm)
    Dmat = (2*Zmat-1)*(Adjmat) - (Gamma[i]-1)/(1+Gamma[i])*Adjmat
    perm = colMeans(Dmat)
    pval[i] = (1+sum(perm>=obs))/(nperm + 1)
  }
  Pmatrix = cbind(Gamma, pval)
  colnames(Pmatrix) = c("Gamma", "P-value")
  if(Changepoint == T)
  {
    proceed = PermTIV(Ytrue, Dtrue, null, alphatrue, alternative, Gamma=1, nperm, Changepoint = F, SensitivityInterval = F)$pval <= alpha
    change = 1
    if(proceed)
    {
      change = uniroot(PermTChangepoint, interval = c(1, 30), Y = Ytrue, D = Dtrue, null = null, alpha = alphatrue, alternative = alternative, nperm = nperm, extendInt = "upX")$root
    } 	
  }	
  if(SensitivityInterval == T)
  {
    lb = rep(-Inf, length(Gamma))
    ub = rep(Inf, length(Gamma))
    for(i in 1:length(Gamma))
    {		
      #UB = uniroot(ERBoundFinder, Ytrue, Dtrue, Gamma[i], interval = c(mean(Ytrue)/mean(Dtrue), mean(Ytrue)/mean(Dtrue)+4*sd(Ytrue)/sqrt(I)), extendInt = "yes")$root
      #LB = -uniroot(ERBoundFinder, Ytrue, -Dtrue, Gamma[i], interval = c(mean(Ytrue)/mean(Dtrue)-4*sd(Ytrue)/sqrt(I), -mean(Ytrue)/mean(Dtrue)), extendInt = "yes")$root
      UB = mean(Ytrue)/mean(Dtrue) + 4*sd(Ytrue)/sqrt(I)
      LB = mean(Ytrue)/mean(Dtrue) - 4*sd(Ytrue)/sqrt(I)
      SUB = Inf
      SLB = -Inf
      if(alternative == "greater")
      {
        SLB = uniroot(PermTSI, interval = c(UB-4*sd(Ytrue)/sqrt(I), UB), extendInt = "upX", Gamma = Gamma[i], Y=Ytrue, D = Dtrue, alternative = "greater", alpha = alpha, nperm = nperm)$root
      }
      if(alternative == "less")
      {
        SUB = uniroot(PermTSI, interval = c(LB, LB + 4*sd(Ytrue)/sqrt(I)), extendInt = "downX", Gamma = Gamma[i], Y=Ytrue, D = Dtrue, alternative = "less", alpha = alpha, nperm = nperm)$root
      }
      if(alternative == "two.sided")
      {
        SLB = uniroot(PermTSI, interval = c(UB-4*sd(Ytrue)/sqrt(I), UB), extendInt = "upX", Gamma = Gamma[i], Y=Ytrue, D = Dtrue, alternative = "greater", alpha = alpha, nperm = nperm)$root
        SUB = uniroot(PermTSI, interval = c(LB, LB+4*sd(Ytrue)/sqrt(I)), extendInt = "downX", Gamma = Gamma[i], Y=Ytrue, D = Dtrue, alternative = "less", alpha = alpha, nperm = nperm)$root
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
PermTChangepoint = function(Gamma, Y, D, null, alternative, alpha, nperm)
{
  alphachange = alpha
  if(alternative=="two.sided")
  {
    alphachange = alpha/2
  }
  PermTIV(Y, D, null, alpha, alternative, Gamma, nperm, Changepoint = F, SensitivityInterval = F)$pval - alphachange 
}

PermTSI = function(null,  Gamma, Y, D, alternative, alpha, nperm)
{
  PermTIV(Y, D, null, alpha, alternative, Gamma, nperm, Changepoint = F, SensitivityInterval = F)$pval - alpha  
  
}

colVars <- function(x) {
  N = nrow(x)
  (colSums(x^2) - colSums(x)^2/N) / (N-1)
}

