
EffModTestPoP = function(Y, D, index = NULL, beta = 0.01,  nperm = 1000)
{
   CIbeta = ERSensitivityPoP(Y, D, index = index, null = 0, alternative = "two.sided", alpha = beta, Gamma = 1, nperm = nperm, SensitivityInterval = T, Changepoint = F)$SensitivityInterval[-1]
   laminit = sum(Y)/sum(D)
   
   pval = -optim(laminit, Frand,  Y = Y, D = D, index=index, nperm = nperm,method = "L-BFGS-B", lower = CIbeta[1], upper = CIbeta[2])$value + beta
   pval   
}
  

Frand = function(lambda, Y, D, index, nperm = 1000)
{
  Adjust = Y - lambda*D
  I = length(Adjust)
  D2 = (Adjust)
  SSEPoP = SEPoP(D2, index)
  SSE = (I-1)*var(D2)
  obs = (SSE - SSEPoP)/(SSEPoP)
  Adjmat = matrix(abs(Adjust), I, nperm)
  Zmat = matrix(runif(I*nperm) <1/2, I, nperm)
  Dmat = (2*Zmat-1)*(Adjmat)
  SSEPoPvec = SEPoPMat(Dmat, index)
  SSEvec = colVars(Dmat)*(I-1)
  perm = (SSEvec - SSEPoPvec)/(SSEPoPvec)
  pval = (1+sum(perm>=obs))/(nperm + 1)
  -pval
}

#EffModTestPoP(Ydobs, Ddobs, index = ind)

