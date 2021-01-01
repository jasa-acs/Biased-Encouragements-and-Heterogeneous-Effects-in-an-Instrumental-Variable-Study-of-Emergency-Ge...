
 tsls_sum <- function(obj, null=0){
  model.summary = summary(obj)
  coefficientTable = model.summary$coefficients
  dfTest = model.summary$df[2]
  estValue = as.numeric(coefficientTable[2,1])
  stdError = as.numeric(coefficientTable[2,2])
  tValue = abs( (estValue - null) / stdError)
  pvalue = 2*(1 - pt(tValue,dfTest))
  confInt = as.numeric(confint(obj)[2,])
  out <- list(pe = estValue, ci = confInt, pval = pvalue)
  }   


gen_eff_ratio <- function(R, D, Z, matchedNumber, alphaLevel = 0.05, null=0){  
  #if(missing(matchedNumber)) matchedNumber = matching(Z,X,"full")
  	
  # Deal with non-matched individuals #
  R.matchedIndiv = R[!is.na(matchedNumber)]
  D.matchedIndiv = D[!is.na(matchedNumber)]
  #Z.matchedIndiv = Z[!is.na(matchedNumber)]
  matchedNumber.matchedIndiv = matchedNumber[!is.na(matchedNumber)]

  # Sort individuals in ascending matching order
  matchedNumber.sortedIndex = sort(matchedNumber.matchedIndiv,index.return=TRUE)$ix
  R.sorted = R[matchedNumber.sortedIndex]
  D.sorted = D[matchedNumber.sortedIndex]
  Z.sorted = Z[matchedNumber.sortedIndex]
  matchedNumber.sorted = matchedNumber.matchedIndiv[matchedNumber.sortedIndex]
  
  # Calculate the size of each matched set and the corresponding weights
  ni = tabulate(matchedNumber.sorted); 
  ni = ni[!(ni == 0)] #this cleans up a flaw with matching algorithm 
                    #where some matched numbers (e.g. 1:I) are not used between
					#1 to I
  I = length(ni)
  wi = ni^2 / (ni - 1)  #Actual formula: ni^2 /(mi * (ni - mi))

  # Calculate Vi, Gi, and Hi
  # Formula
  # Vi = Gi - lambda * Hi
  # Gi = wi*(Zij - Zimean)(Rij - Rimean) = wi*(ZijRij - ni*Zmean.*Rmean.)
  # Hi = wi*(Zij - Zimean)(Dij - Dimean) = wi*(ZijDij - ni*Zmean.*Dmean.)

  ni.cumsum = cumsum(ni)
  RZ.cumsum = cumsum(R.sorted*Z.sorted); 
  DZ.cumsum = cumsum(D.sorted*Z.sorted)
  R.cumsum = cumsum(R.sorted); D.cumsum = cumsum(D.sorted); Z.cumsum = cumsum(Z.sorted)
  R.means = (R.cumsum[ni.cumsum] - c(0,R.cumsum[ni.cumsum[-length(ni.cumsum)]]))/ni
  Z.means = (Z.cumsum[ni.cumsum] - c(0,Z.cumsum[ni.cumsum[-length(ni.cumsum)]]))/ni
  D.means = (D.cumsum[ni.cumsum] - c(0,D.cumsum[ni.cumsum[-length(ni.cumsum)]]))/ni
  
  RZ.group = RZ.cumsum[ni.cumsum] - c(0,RZ.cumsum[ni.cumsum[-length(ni.cumsum)]])
  DZ.group = DZ.cumsum[ni.cumsum] - c(0,DZ.cumsum[ni.cumsum[-length(ni.cumsum)]])
  
  if(length(R.means) != length(Z.means)) stop("boo!")
  if(length(RZ.group) != length(Z.means)) stop("boo!!")
  if(length(DZ.group) != length(Z.means)) stop("boo!!!")

  Gi = wi * RZ.group - wi * ni * R.means * Z.means
  Hi = wi * DZ.group - wi * ni * D.means * Z.means
  
  # Compute the point estimate 
  pointEst = sum(Gi) / sum(Hi)
  
  # Compute the p-value under the null
  ViNull = Gi - null * Hi
  testStatNull = mean(ViNull) / sqrt(1/(I *(I-1)) * sum( (ViNull - mean(ViNull))^2) )
  pvalue = 2* (1 - pnorm(abs(testStatNull)))
  
# Compute the quadratic terms for CI
  # A*lambda^2 + B*lambda + C = 0 (technically, it's A*lambda^2 + B*lambda + C <= 0)
  q = qnorm(1 - alphaLevel/2)
  A = sum(Hi)^2/I^2 - q^2 /(I * (I-1)) * sum( (Hi - mean(Hi))^2)
  B = -2 * (sum(Hi) * sum(Gi) / I^2 - q^2 / (I * (I-1)) * sum( (Hi - mean(Hi)) * (Gi - mean(Gi))))
  C = sum(Gi)^2/I^2 - q^2/(I * (I-1)) * sum( (Gi - mean(Gi))^2)

  detQuad = round(B^2 - 4*A*C,9) #6 is set for numerical accuracy
  if( detQuad <= 0) {
    if(A < 0) {
	  cis = matrix(c(-Inf,Inf),1,2)
	} else {
	  cis = matrix(c(NA,NA),1,2)
	} 
  }
  if(detQuad > 0) {
    if(A < 0) {
	  up.ci = (-B - sqrt(detQuad))/ (2*A) 
      low.ci = (-B + sqrt(detQuad)) / (2*A) 
	  cis = matrix(c(-Inf,low.ci,up.ci,Inf),2,2,byrow=TRUE)
	} else {
	  low.ci = (-B - sqrt(detQuad))/ (2*A) 
      up.ci = (-B + sqrt(detQuad)) / (2*A) 
	  cis = matrix(c(low.ci,up.ci),1,2)
	}
  }

  
  out <- list(pe=pointEst, ci.l = cis[1], ci.u = cis[2], pval= pvalue)
 }
 
 
 sensInt <- function(Gamma,Z,R,pairmatchvec) {
  matchedSets = unique(pairmatchvec[!is.na(pairmatchvec)]); 
  I <- length(matchedSets)
  cs.plus = rep(0,length(I))
  ps.plus = rep(0,length(I))
  ps.minus = rep(0,length(I))
  testStat = rep(0,length(I))
  ds= rep(0,length(I))
  for(i in 1:length(matchedSets)) { 
    index = matchedSets[i]
    matchedSetID = which(pairmatchvec == index)
    ni = length(matchedSetID); 
    mi = sum(Z[matchedSetID])	  
    cs.plus[i] = sum(R[matchedSetID])
    ps.plus[i] = Gamma * cs.plus[i]/(Gamma*cs.plus[i] + ni - cs.plus[i])
    ps.minus[i] = cs.plus[i]/(cs.plus[i] + (ni - cs.plus[i])*Gamma)
    ds[i] = ni/(mi*(ni - mi))
    testStat[i] = ds[i] * sum(Z[matchedSetID] * R[matchedSetID])
  }
  T = sum(testStat)
  lowerBound = (T - sum(ds * ps.plus)) /(sqrt(sum(ds^2 * ps.plus * (1 - ps.plus))))
  upperBound = (T - sum(ds * ps.minus)) /(sqrt(sum(ds^2 * ps.minus * (1 - ps.minus))))
  return(c(lowerBound,upperBound))
}

