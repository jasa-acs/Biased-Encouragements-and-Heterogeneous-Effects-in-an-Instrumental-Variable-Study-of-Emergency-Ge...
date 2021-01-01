
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
  vv = tapply(L, index, FUN=var)*ni
  se = (1/npair)*sqrt(sum(vv))
  se
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



