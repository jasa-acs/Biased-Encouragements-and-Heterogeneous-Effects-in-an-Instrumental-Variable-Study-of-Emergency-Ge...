
nsim = 10000
source("./rfunctions/ERSensitivityHet.R")
source("./rfunctions/ERSensitivityPoP.R")


simSettings = letters[1:16]
EstSize = matrix(0, 16, 3)
CILength = matrix(0, 16, 3)
#settings a-p are each of the rows in Table 4
for(ii in 1:length(simSettings))
{
  n= 0
  K = 0
  a = 0
    if(simSettings[ii] == "a")
    {
      K = 5
      n=100
      a = 1
    }
    if(simSettings[ii] == "b")
    {
      K = 10
      n=100
      a = 1
    }  
    if(simSettings[ii] == "c")
    {
      K = 5
      n=300
      a = 1
    }
    if(simSettings[ii] == "d")
    {
      K = 10
      n=300
      a = 1
    }  
    if(simSettings[ii] == "e")
    {
      K = 5
      n=1000
      a = 1
    }
    if(simSettings[ii] == "f")
    {
      K = 10
      n=1000
      a = 1
    }  
    if(simSettings[ii] == "g")
    {
      K = 5
      n=2000
      a = 1
    }
    if(simSettings[ii] == "h")
    {
      K = 10
      n=2000
      a = 1
    }
    
    if(simSettings[ii] == "i")
    {
      K = 5
      n=100
      a = 2
    }
    if(simSettings[ii] == "j")
    {
      K = 10
      n=100
      a = 2
    }  
    if(simSettings[ii] == "k")
    {
      K = 5
      n=300
      a = 2
    }
    if(simSettings[ii] == "l")
    {
      K = 10
      n=300
      a = 2
    }
    
    if(simSettings[ii] == "m")
    {
      K = 5
      n=1000
      a = 2
    }
    if(simSettings[ii] == "n")
    {
      K = 10
      n=1000
      a = 2
    }
    if(simSettings[ii] == "o")
    {
      K = 5
      n=2000
      a = 2
    }
    if(simSettings[ii] == "p")
    {
      K = 10
      n=2000
      a = 2
    }
    

B=n
I = 2*B
block.ind = rep(1:B, each = 2)
rejstandard = rep(0, nsim)
rejlin = rep(0, nsim)
rejPoP = rep(0, nsim)
SIstandard = matrix(0, nsim, 2)
SIlin = SIstandard
SIPoP = SIstandard

for(i in 1:nsim)
{
  X = matrix(runif(B*K), B, K)
  mu = 1*(10*sin(pi*X[block.ind,1]*X[block.ind,2]) + 20*(X[block.ind,3] - .5)^2 + 10*exp((X[block.ind,3])) + 5*(X[block.ind,5]-.5)^3)
  assign = rep(0, I)
  X2 = cbind(sin(pi*X[,1]*X[,2]), (X[,3] - .5)^2, exp((X[,4])), (X[,5]-.5)^3)
  assign = rep(0, I)
  epsilon = rnorm(I, 0, 1)
  piC = .58
  piA = (1-piC)/2
  piN = piA
  pS = c(piC + piN*piA, 1 - piC - 2*piA*piN, piA*piN)
  status = sample(c("A", "C", "N"), 2*B, replace = T, prob = c(piA, piC, piN))
  D = matrix(0, 2*B, 2)
  D[status=="A",] = 1
  D[status=="N",] = 0
  D[status=="C", 1] = 0
  D[status=="C", 2] = 1
  dc = D[,1]
  dt = D[,2]
  rttemp = a*(10*sin(pi*X[block.ind,1]*X[block.ind,2]) + 20*(X[block.ind,3] - .5)^2 + 10*exp((X[block.ind,4]))+ 5*(X[block.ind,5]-.5)^3+epsilon) 
  rctemp= 10*sin(pi*X[block.ind,1]*X[block.ind,2]) + 20*(X[block.ind,3] - .5)^2   +  10*exp((X[block.ind,4])) +  5*(X[block.ind,5]-.5)^3+epsilon 
  rt = rctemp*(1-dt) + rttemp*dt
  rc = rctemp*(1-dc) + rttemp*dc
  lambda = sum(rt-rc)/sum(dt-dc)
  
  YT = rep(0, B)
  YC = rep(0, B)
  DT = rep(0, B)
  DC = rep(0, B)
  zetaT = DT
  zetaC = DC
  for(j in 1:B)
  {
    ind = which(block.ind == j)
    YT[j] = (rt[ind[1]] - rc[ind[2]])
    YC[j] = (rt[ind[2]] - rc[ind[1]])
    DT[j] = (dt[ind[1]] - dc[ind[2]])
    DC[j] = (dt[ind[2]] - dc[ind[1]])
  }

  probs = 1/2
  Z1 = (runif(B) <= probs)
  Ydobs = YT*Z1 + YC*(1-Z1)
  Ddobs = DT*Z1 + DC*(1-Z1)
  sensstandard = ERSensitivityHet(Ydobs, Ddobs, X = NULL, null = lambda, alpha = 0.1, alternative = "two.sided", Gamma=1, nperm = 1000, Changepoint = F, SensitivityInterval = T)
  rejstandard[i] =  sensstandard$pval <= alph
  SIstandard[i,] = sensstandard$SensitivityInterval[-1]
  senslin = ERSensitivityHet(Ydobs, Ddobs, X = X, null = lambda, alpha = 0.1, alternative = "two.sided", Gamma=1, nperm = 1000, Changepoint = F, SensitivityInterval = T)
  rejlin[i] = senslin$pval <= alph
  SIlin[i,] = senslin$SensitivityInterval[-1]
  ind = indexPoP(X, X)
  sensPoP = ERSensitivityPoP(Ydobs, Ddobs, index = ind, null = lambda, alpha = 0.1, alternative = "two.sided", Gamma=1, nperm = 1000, Changepoint = F, SensitivityInterval = T)
  rejPoP[i] = sensPoP$pval <= alph
  SIPoP[i,] = sensPoP$SensitivityInterval[-1]
  if(i%%50==0)
  {
    print(c(ii, i))
  }
}

CILength[ii,] = c(mean(SIstandard[,2] - SIstandard[,1]), mean(SIlin[,2] - SIlin[,1]), mean(SIPoP[,2] - SIPoP[,1]))
EstSize[ii,] = colMeans(cbind(rejstandard, rejlin, rejPoP))
}


save(CILength, EstSize, file="./routput/Sim1.RData")




