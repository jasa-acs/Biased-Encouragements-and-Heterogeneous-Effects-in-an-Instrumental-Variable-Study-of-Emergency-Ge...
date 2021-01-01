library(rmutil)

design.sens.IV = function(m, sigma2, piC, piA = (1-piC)/2, piN = (1-piC)/2)
{
  if(round(piA+piC+piN,2) != 1)
  {
    stop("pi must be probability dist")
  }
  s = sqrt(sigma2)
  Ezeta = piC*m
  Eabs = (piC+2*piA*piN)*(s*sqrt(2/pi)*exp(-m^2/(2*s^2)) + m*(1 - 2*pnorm(-m/s))) + (1- piC - 2*piA*piN)*(s*sqrt(2/pi)*exp(-0^2/(2*s^2)) + 0*(1 - 2*pnorm(-0/s)))
  
  GP = (Eabs+Ezeta)/(Eabs - Ezeta)
  GP
}
integrand = function(y,m, b)
{
  abs(y)*dlaplace(y=y,m=m, s=b)
}
Eabslaplace = function(m, sigma2)
{
  integrate(integrand, lower = -Inf, upper = Inf, m = m, b = sqrt(sigma2/2))$value
}

design.sens.IV.Lap = function(m, sigma2, piC, piA = (1-piC)/2, piN = (1-piC)/2)
{
  if(round(piA+piC+piN,2) != 1)
  {
    stop("pi must be probability dist")
  }
  s = sqrt(sigma2)
  Ezeta = piC*m
  Eabs = (piC+2*piA*piN)*(Eabslaplace(m, sigma2)) + (1- piC - 2*piA*piN)*Eabslaplace(0, sigma2)
  
  GP = (Eabs+Ezeta)/(Eabs - Ezeta)
  GP
}


design.sens.P = function(m, sigma2)
{
  s = sqrt(sigma2)
  GP = ((s*sqrt(2/pi)*exp(-m^2/(2*s^2)) + m*(1 - 2*pnorm(-m/s))) + m)/((s*sqrt(2/pi)*exp(-m^2/(2*s^2)) + m*(1 - 2*pnorm(-m/s))) - m)
  GP
}

design.sens.P(1/2, 1/2)
design.sens.IV(1/2, 1/2, 1)
design.sens.IV(1/2, 1/2, .5)

# Table 6 Row 1
design.sens.IV(6.8, 25.3^2, 1)
design.sens.IV(6.8, 25.3^2, .75)
design.sens.IV(6.8, 25.3^2, .58)
design.sens.IV(6.8, 25.3^2, .5)
design.sens.IV(6.8, 25.3^2, .25)
design.sens.IV(6.8, 25.3^2, .1)

design.sens.IV.Lap(6.8, 25.3^2, 1)
design.sens.IV.Lap(6.8, 25.3^2, .75)
design.sens.IV.Lap(6.8, 25.3^2, .58)
design.sens.IV.Lap(6.8, 25.3^2, .5)
design.sens.IV.Lap(6.8, 25.3^2, .25)
design.sens.IV.Lap(6.8, 25.3^2, .1)

# Table 6 Row 2
design.sens.IV(4.1, 8.9^2, 1)
design.sens.IV(4.1, 8.9^2, .75)
design.sens.IV(4.1, 8.9^2, .58)
design.sens.IV(4.1, 8.9^2, .5)
design.sens.IV(4.1, 8.9^2, .25)
design.sens.IV(4.1, 8.9^2, .1)

design.sens.IV.Lap(4.1, 8.9^2, 1)
design.sens.IV.Lap(4.1, 8.9^2, .75)
design.sens.IV.Lap(4.1, 8.9^2, .58)
design.sens.IV.Lap(4.1, 8.9^2, .5)
design.sens.IV.Lap(4.1, 8.9^2, .25)
design.sens.IV.Lap(4.1, 8.9^2, .1)



#simulation 2 - potential for improved power
lambdaS = 6.8
sigma2S = 25.3^2


Ifactor = 8
lambdaNS = 4.1
sigma2NS = 8.9^2
nsim = 5000
gam = seq(1, 2.5, .05)
power.permt = rep(0, nsim)
power.ate = rep(0, nsim)
power.student = rep(0, nsim)
student = rep(0, length(gam))
var.perm = rep(0, nsim)
var.ate = rep(0, nsim)
SumD = rep(0, nsim)

piC = .58

Ivec = c(200, 400,1000, 2000, 20000)
#powerS = matrix(0, length(gam), length(Ivec))
#powerNS = powerS

for(j in 1:1)
{
I = Ivec[j]
lambda = lambdaS
sigma2 = sigma2S


for(a in 1:length(gam))
{
  g = gam[a]
  p = g/(1+g)

  piA = (1-piC)/2
  piN = piA
  pS = c(piC + piN*piA, 1 - piC - 2*piA*piN, piA*piN)



  for(i in 1:nsim)
  {
    eps = rnorm(I, 0, sqrt(sigma2))


    S = sample(c(1,0,-1), size = I, replace = T, prob = pS)

    Ydiff = eps + S*lambda
    ss = ERSensitivity(Y=Ydiff, D = S, null = 0, Gamma=g, Changepoint=F, SensitivityInterval=F, nperm = 2000)
    power.student[i] = (ss$pval <=.05)
    if(i%%50 == 0)
    {
      print(c(j,a,i,1))
    }
  }
  student[a] = mean(power.student)
}
powerS[,j] = student

I=floor(I/Ifactor)
lambda = lambdaNS
sigma2 = sigma2NS

for(a in 1:length(gam))
{
  g = gam[a]
  p = g/(1+g)

  piA = (1-piC)/2
  piN = piA
  pS = c(piC + piN*piA, 1 - piC - 2*piA*piN, piA*piN)



  for(i in 1:nsim)
  {
    eps = rnorm(I, 0, sqrt(sigma2))


    S = sample(c(1,0,-1), size = I, replace = T, prob = pS)

    Ydiff = eps + S*lambda
    ss = ERSensitivity(Y=Ydiff, D = S, null = 0, Gamma=g, Changepoint=F, SensitivityInterval=F, nperm = 2000)
    power.student[i] = (ss$pval <=.05)
    if(i%%50 == 0)
    {
      print(c(j,a,i,2))
    }
  }
  student[a] = mean(power.student)
}
powerNS[,j] = student

}


dS = design.sens.IV(6.8, 25.3^2, .58)
dNS = design.sens.IV(4.1, 8.9^2, .58)

pdf("./routput/figures/powerSepsis.pdf")
par(mfrow = c(2,2), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1.5,.5), oma= c(0,0,2,0))
gamplot = gam
par(mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1.5,.5))
plot(gamplot, powerS[,1], col = "orange", ylim = c(0,1), type = "l", lwd = 3, ylab = "Power", xlab = expression(paste(Gamma)), main =  "I = 200")
lines(gamplot, powerNS[,1], col = "black", lwd = 3, lty = 2)

plot(gamplot, powerS[,3], col = "orange", ylim = c(0,1), type = "l", lwd = 3, ylab = "Power", xlab = expression(paste(Gamma)), main =  "I = 1000")
lines(gamplot, powerNS[,3], col = "black", lwd = 3, lty = 2)
legend("topright", c(expression(paste("Septic (", lambda, "= 6.8, ", sigma, "=25.3, ","n=I)", sep= "")), expression(paste("Non-septic (", lambda, "= 4.1, ", sigma, "=8.9, ","n=I/8)", sep = ""))), lty = c(1,2,6,3), lwd =2, col = c("orange", "black", "blue"),  bty = "n")

plot(gamplot, powerS[,4], col = "orange", ylim = c(0,1), type = "l", lwd = 3, ylab = "Power", xlab = expression(paste(Gamma)), main =  "I = 2000")
lines(gamplot, powerNS[,4], col = "black", lwd = 3, lty = 2)

plot(gamplot, powerS[,5], col = "orange", ylim = c(0,1), type = "l", lwd = 3, ylab = "Power", xlab = expression(paste(Gamma)), main =  "I = 10000")
lines(gamplot, powerNS[,5], col = "black", lwd = 3, lty = 2)

dev.off()


