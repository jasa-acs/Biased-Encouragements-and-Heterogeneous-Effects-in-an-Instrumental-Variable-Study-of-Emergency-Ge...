###########
#sensEffectRatio
###########
sensEffectRatio = function(index, treatment, outcome, dose, null=0, DE = "both", MO = T, ER = T, alternative = "two.sided", alpha = 0.05, Gamma.vec = 1, calculate.pval = T, continuous.relax = F)
{
  PVAL = calculate.pval
  require(gurobi)
  require(Matrix)
  
  
  
  ns = table(index)
  ms = table(index[treatment*1==1]) 
  N.total = sum(ns)
  
  
  nostratum = length(unique(index))
  treatment = 1*treatment
  treatment = (1*treatment == 1)
  gur = suppressWarnings(require(gurobi))
  
  
  vec.111 = tapply((treatment)*outcome*(dose), index, sum)
  vec.110 = tapply((treatment)*outcome*(!dose), index, sum)
  vec.010 = tapply((!treatment)*outcome*(!dose), index, sum)
  vec.101 = tapply((treatment)*(!outcome)*dose, index, sum)
  vec.100 = tapply((treatment)*(!outcome)*(!dose), index, sum)
  vec.001 = tapply((!treatment)*(!outcome)*dose, index, sum)
  vec.011 = tapply((!treatment)*(outcome)*dose, index, sum)
  vec.000 = tapply((!treatment)*(!outcome)*(!dose), index, sum)
  V = cbind(vec.111, vec.110, vec.101, vec.100, vec.011, vec.010,vec.001, vec.000)
  id = apply(V, 1, paste, collapse = "-")
  
  num.id =  sort.new(xtfrm(id))
  nosymm = length(unique(num.id))
  cc = table(num.id)
  bb = tapply(ns, num.id, mean)
  N.sr = sum(cc-1)
  m.011 = vec.011+1
  m.010 = vec.010+1
  m.001 = vec.001+1
  m.000 = vec.000+1
  m.101 = vec.101 + 1
  m.100 = vec.100 + 1
  m.111 = vec.111 + 1
  m.110 = vec.110 + 1
  
  
  mult.000 = tapply(m.000, num.id, mean)
  mult.001 = tapply(m.001, num.id, mean)
  mult.010 = tapply(m.010, num.id, mean)
  mult.011 = tapply(m.011, num.id, mean)
  mult.110 = tapply(m.110, num.id, mean)
  mult.111 = tapply(m.111, num.id, mean)
  mult.101 = tapply(m.101, num.id, mean)
  mult.100 = tapply(m.100, num.id, mean)
  
  wATE.per.strat = (ns)*(tapply((treatment)*outcome, index, sum)/ms - tapply((1-treatment)*outcome, index, sum)/(ns-ms))
  ATE.est = sum(wATE.per.strat)
  wATE2.per.strat = (ns)*(tapply((treatment)*dose, index, sum)/ms - tapply((1-treatment)*dose, index, sum)/(ns-ms))
  ATE2.est = sum(wATE2.per.strat)
  
  lambdahat = ATE.est/ATE2.est
  #Gamma.vec = round(Gamma.vec, 2)
  null = round(null, 2)
  
  
  wTE.per.strat = (tapply((treatment)*(outcome - null*dose), index, sum)/ms - tapply((1-treatment)*(outcome - null*dose), index, sum)/(ns-ms))
  
  
  
  TE.est = sum(wTE.per.strat)
  max.e = (TE.est > 0)
  
  pvalvec = rep(0, length(Gamma.vec))
  Rejectvec = rep(0, length(Gamma.vec))
  kappavec = rep(0, length(Gamma.vec))
  SE = 0
  if(ER == F)
  {
    mult.po = (mult.111)*mult.110*(mult.101)*mult.100*mult.001*(mult.000)*(mult.010)*mult.011
    mult.do = (mult.111)*mult.110*(mult.101)*mult.100*mult.001*(mult.000)*(mult.010)*mult.011
    if(DE == "nonpositive")
    {
      mult.po = mult.100*mult.101*mult.010*mult.011
    }
    if(DE=="nonnegative")
    {
      mult.po = mult.111*mult.110*mult.000*mult.001
    }
    if(MO == T)
    {
      mult.do = mult.111*mult.101*mult.000*mult.010
    }
    mult.all = mult.po*mult.do
    ns.type = tapply(ns, num.id, mean)
    N.vars = sum(mult.all*(ns.type-1))
    index.symm = rep(1:nosymm, mult.all*(ns.type-1))
    n.types = (mult.all)
    n.po = mult.po
    n.do= mult.do
    n.per = cc
    Diff = rep(0, N.vars)
    Diff2 = Diff
    #V.list = vector("list", N.vars)
    row.ind = rep(0, 4*N.vars+1)
    col.ind = row.ind
    values = row.ind
    b = rep(0, nosymm + 2)
    for(kk in 1:nosymm)
    {
      row.ind[which(index.symm==kk)] = rep(kk, (ns.type[kk]-1)*n.types[kk])
      col.ind[which(index.symm==kk)] = which(index.symm==kk)  
      values[which(index.symm==kk)] = rep(1, (ns.type[kk]-1)*(n.types[kk]))
      b[kk] = n.per[kk]
    }
    row.ind[(N.vars+1):(2*N.vars)] = rep(nosymm + 1, N.vars)
    col.ind[(N.vars+1):(2*N.vars)] = 1:N.vars
    row.ind[(2*N.vars+1):(3*N.vars)] = rep(nosymm + 2, N.vars)
    col.ind[(2*N.vars+1):(3*N.vars)] = 1:N.vars
    row.ind[(3*N.vars+1):(4*N.vars+1)] = rep(nosymm + 3, N.vars+1)
    col.ind[(3*N.vars+1):(4*N.vars+1)] = 1:(N.vars+1)
    PM = rep(0, N.vars)
    PV = rep(0, N.vars)
    
    
    
    #Gamma.vec = seq(1.2, 1.22, by = .001)
    
    zscore = rep(0, length(Gamma.vec))
    
    for(ee in 1:length(Gamma.vec))
    {
      Gamma.sens = Gamma.vec[ee]
      
      for(kk in 1:nosymm)
      {
        i = which(num.id==kk)[1]
        symmgroup = which(index.symm == kk)
        ind = which(index==i)
        treatstrat = treatment[ind]
        outstrat = outcome[ind]
        dosestrat = dose[ind]
        dosesymm = c(rep(0, sum(treatstrat==F&dosestrat==F)), rep(1, sum(treatstrat==F& dosestrat==T)), rep(0, sum(treatstrat==T&dosestrat==F)), rep(1, sum(treatstrat==T& dosestrat==T)))
        
        outsymm = c(sort(outstrat[treatstrat==F& dosestrat==F]), sort(outstrat[treatstrat==F& dosestrat==T]), sort(outstrat[treatstrat==T&dosestrat == F]), sort(outstrat[treatstrat==T&dosestrat == T]))
        PO = matrix(0, n.po[kk], ns[i])
        DO = matrix(0, n.do[kk],ns[i])
        
        mt.011 = mult.011[kk]
        mt.001 = mult.001[kk]
        mt.000 = mult.000[kk]
        mt.010 = mult.010[kk]
        mt.110 = mult.110[kk]
        mt.100 = mult.100[kk]
        mt.111 = mult.111[kk]
        mt.101 = mult.101[kk]
        p.011 = mt.011
        p.001 = mt.001
        p.000 = mt.000
        p.010 = mt.010
        p.110 = mt.110
        p.100 = mt.100
        p.111 = mt.111
        p.101 = mt.101
        
        T.000 = matrix(1, mt.000-1, mt.000-1)
        T.000[lower.tri(T.000)] = 0
        T.000 = rbind(T.000, c(rep(0, mt.000-1)) )
        T.001 = matrix(1, mt.001-1, mt.001-1)
        T.001[lower.tri(T.001)] = 0
        T.001 = rbind(T.001, c(rep(0, mt.001-1)) )
        T.010 = matrix(1, mt.010-1, mt.010-1)
        T.010[lower.tri(T.010)] = 0
        T.010 = rbind(T.010, c(rep(0, mt.010-1)) )
        T.011 = matrix(1, mt.011-1, mt.011-1)
        T.011[lower.tri(T.011)] = 0
        T.011 = rbind(T.011, c(rep(0, mt.011-1)))
        T.100 = matrix(1, mt.100-1, mt.100-1)
        T.100[lower.tri(T.100)] = 0
        T.100 = rbind(T.100, c(rep(0, mt.100-1)) )
        T.101 = matrix(1, mt.101-1, mt.101-1)
        T.101[lower.tri(T.101)] = 0
        T.101 = rbind(T.101, c(rep(0, mt.101-1)) )
        T.110 = matrix(1, mt.110-1, mt.110-1)
        T.110[lower.tri(T.110)] = 0
        T.110 = rbind(T.110, c(rep(0, mt.110-1)) )
        T.111 = matrix(1, mt.111-1, mt.111-1)
        T.111[lower.tri(T.111)] = 0
        T.111 = rbind(T.111, c(rep(0, mt.111-1)))
        
        if(DE == "nonnegative")
        {
          T.100 = matrix(0, 1, mt.100-1)
          p.100 = 1
          T.101 = matrix(0, 1, mt.101-1)
          p.101 = 1
          T.011 = matrix(1, 1, mt.011-1)
          p.011 = 1
          T.010 = matrix(1, 1, mt.010-1)
          p.010 = 1
        }
        if(DE == "nonpositive")
        {
          T.111 = matrix(1, 1, mt.111-1)
          p.111 = 1
          T.110 = matrix(1, 1, mt.110-1)
          p.110 = 1
          T.001 = matrix(0, 1, mt.001-1)
          p.001 = 1
          T.000 = matrix(0, 1, mt.000-1)
          p.000 = 1
        }
        count = 1
        for(ll in 1:p.000)
        {
          for(mm in 1:p.010)
          {
            for(uu in 1:p.001)
            {
              for(vv in 1:p.011)
              {
                for(ww in 1:p.100)
                {
                  for(xx in 1:p.110)
                  {
                    for(yy in 1:p.101)
                    {
                      for(zz in 1:p.111)	
                      {
                        tempvec = c(T.000[ll,], T.010[mm,], T.001[uu,], T.011[vv,], T.100[ww,], T.110[xx,], T.101[yy,], T.111[zz,])
                        tempvec = tempvec[!is.na(tempvec)]
                        PO[count,] = tempvec
                        count = count+1
                      }
                    }
                  }
                }
                
              }
              
            }
          }
          
        }
        count =1 
        p.011 = mt.011
        p.001 = mt.001
        p.000 = mt.000
        p.010 = mt.010
        p.110 = mt.110
        p.100 = mt.100
        p.111 = mt.111
        p.101 = mt.101
        
        T.000 = matrix(1, mt.000-1, mt.000-1)
        T.000[lower.tri(T.000)] = 0
        T.000 = rbind(T.000, c(rep(0, mt.000-1)) )
        T.001 = matrix(1, mt.001-1, mt.001-1)
        T.001[lower.tri(T.001)] = 0
        T.001 = rbind(T.001, c(rep(0, mt.001-1)) )
        T.010 = matrix(1, mt.010-1, mt.010-1)
        T.010[lower.tri(T.010)] = 0
        T.010 = rbind(T.010, c(rep(0, mt.010-1)) )
        T.011 = matrix(1, mt.011-1, mt.011-1)
        T.011[lower.tri(T.011)] = 0
        T.011 = rbind(T.011, c(rep(0, mt.011-1)))
        T.100 = matrix(1, mt.100-1, mt.100-1)
        T.100[lower.tri(T.100)] = 0
        T.100 = rbind(T.100, c(rep(0, mt.100-1)) )
        T.101 = matrix(1, mt.101-1, mt.101-1)
        T.101[lower.tri(T.101)] = 0
        T.101 = rbind(T.101, c(rep(0, mt.101-1)) )
        T.110 = matrix(1, mt.110-1, mt.110-1)
        T.110[lower.tri(T.110)] = 0
        T.110 = rbind(T.110, c(rep(0, mt.110-1)) )
        T.111 = matrix(1, mt.111-1, mt.111-1)
        T.111[lower.tri(T.111)] = 0
        T.111 = rbind(T.111, c(rep(0, mt.111-1)))
        
        if(MO == T)
        {
          T.100 = matrix(0, 1, mt.100-1)
          p.100 = 1
          T.110 = matrix(0, 1, mt.110-1)
          p.110 = 1
          T.011 = matrix(1, 1, mt.011-1)
          p.011 = 1
          T.001 = matrix(1, 1, mt.001-1)
          p.001 = 1
        }
        for(ll in 1:p.000)
        {
          for(mm in 1:p.010)
          {
            for(uu in 1:p.001)
            {
              for(vv in 1:p.011)
              {
                for(ww in 1:p.100)
                {
                  for(xx in 1:p.110)
                  {
                    for(yy in 1:p.101)
                    {
                      for(zz in 1:p.111)	
                      {
                        tempvec = c(T.000[ll,], T.010[mm,], T.001[uu,], T.011[vv,], T.100[ww,], T.110[xx,], T.101[yy,], T.111[zz,])
                        tempvec = tempvec[!is.na(tempvec)]
                        DO[count,] = tempvec
                        count = count+1
                      }
                    }
                  }
                }
                
              }
              
            }
          }
          
        }
        
        
        
        
        
        count = 1
        
        for(jj in 1:n.po[kk])
        {
          for(ll in 1:n.do[kk])
          {
            ind.jj = (((count-1)*(ns[i]-1))+1):(count*(ns[i]-1))
            po.symm = PO[jj,]
            do.symm = DO[ll,]
            treatsymm = c(rep(F, ns[i]-ms[i]), rep(T, ms[i]))
            outcontrol = outsymm*(1-treatsymm) + po.symm*(treatsymm)
            outtreat = outsymm*(treatsymm) + po.symm*(1-treatsymm) 
            dosecontrol = dosesymm*(1-treatsymm) + do.symm*(treatsymm)
            dosetreat = dosesymm*(treatsymm) + do.symm*(1-treatsymm) 
            
            sum.cont = sum(outcontrol - null*dosecontrol)/(ns[i]-1)
            Q = (outtreat + outcontrol/(ns[i]-1) - null*dosetreat - null*dosecontrol/(ns[i]-1)- sum.cont)*ns[i]
            if(sum(treatstrat)>1)
            {
              sum.cont = sum(outtreat - null*dosetreat)/(ns[i]-1)
              Q = -(outtreat/(ns[i]-1) + outcontrol - null*dosetreat/(ns[i]-1) - null*dosecontrol- sum.cont)*ns[i]
              
            }
            
            
            
            qi = Q*max.e - Q*(!max.e)
            ord = order(qi)
            qi.sort = sort(qi)
            
            eta = diff(outcontrol+outtreat - null*(dosecontrol+dosetreat))/2
            taubar = mean(outtreat-outcontrol - null*(dosetreat-dosecontrol))
            
            mu = rep(0, length(ind)-1)
            sigma2 = rep(0, length(ind)-1)
            theta = Gamma.sens/(1+Gamma.sens)
            
            for(j in 1:(length(ind)-1))
            {
              mu[j] = (2*theta-1)*abs(eta) + taubar - (2*theta-1)*(theta*abs(taubar + abs(eta)) + (1-theta)*abs(taubar - abs(eta)))
              sigma2[j] = theta*(1-theta)*(2*abs(eta) - (2*theta-1)*(abs(taubar + abs(eta))-  abs(taubar - abs(eta))))^2
            }
            mu[abs(mu) < 1e-8] = 0
            sigma2[sigma2 < 1e-8] = 0
            
            
            PM[symmgroup[ind.jj]] = mu*(max.e) - mu*(!max.e)
            PV[symmgroup[ind.jj]] = (sigma2)
            
            Diff[symmgroup[ind.jj]] = sum(outtreat -outcontrol - (null*(dosetreat - dosecontrol)))
            Diff2[symmgroup[count]] = sum(((dosetreat - dosecontrol)))
            count = count+1
          }
        }
      }
      
      
      values[(N.vars+1):(2*N.vars)] = Diff
      values[(2*N.vars+1):(3*N.vars)] = Diff2
      values[(3*N.vars+1):(4*N.vars+1)] = c(-PM, 1)
      b[nosymm+1] = 0
      b[nosymm+2] = 1
      b[nosymm+3] = 0
      alpha.opt = alpha
      if(alternative != "two.sided")
      {
        alpha.opt = 2*alpha
      }
      
      const.dir = c(rep("=", nosymm+1), ">=", "=")
      model = list()
      # if(Gamma.sens==1)
      # {
      #   model$A = sparseMatrix(row.ind[1:(3*N.vars)], col.ind[1:(3*N.vars)], x=values[1:(3*N.vars)])
      #   model$obj = c(PV)
      #   model$sense = const.dir[1:(nosymm+2)]
      #   model$rhs = b[1:(nosymm+2)]
      #   model$vtype = c(rep("I", N.vars))
      #   if(continuous.relax == T){model$vtype = c(rep("C", N.vars))}
      #   
      #   
      #   model$modelsense = "max"
      #   
      #   
      #   solm = gurobi(model, params = list(OutputFlag = 0))
      #   zed = (TE.est/sqrt(solm$objval))
      #   kappavec[ee] = (TE.est)^2 - qchisq(1-alpha.opt, 1)*solm$objval
      #   SE = sqrt(solm$objval)
      #   x = solm$x[1:N.vars]
      #   tstat = zed
      #   pval = 0
      #   if(alternative == "two.sided")
      #   {
      #     pval = 2*pnorm(-abs(tstat))
      #   }
      #   if(alternative == "greater")
      #   {
      #     pval = 1 - pnorm((tstat))
      #   }
      #   if(alternative == "less")
      #   {
      #     pval = pnorm((tstat))
      #   }
      #   Reject = (pval < alpha)
      # }
      # if(Gamma.sens != 1)
      # {
        diff = 200
        kappa = qchisq(1-alpha.opt, 1)
        while(diff > 1e-8)
        {
          th = Gamma.sens/(1+Gamma.sens)
          TE.est.Gamma = TE.est - sum((2*th-1)*abs(wTE.per.strat))
          Plin = -2*TE.est.Gamma*PM - kappa*PV 
          rowind.q =  1:(N.vars+1)
          colind.q = 1:(N.vars+1)
          values.q = c(rep(0, N.vars),1)
          Q = sparseMatrix(rowind.q, colind.q, x=values.q)
          model$A = sparseMatrix(row.ind, col.ind, x=values)
          model$obj = c(Plin,0)
          model$Q = Q
          model$sense = const.dir
          model$rhs = b
          
          
          model$vtype = c(rep("I", N.vars), "C")
          if(continuous.relax == T){model$vtype = c(rep("C", N.vars), "C")}
          model$lb = c(rep(0, N.vars), -Inf)
          
          
          model$modelsense = "min"
          
          
          solm = gurobi(model, params = list(OutputFlag = 0))
          x = solm$x[1:N.vars]
          kappa.new = (TE.est.Gamma - sum(PM*x))^2/sum(PV*x)
          diff = abs(kappa.new - kappa)
          pval = 0
          if(PVAL == F)
          {
            diff = 0
            Reject = (kappa.new > kappa)
          }
          kappavec[ee] = (TE.est.Gamma - sum(PM*x))^2 - qchisq(1-alpha.opt, 1)*sum(PV*x)
          kappa = kappa.new
        }
        zed = sqrt((TE.est.Gamma - sum(PM*x))^2/sum(PV*x))
        
        if(alternative == "less")
        {
          zed = -zed
        }
        zscore[ee] = zed
        tstat = zed
        
        
        if(alternative == "two.sided")
        {
          pval = 2*pnorm(-abs(tstat))
        }
        if(alternative == "greater")
        {
          pval = 1 - pnorm((tstat))
        }
        if(alternative == "less")
        {
          pval = pnorm((tstat))
        }
        
        
        
        
        
        if(sign(TE.est.Gamma - sum(PM*x))!=sign(TE.est))
        {
          Reject = F
          pval = 0.5
          kappavec[ee] = -10
          
          if(alternative == "two.sided")
          {
            pval = 1
            kappavec[ee] = -10
          }
        }
        
        if(alternative == "greater" & sign(TE.est.Gamma - sum(PM*x)) < 0)
        {
          pval = .5
          kappavec[ee] = -10
        }
        if(alternative == "less" & sign(TE.est.Gamma - sum(PM*x)) > 0)
        {
          pval = .5
          kappavec[ee] = -10
        }
        Reject = (pval < alpha)
        
      }
      pvalvec[ee] = pval
      Rejectvec[ee] = Reject  
    #}
  }
  if(ER == T)
  {
    mult.all = (mult.101*(mult.101+1)/2)*(mult.111*(mult.111+1)/2)*(mult.010*(mult.010+1)/2)*(mult.000*(mult.000+1)/2)*(mult.100*(mult.100+1)/2)*(mult.110*(mult.110+1)/2)*(mult.011*(mult.011+1)/2)*(mult.001*(mult.001+1)/2)
    if(MO == T)
    {
      mult.all = (mult.101*(mult.101+1)/2)*(mult.111*(mult.111+1)/2)*(mult.010*(mult.010+1)/2)*(mult.000*(mult.000+1)/2)
      if(DE == "nonpositive")
      {
        mult.all = (mult.111)*(mult.101*(mult.101+1)/2)*(mult.000)*(mult.010*(mult.010+1)/2)
      }
      if(DE=="nonnegative")
      {
        mult.all = (mult.101)*(mult.111*(mult.111+1)/2)*(mult.010)*(mult.000*(mult.000+1)/2)
      }
    }
    if(MO == F)
    {
      mult.all = (mult.101*(mult.101+1)/2)*(mult.111*(mult.111+1)/2)*(mult.010*(mult.010+1)/2)*(mult.000*(mult.000+1)/2)*(mult.100*(mult.100+1)/2)*(mult.110*(mult.110+1)/2)*(mult.011*(mult.011+1)/2)*(mult.001*(mult.001+1)/2)
      if(DE == "nonpositive")
      {
        mult.all = (mult.111)*(mult.101*(mult.101+1)/2)*(mult.000)*(mult.010*(mult.010+1)/2)*mult.001*(mult.011*(mult.011+1)/2)*mult.110*(mult.100*(mult.100+1)/2)
      }
      if(DE=="nonnegative")
      {
        mult.all = (mult.101)*(mult.111*(mult.111+1)/2)*(mult.010)*(mult.000*(mult.000+1)/2)*mult.011*(mult.001*(mult.001+1)/2)*mult.100*(mult.110*(mult.110+1)/2)
      }
    }
    
    ns.type = tapply(ns, num.id, mean)
    N.vars = sum(mult.all*(ns.type-1))
    index.symm = rep(1:nosymm, mult.all*(ns.type-1))
    n.types = (mult.all)
    mult.po= mult.all
    mult.do = mult.all
    n.po = mult.po
    n.do= mult.do
    n.per = cc
    Diff = rep(0, N.vars)
    Diff2 = Diff
    #V.list = vector("list", N.vars)
    row.ind = rep(0, 3*N.vars+1)
    col.ind = row.ind
    values = row.ind
    b = rep(0, nosymm + 2)
    for(kk in 1:nosymm)
    {
      row.ind[which(index.symm==kk)] = rep(kk, (ns.type[kk]-1)*n.types[kk])
      col.ind[which(index.symm==kk)] = which(index.symm==kk)  
      values[which(index.symm==kk)] = rep(1, (ns.type[kk]-1)*(n.types[kk]))
      b[kk] = n.per[kk]
    }
    row.ind[(N.vars+1):(2*N.vars)] = rep(nosymm + 1, N.vars)
    col.ind[(N.vars+1):(2*N.vars)] = 1:N.vars
    row.ind[(2*N.vars+1):(3*N.vars)] = rep(nosymm + 2, N.vars)
    col.ind[(2*N.vars+1):(3*N.vars)] = 1:N.vars
    row.ind[(3*N.vars+1):(4*N.vars+1)] = rep(nosymm + 3, N.vars+1)
    col.ind[(3*N.vars+1):(4*N.vars+1)] = 1:(N.vars+1)
    PM = rep(0, N.vars)
    PV = rep(0, N.vars)
    
    
    zscore = rep(0, length(Gamma.vec))
    
    for(ee in 1:length(Gamma.vec))
    {
      Gamma.sens = Gamma.vec[ee]
      
      for(kk in 1:nosymm)
      {
        i = which(num.id==kk)[1]
        symmgroup = which(index.symm == kk)
        ind = which(index==i)
        treatstrat = treatment[ind]
        outstrat = outcome[ind]
        dosestrat = dose[ind]
        dosesymm = c(rep(0, sum(treatstrat==F&dosestrat==F)), rep(1, sum(treatstrat==F& dosestrat==T)), rep(0, sum(treatstrat==T&dosestrat==F)), rep(1, sum(treatstrat==T& dosestrat==T)))
        
        outsymm = c(sort(outstrat[treatstrat==F& dosestrat==F]), sort(outstrat[treatstrat==F& dosestrat==T]), sort(outstrat[treatstrat==T&dosestrat == F]), sort(outstrat[treatstrat==T&dosestrat == T]))
        PO = matrix(0, n.po[kk], ns[i])
        DO = matrix(0, n.po[kk],ns[i])
        
        mt.011 = mult.011[kk]
        mt.001 = mult.001[kk]
        mt.000 = mult.000[kk]
        mt.010 = mult.010[kk]
        mt.110 = mult.110[kk]
        mt.100 = mult.100[kk]
        mt.111 = mult.111[kk]
        mt.101 = mult.101[kk]
        p.011 = mt.011
        p.001 = mt.001
        p.000 = mt.000
        p.010 = mt.010
        p.110 = mt.110
        p.100 = mt.100
        p.111 = mt.111
        p.101 = mt.101
        
        D.000 = matrix(1, mt.000-1, mt.000-1)
        D.000[lower.tri(D.000)] = 0
        D.000 = rbind(D.000, c(rep(0, mt.000-1)) )
        D.010 = matrix(1, mt.010-1, mt.010-1)
        D.010[lower.tri(D.010)] = 0
        D.010 = rbind(D.010, c(rep(0, mt.010-1)) )
        D.001 = matrix(1, mt.001-1, mt.001-1)
        D.001[lower.tri(D.001)] = 0
        D.001 = rbind(D.001, c(rep(0, mt.001-1)) )
        D.011 = matrix(1, mt.011-1, mt.011-1)
        D.011[lower.tri(D.011)] = 0
        D.011 = rbind(D.011, c(rep(0, mt.011-1)))
        D.100 = matrix(1, mt.100-1, mt.100-1)
        D.100[lower.tri(D.100)] = 0
        D.100 = rbind(D.100, c(rep(0, mt.100-1)) )
        D.110 = matrix(1, mt.110-1, mt.110-1)
        D.110[lower.tri(D.110)] = 0
        D.110 = rbind(D.110, c(rep(0, mt.110-1)) )
        D.101 = matrix(1, mt.101-1, mt.101-1)
        D.101[lower.tri(D.101)] = 0
        D.101 = rbind(D.101, c(rep(0, mt.101-1)) )
        
        D.111 = matrix(1, mt.111-1, mt.111-1)
        D.111[lower.tri(D.111)] = 0
        D.111 = rbind(D.111, c(rep(0, mt.111-1)))
        
        if(MO == T)
        {
          D.100 = matrix(0, 1, mt.100-1)
          p.100 = 1
          D.110 = matrix(0, 1, mt.110-1)
          p.110 = 1
          D.011 = matrix(1, 1, mt.011-1)
          p.011 = 1
          D.001 = matrix(1, 1, mt.001-1)
          p.001 = 1
        }
        count = 1
        for(ll in 1:p.000)
        {
          for(mm in 1:p.010)
          {
            for(uu in 1:p.001)
            {
              for(vv in 1:p.011)
              {
                for(ww in 1:p.100)
                {
                  for(xx in 1:p.110)
                  {
                    for(yy in 1:p.101)
                    {
                      for(zz in 1:p.111)	
                      {
                        T.000 = matrix(1, p.000-ll, mt.000-ll)
                        T.000[lower.tri(T.000)] = 0
                        T.000 = rbind(T.000, c(rep(0, mt.000-ll)))
                        F.000 = matrix(0, (p.000-ll+1), (ll-1))
                        T.000 = cbind(F.000, T.000)
                        
                        T.010 = matrix(1, p.010-mm, mt.010-mm)
                        T.010[lower.tri(T.010)] = 0
                        T.010 = rbind(T.010, c(rep(0, mt.010-mm)))
                        F.010 = matrix(1, (p.010-mm+1), (mm-1))
                        T.010 = cbind(F.010, T.010)
                        
                        T.001 = matrix(1, uu-1, uu-1)
                        T.001[lower.tri(T.001)] = 0
                        T.001 = rbind(T.001, c(rep(0, uu-1)))
                        F.001 = matrix(0, (uu), (mt.001-uu))
                        T.001 = cbind(F.001, T.001)
                        if(MO == T)
                        {
                          T.001 = matrix(0, 1, mt.001-1)
                        }
                        
                        T.011 = matrix(1, vv-1, vv-1)
                        T.011[lower.tri(T.011)] = 0
                        T.011 = rbind(T.011, c(rep(0, vv-1)))
                        F.011 = matrix(1, (vv), (mt.011-vv))
                        T.011 = cbind(F.011, T.011)
                        if(MO == T)
                        {
                          T.011 = matrix(1, 1, mt.011-1)
                        }
                        
                        
                        T.100 = matrix(1, p.100-ww, mt.100-ww)
                        T.100[lower.tri(T.100)] = 0
                        T.100 = rbind(T.100, c(rep(0, mt.100-ww)))
                        F.100 = matrix(0, (p.100-ww+1), (ww-1))
                        T.100 = cbind(F.100, T.100)
                        if(MO == T)
                        {
                          T.100 = matrix(0, 1,mt.100-1)
                        }
                        
                        T.110 = matrix(1, p.110-xx, mt.110-xx)
                        T.110[lower.tri(T.110)] = 0
                        T.110 = rbind(T.110, c(rep(0, mt.110-xx)))
                        F.110 = matrix(1, (p.110-xx+1), (xx-1))
                        T.110 = cbind(F.110, T.110)
                        if(MO == T)
                        {
                          T.110 = matrix(1, 1, mt.110-1)
                        }
                        
                        T.101 = matrix(1, yy-1, yy-1)
                        T.101[lower.tri(T.101)] = 0
                        T.101 = rbind(T.101, c(rep(0, yy-1)))
                        F.101 = matrix(0, (yy), (mt.101-yy))
                        T.101 = cbind(F.101, T.101)
                        
                        T.111 = matrix(1, zz-1, zz-1)
                        T.111[lower.tri(T.111)] = 0
                        T.111 = rbind(T.111, c(rep(0, zz-1)))
                        F.111 = matrix(1, (zz), (mt.111-zz))
                        T.111 = cbind(F.111, T.111) 
                        
                        
                        
                        
                        
                        if(DE == "nonnegative")
                        {
                          T.100 = matrix(0, 1, mt.100-1)
                          T.101 = matrix(0, 1, mt.101-1)
                          T.011 = matrix(1, 1, mt.011-1)
                          T.010 = matrix(1, 1, mt.010-1)
                        }
                        if(DE == "nonpositive")
                        {
                          T.111 = matrix(1, 1, mt.111-1)
                          T.110 = matrix(1, 1, mt.110-1)
                          T.001 = matrix(0, 1, mt.001-1)
                          T.000 = matrix(0, 1, mt.000-1)
                          
                        }
                        
                        
                        
                        dosevec = c(D.000[ll,], D.010[mm,], D.001[uu,], D.011[vv,], D.100[ww,], D.110[xx,], D.101[yy,], D.111[zz,])
                        dosevec = dosevec[!is.na(dosevec)]
                        for(lll in 1:nrow(T.000))
                        {
                          for(mmm in 1:nrow(T.010))
                          {
                            for(uuu in 1:nrow(T.001))
                            {
                              for(vvv in 1:nrow(T.011))
                              {
                                for(www in 1:nrow(T.100))
                                {
                                  for(xxx in 1:nrow(T.110))
                                  {
                                    for(yyy in 1:nrow(T.101))
                                    {
                                      for(zzz in 1:nrow(T.111))	
                                      {
                                        DO[count,] = dosevec
                                        tempvec = c(T.000[lll,], T.010[mmm,], T.001[uuu,], T.011[vvv,], T.100[www,], T.110[xxx,], T.101[yyy,], T.111[zzz,])
                                        tempvec = tempvec[!is.na(tempvec)]
                                        PO[count,] = tempvec
                                        count = count+1
                                      }}}}}}}}}}}}}}}}             
        
        
        
        count = 1
        
        for(jj in 1:n.po[kk])
        {
          
          ind.jj = (((jj-1)*(ns[i]-1))+1):(jj*(ns[i]-1))
          po.symm = PO[jj,]
          do.symm = DO[jj,]
          treatsymm = c(rep(F, ns[i]-ms[i]), rep(T, ms[i]))
          outcontrol = outsymm*(1-treatsymm) + po.symm*(treatsymm)
          outtreat = outsymm*(treatsymm) + po.symm*(1-treatsymm) 
          dosecontrol = dosesymm*(1-treatsymm) + do.symm*(treatsymm)
          dosetreat = dosesymm*(treatsymm) + do.symm*(1-treatsymm) 
          sum.cont = sum(outcontrol - null*dosecontrol)/(ns[i]-1)
          Q = (outtreat + outcontrol/(ns[i]-1) - null*dosetreat - null*dosecontrol/(ns[i]-1)- sum.cont)*ns[i]
          if(sum(treatstrat)>1)
          {
            sum.cont = sum(outtreat - null*dosetreat)/(ns[i]-1)
            Q = -(outtreat/(ns[i]-1) + outcontrol - null*dosetreat/(ns[i]-1) - null*dosecontrol- sum.cont)*ns[i]
            
          }
          
          
          
          qi = Q*max.e - Q*(!max.e)
          ord = order(qi)
          qi.sort = sort(qi)
          
          
          eta = diff(outcontrol+outtreat - null*(dosecontrol+dosetreat))/2
          taubar = mean(outtreat-outcontrol - null*(dosetreat-dosecontrol))
          
          mu = rep(0, length(ind)-1)
          sigma2 = rep(0, length(ind)-1)
          theta = Gamma.sens/(1+Gamma.sens)
          
          for(j in 1:(length(ind)-1))
          {
            mu[j] = (2*theta-1)*abs(eta) + taubar - (2*theta-1)*(theta*abs(taubar + abs(eta)) + (1-theta)*abs(taubar - abs(eta)))
            sigma2[j] = theta*(1-theta)*(2*abs(eta) - (2*theta-1)*(abs(taubar + abs(eta))-  abs(taubar - abs(eta))))^2
          }
          mu[abs(mu) < 1e-8] = 0
          sigma2[sigma2 < 1e-8] = 0
          
          
          PM[symmgroup[ind.jj]] = mu*(max.e) - mu*(!max.e)
          PV[symmgroup[ind.jj]] = (sigma2)
          Diff[symmgroup[ind.jj]] = sum(outtreat -outcontrol - (null*(dosetreat - dosecontrol)))
          
          Diff2[symmgroup[count]] = sum(((dosetreat - dosecontrol)))
          count = count+1
          
        }
      }
      
      
      values[(N.vars+1):(2*N.vars)] = Diff
      values[(2*N.vars+1):(3*N.vars)] = Diff2
      values[(3*N.vars+1):(4*N.vars+1)] = c(-PM, 1)
      b[nosymm+1] = 0
      b[nosymm+2] = 1
      b[nosymm+3] = 0  
      alpha.opt = alpha
      if(alternative != "two.sided")
      {
        alpha.opt = 2*alpha
      }
      
      const.dir = c(rep("=", nosymm+1), ">=", "=")
      model = list()
      if(Gamma.sens==1)
      {
        model$A = sparseMatrix(row.ind[1:(3*N.vars)], col.ind[1:(3*N.vars)], x=values[1:(3*N.vars)])
        model$obj = c(PV)
        model$sense = const.dir[1:(nosymm+2)]
        model$rhs = b[1:(nosymm+2)]
        model$vtype = c(rep("I", N.vars))
        if(continuous.relax == T){model$vtype = c(rep("C", N.vars))}


        model$modelsense = "max"


        solm = gurobi(model, params = list(OutputFlag = 0))
        zed = (TE.est/sqrt(solm$objval))
        x = solm$x[1:N.vars]
        SE = sqrt(solm$objval)
        kappavec[ee] = (TE.est)^2 - qchisq(1-alpha.opt, 1)*solm$objval
        tstat = zed
        pval = 0
        if(alternative == "two.sided")
        {
          pval = 2*pnorm(-abs(tstat))
        }
        if(alternative == "greater")
        {
          pval = 1 - pnorm((tstat))
        }
        if(alternative == "less")
        {
          pval = pnorm((tstat))
        }
        Reject = (pval < alpha)
      }
      if(Gamma.sens != 1)
      {
        diff = 200
        kappa = qchisq(1-alpha.opt, 1)
        while(diff > 1e-8)
        {
          th = Gamma.sens/(1+Gamma.sens)
          TE.est.Gamma = TE.est - sum((2*th-1)*abs(wTE.per.strat))
          Plin = -2*TE.est.Gamma*PM - kappa*PV 
          rowind.q =  1:(N.vars+1)
          colind.q = 1:(N.vars+1)
          values.q = c(rep(0, N.vars),1)
          Q = sparseMatrix(rowind.q, colind.q, x=values.q)
          model$A = sparseMatrix(row.ind, col.ind, x=values)
          model$obj = c(Plin,0)
          model$Q = Q
          model$sense = const.dir
          model$rhs = b
          
          
          model$vtype = c(rep("I", N.vars), "C")
          if(continuous.relax == T){model$vtype = c(rep("C", N.vars+1))}
          model$lb = c(rep(0, N.vars), -Inf)
          
          
          model$modelsense = "min"
          
          
          solm = gurobi(model, params = list(OutputFlag = 0))
          x = solm$x[1:N.vars]
          kappa.new = (TE.est.Gamma - sum(PM*x))^2/sum(PV*x)
          kappavec[ee] = (TE.est.Gamma - sum(PM*x))^2 - qchisq(1-alpha.opt, 1)*sum(PV*x)
          diff = abs(kappa.new - kappa)
          pval = 0
          if(PVAL == F)
          {
            diff = 0
            Reject = kappa.new > kappa
          }
          kappa = kappa.new
        }
        zed = sqrt((TE.est.Gamma - sum(PM*x))^2/sum(PV*x))
        
        if(alternative == "less")
        {
          zed = -zed
        }
        zscore[ee] = zed
        tstat = zed
        
        
        if(alternative == "two.sided")
        {
          pval = 2*pnorm(-abs(tstat))
        }
        if(alternative == "greater")
        {
          pval = 1 - pnorm((tstat))
        }
        if(alternative == "less")
        {
          pval = pnorm((tstat))
        }
       
        
        if(sign(TE.est.Gamma - sum(PM*x))!=sign(TE.est))
        {
          Reject = F
          pval = 0.5
          kappavec[ee] = -10
          
          if(alternative == "two.sided")
          {
            pval = 1
            kappavec[ee] = -10
          }
        }
        
        if(alternative == "greater" & sign(TE.est.Gamma - sum(PM*x)) < 0)
        {
          pval = .5
          kappavec[ee] = -10
        }
        if(alternative == "less" & sign(TE.est.Gamma - sum(PM*x)) > 0)
        {
          pval = .5
          kappavec[ee] = -10
        }
        Reject = (pval < alpha)
      }
      pvalvec[ee] = pval
      Rejectvec[ee] = Reject 
    } 
}
  if(PVAL == F)
  {
    return(list(Gamma.vec = Gamma.vec, Reject = Rejectvec, lambdahat = lambdahat, kappa = kappavec))
  }
  if(PVAL == T)
  {
    return(list(Gamma.vec = Gamma.vec, pval = pvalvec, lambdahat = lambdahat))
  }
  
  
}  

sensEffectRatio2  = function(Gamma.vec, index, treatment, outcome, dose, null = 0, DE = "both", MO = T, ER = T, alternative = "two.sided", alpha = 0.05, continuous.relax = F)
{
  sensEffectRatio(index, treatment, outcome, dose, null, DE, MO, ER, alternative, alpha, Gamma.vec = Gamma.vec, calculate.pval = F, continuous.relax)$kappa
}