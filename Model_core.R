ModelU_core <- function(daysmax = 30,
                        td=1/4, NR=10,
                        S0, E0, I0, R0, D0,
                        betaM,
                        epsilonrange = c(0.19, 0.21),
                        gammarange = c(1/8, 1/8),
                        manrange = c(0.6, 0.8),
                        rrange = c(0.2, 0.3),
                        w = 0.10,
                        mitfactor = matrix(1, ncol(betaM), ncol(betaM))){
  
  tmax = daysmax/td
  betaMR = w*mitfactor*betaM
  
  SR = matrix(0, NR, tmax)
  ER = matrix(0, NR, tmax)
  IR = matrix(0, NR, tmax)
  RR = matrix(0, NR, tmax)
  HR = matrix(0, NR, tmax)
  ICUR = matrix(0, NR, tmax)
  DR = matrix(0, NR, tmax)
  
  
  I1R = matrix(0, NR, tmax)
  I2R = matrix(0, NR, tmax)
  I3R = matrix(0, NR, tmax)
  I4R = matrix(0, NR, tmax)
  I5R = matrix(0, NR, tmax)
  I6R = matrix(0, NR, tmax)
  I7R = matrix(0, NR, tmax)
  I8R = matrix(0, NR, tmax)
  I9R = matrix(0, NR, tmax)
  I10R = matrix(0, NR, tmax)
  I11R = matrix(0, NR, tmax)
  I12R = matrix(0, NR, tmax)
  I13R = matrix(0, NR, tmax)
  I14R = matrix(0, NR, tmax)
  I15R = matrix(0, NR, tmax)
  I16R = matrix(0, NR, tmax)
  I17R = matrix(0, NR, tmax)
  I18R = matrix(0, NR, tmax)
  
  
  #M contains the basic forecasting information for every random iteration
  M = matrix(0, NR, 4)
  
  #PR stores the randomly drawn parameter values
  PR = matrix(0, NR, 3)
  colnames(PR) = c("gamma", "man", "r")
  
  
  #loop for random trajectories
  for (nr in 1:NR){
    print(nr)
    
    
    #random draw of parameters from given range
    gammanr = runif(1, min(gammarange), max(gammarange))*td
    epsilonr = runif(1, min(epsilonrange), max(epsilonrange))*td
    mannr = runif(1, min(manrange), max(manrange))
    rnr = runif(1, min(rrange), max(rrange))
    
    
    PR[nr,] = c(gammanr/td, mannr, rnr)

    
    
    #initialization
    S = rep(0, tmax)
    I = rep(0, tmax)
    E = rep(0, tmax)
    R = rep(0, tmax)
    H = rep(0, tmax)
    D = rep(0, tmax)
    ICU = rep(0, tmax)
    
    
    
    Sa = matrix(0, tmax, ngroups)
    Ea = matrix(0, tmax, ngroups)
    Ia = matrix(0, tmax, ngroups)
    Ra = matrix(0, tmax, ngroups)
    Ha = matrix(0, tmax, ngroups)
    Da = matrix(0, tmax, ngroups)
    
    
    #initial values
    Sa[1, ] = S0
    Ea[1, ] = E0
    Ia[1, ] = I0
    for (a in 1:ngroups){
      Ha[1,a] = ha[a]*Ia[1,a]
    }
    Ra[1, ] = R0
    Da[1, ] = D0
    
    
    
    S[1] = sum(Sa[1,])
    E[1] = sum(Ea[1,])
    I[1] = sum(Ia[1,])
    R[1] = sum(Ra[1,])
    H[1] = sum(Ha[1,])
    D[1] = sum(Da[1,])
    ICU[1] = sum(ICUp*Ha[1,])
    
    
    SR[nr, 1] = S[1]
    ER[nr, 1] = E[1]
    IR[nr, 1] = I[1]
    RR[nr, 1] = R[1]
    HR[nr, 1] = H[1]
    ICUR[nr, 1] = ICU[1]
    DR[nr, 1] = D[1]
    
    
    
    
    tseq = 2:tmax
    for (t in tseq){
      for (a in 1:ngroups){
        Ea[t,a] = Ea[t-1,a] - mannr*epsilonr*Ea[t-1,a] - (1-mannr)*epsilonr*Ea[t-1,a]
        Sa[t,a] = Sa[t-1,a]
        for (b in 1:ngroups){
          betaab = betaMR[a,b]
          free_Ib = (1-rnr)*(1-ha[b])*Ia[t-1,b]
          new_exposedab = betaab*Sa[t-1,a]*free_Ib/Sa[1,b]
          Sa[t,a] = Sa[t,a] - new_exposedab
          Ea[t,a] = Ea[t,a] + new_exposedab
        }
        Ia[t,a] = Ia[t-1,a] + mannr*epsilonr*Ea[t-1,a] - gammanr*Ia[t-1,a]
        Ra[t,a] = Ra[t-1,a] + (1-mannr)*epsilonr*Ea[t-1,a] + (1-dI[a])*gammanr*Ia[t-1,a]
        Da[t,a] = Da[t-1,a] + dI[a]*gammanr*Ia[t-1,a]
        Ha[t,a] = ha[a]*Ia[t,a]
      }
      I1R[nr, t] = Ia[t,1]
      I2R[nr, t] = Ia[t,2]
      I3R[nr, t] = Ia[t,3]
      I4R[nr, t] = Ia[t,4]
      I5R[nr, t] = Ia[t,5]
      I6R[nr, t] = Ia[t,6]
      I7R[nr, t] = Ia[t,7]
      I8R[nr, t] = Ia[t,8]
      I9R[nr, t] = Ia[t,9]
      I10R[nr, t] = Ia[t,10]
      I11R[nr, t] = Ia[t,11]
      I12R[nr, t] = Ia[t,12]
      I13R[nr, t] = Ia[t,13]
      I14R[nr, t] = Ia[t,14]
      I15R[nr, t] = Ia[t,15]
      I16R[nr, t] = Ia[t,16]
      I17R[nr, t] = Ia[t,17]
      I18R[nr, t] = Ia[t,18]
      
      
      S[t] = sum(Sa[t,])
      E[t] = sum(Ea[t,])
      I[t] = sum(Ia[t,])
      R[t] = sum(Ra[t,])
      H[t] = sum(Ha[t,])
      ICU[t] = sum(ICUp*Ha[t,])
      D[t] = sum(Da[t,])
      
      
      SR[nr, t] = S[t]
      ER[nr, t] = E[t]
      IR[nr, t] = I[t]
      RR[nr, t] = R[t]
      HR[nr, t] = H[t]
      ICUR[nr, t] = ICU[t]
      DR[nr, t] = D[t]
      
    }
    
    

    #deaths
    deatha = Da[tmax,]
    death_total = sum(deatha)
    
    
    for (t in 1:tmax){
      if (I[t] == max(I)){
        Imax = t*td
      }
    }
    M[nr,] = c(nr, Imax, max(ICU), death_total)
  }
  
  
  I1R[, 1] = Ia[1,1]
  I2R[, 1] = Ia[1,2]
  I3R[, 1] = Ia[1,3]
  I4R[, 1] = Ia[1,4]
  I5R[, 1] = Ia[1,5]
  I6R[, 1] = Ia[1,6]
  I7R[, 1] = Ia[1,7]
  I8R[, 1] = Ia[1,8]
  I9R[, 1] = Ia[1,9]
  I10R[, 1] = Ia[1,10]
  I11R[, 1] = Ia[1,11]
  I12R[, 1] = Ia[1,12]
  I13R[, 1] = Ia[1,13]
  I14R[, 1] = Ia[1,14]
  I15R[, 1] = Ia[1,15]
  I16R[, 1] = Ia[1,16]
  I17R[, 1] = Ia[1,17]
  I18R[, 1] = Ia[1,18]
  
  
  SR[, 1] = S[1]
  ER[, 1] = E[1]
  IR[, 1] = I[1]
  RR[, 1] = R[1]
  HR[, 1] = H[1]
  ICUR[, 1] = ICU[1]
  DR[, 1] = D[1]
  
  
  
  M = round(M, 2)
  colnames(M) = c("Iteration", "Imax", "ICUmax", "Deaths")
  
  Sm = colMeans(SR)
  Em = colMeans(ER)
  Im = colMeans(IR)
  Rm = colMeans(RR)
  Hm = colMeans(HR)
  Dm = colMeans(DR)
  ICUm = colMeans(ICUR)
  
  forecast_horizon = daysmax
  forecast = matrix(0, forecast_horizon, 4)
  tseqr = seq(1/td, tmax, 1/td)
  forecast[,1] = (Im+Rm+Dm)[tseqr]
  forecast[,2] = Im[tseqr]
  forecast[,3] = ICUm[tseqr]
  forecast[,4] = Dm[tseqr]
  colnames(forecast) = c("Icum", "Iactive", "ICU", "Deaths")
  rownames(forecast) = 1:daysmax
  
  forecast = round(forecast, 2)
  
  
  plot(c(1, tseq), Rm, log='y', typ="l", col=rgb(0, 0.6, 0.0), xaxt="n", main="Overall course",
       ylim=c(1, max(Em, Im, Rm, Dm)))
  lines(c(1, tseq), Em, col=linecolE)
  lines(c(1, tseq), Hm, col=linecolH)
  lines(c(1, tseq), ICUm, col=linecolICU)
  lines(c(1, tseq), Im, col=linecolI)
  lines(c(1, tseq), Dm, col=linecolD)
  


  result = list(Sm = Sm, Em = Em, Im = Im, Rm = Rm, Hm = Hm, Dm = Dm, ICUm = ICUm,
                SR = SR, ER = ER, IR = IR, RR = RR, HR = HR, DR = DR, ICUR = ICUR,
                I1R = I1R, I2R = I2R, I3R = I3R, I4R = I4R, I5R = I5R, I6R = I6R,
                I7R = I7R, I8R = I8R, I9R = I9R, I10R = I10R, I11R = I11R, I12R = I12R,
                I13R = I13R, I14R = I14R, I15R = I15R, I16R = I16R, I17R = I17R, I18R = I18R,
                forecast = forecast,
                M = M,
                PR = PR)
  return(result)
}
