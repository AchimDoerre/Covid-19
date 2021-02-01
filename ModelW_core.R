ModelW_core <- function(days,
                        td=1/4, NR=10,
                        Icum0 = rep(0, 18),
                        S0, E0, I0, R0, D0,
                        CM = CM,
                        w = 0.1*td,
                        mitfactor = matrix(0, 18, 18),
                        general_compliance_drift = rep(1, length(days)/td),
                        plot_overall_course = FALSE){
  
  daysmax = length(days)
  tmax = daysmax/td
  
  
  SR = matrix(0, NR, tmax)
  ER = matrix(0, NR, tmax)
  IR = matrix(0, NR, tmax)
  RR = matrix(0, NR, tmax)
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
  
  
  D1R = matrix(0, NR, tmax)
  D2R = matrix(0, NR, tmax)
  D3R = matrix(0, NR, tmax)
  D4R = matrix(0, NR, tmax)
  D5R = matrix(0, NR, tmax)
  D6R = matrix(0, NR, tmax)
  D7R = matrix(0, NR, tmax)
  D8R = matrix(0, NR, tmax)
  D9R = matrix(0, NR, tmax)
  D10R = matrix(0, NR, tmax)
  D11R = matrix(0, NR, tmax)
  D12R = matrix(0, NR, tmax)
  D13R = matrix(0, NR, tmax)
  D14R = matrix(0, NR, tmax)
  D15R = matrix(0, NR, tmax)
  D16R = matrix(0, NR, tmax)
  D17R = matrix(0, NR, tmax)
  D18R = matrix(0, NR, tmax)
  
  
  P = matrix(0, NR, 3)
  
  
  PR = read.table("ModelW_parameters.txt", header=TRUE)
  
  
  for (nr in 1:NR){
    print(nr)
    
    
    gammanr = PR[nr,"gammanr"]*td
    epsilonnr = PR[nr,"epsilonnr"]*td
    mannr = PR[nr,"mannr"]
    rnr = PR[nr,"rnr"]
    
    
    P[nr,] = c(gammanr/td, mannr, rnr)
    colnames(P) = c("gamma", "man", "r")
    
    
    
    #initialization
    S = rep(0, tmax)
    I = rep(0, tmax)
    E = rep(0, tmax)
    R = rep(0, tmax)
    D = rep(0, tmax)

    
    Sa = matrix(0, tmax, ngroups)
    Ea = matrix(0, tmax, ngroups)
    Ia = matrix(0, tmax, ngroups)
    Ra = matrix(0, tmax, ngroups)
    Da = matrix(0, tmax, ngroups)
    
    
    #initial values
    Sa[1, ] = S0
    Ea[1, ] = E0
    Ia[1, ] = I0
    Ra[1, ] = R0
    Da[1, ] = D0
    
    
    
    S[1] = sum(Sa[1,])
    E[1] = sum(Ea[1,])
    I[1] = sum(Ia[1,])
    R[1] = sum(Ra[1,])
    D[1] = sum(Da[1,])

    
    SR[nr, 1] = S[1]
    ER[nr, 1] = E[1]
    IR[nr, 1] = I[1]
    RR[nr, 1] = R[1]
    DR[nr, 1] = D[1]
    
    

    
    tseq = 2:tmax
    for (t in tseq){
      betaM = general_compliance_drift[t]*w*(mitfactor*CM)
      #print(t)
      for (a in 1:ngroups){
        Ea[t,a] = Ea[t-1,a] - mannr*epsilonnr*Ea[t-1,a] - (1-mannr)*epsilonnr*Ea[t-1,a]
        Sa[t,a] = Sa[t-1,a]
        for (b in 1:ngroups){
          betaab = betaM[a,b]
          free_Ib = (1-rnr)*(1-ha[b])*Ia[t-1,b]
          new_exposedab = betaab*Sa[t-1,a]*free_Ib/Sa[1,b]
          Sa[t,a] = Sa[t,a] - new_exposedab
          Ea[t,a] = Ea[t,a] + new_exposedab
        }
        Ia[t,a] = Ia[t-1,a] + mannr*epsilonnr*Ea[t-1,a] - gammanr*Ia[t-1,a]
        Ra[t,a] = Ra[t-1,a] + (1-mannr)*epsilonnr*Ea[t-1,a] + (1-dI[a])*gammanr*Ia[t-1,a]
        Da[t,a] = Da[t-1,a] + dI[a]*gammanr*Ia[t-1,a]
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
      
      D1R[nr, t] = Da[t,1]
      D2R[nr, t] = Da[t,2]
      D3R[nr, t] = Da[t,3]
      D4R[nr, t] = Da[t,4]
      D5R[nr, t] = Da[t,5]
      D6R[nr, t] = Da[t,6]
      D7R[nr, t] = Da[t,7]
      D8R[nr, t] = Da[t,8]
      D9R[nr, t] = Da[t,9]
      D10R[nr, t] = Da[t,10]
      D11R[nr, t] = Da[t,11]
      D12R[nr, t] = Da[t,12]
      D13R[nr, t] = Da[t,13]
      D14R[nr, t] = Da[t,14]
      D15R[nr, t] = Da[t,15]
      D16R[nr, t] = Da[t,16]
      D17R[nr, t] = Da[t,17]
      D18R[nr, t] = Da[t,18]
      
      
      S[t] = sum(Sa[t,])
      E[t] = sum(Ea[t,])
      I[t] = sum(Ia[t,])
      R[t] = sum(Ra[t,])
      D[t] = sum(Da[t,])
      
      
      SR[nr, t] = S[t]
      ER[nr, t] = E[t]
      IR[nr, t] = I[t]
      RR[nr, t] = R[t]
      DR[nr, t] = D[t]
      
    }
    
    
    #----------------------------------------------
    # Tode
    #----------------------------------------------
    
    deatha = Da[tmax,]
    death_total = sum(deatha)
    
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
  
  
  D1R[, 1] = Da[1,1]
  D2R[, 1] = Da[1,2]
  D3R[, 1] = Da[1,3]
  D4R[, 1] = Da[1,4]
  D5R[, 1] = Da[1,5]
  D6R[, 1] = Da[1,6]
  D7R[, 1] = Da[1,7]
  D8R[, 1] = Da[1,8]
  D9R[, 1] = Da[1,9]
  D10R[, 1] = Da[1,10]
  D11R[, 1] = Da[1,11]
  D12R[, 1] = Da[1,12]
  D13R[, 1] = Da[1,13]
  D14R[, 1] = Da[1,14]
  D15R[, 1] = Da[1,15]
  D16R[, 1] = Da[1,16]
  D17R[, 1] = Da[1,17]
  D18R[, 1] = Da[1,18]
  
  
  SR[, 1] = S[1]
  ER[, 1] = E[1]
  IR[, 1] = I[1]
  RR[, 1] = R[1]
  DR[, 1] = D[1]
  
  
  
  Sm = colMeans(SR)
  Em = colMeans(ER)
  Im = colMeans(IR)
  Rm = colMeans(RR)
  Dm = colMeans(DR)

  
  forecast_horizon = daysmax
  forecast = matrix(0, forecast_horizon, 3)
  tseqr = seq(1/td, tmax, 1/td)
  forecast[,1] = sum(Icum0) + Im[tseqr] + (Rm[tseqr] - sum(R0)) + (Dm[tseqr] - sum(D0))
  forecast[,2] = Im[tseqr]
  forecast[,3] = Dm[tseqr]
  colnames(forecast) = c("Icum", "Iactive", "Deaths")
  rownames(forecast) = days
  
  forecast = round(forecast, 0)
  
  
  if (plot_overall_course){
    options(scipen=5)
    par(mar=c(4,4,2,0.1))
    plot(c(1, tseq), Em, typ="l", col=linecolE, main="Overall (mean) course", xlab="time", ylim=c(1, max(Em, Im, Dm)))
    lines(c(1, tseq), Im, col=linecolI)
    lines(c(1, tseq), Dm, col=linecolD)
    legend("left", legend=c("E", "I", "D"), col=c(linecolE, linecolI, linecolD), lwd=1)
  }
  
  result = list(Sm = Sm, Em = Em, Im = Im, Rm = Rm, Dm = Dm,
                SR = SR, ER = ER, IR = IR, RR = RR, DR = DR,
                I1R = I1R, I2R = I2R, I3R = I3R, I4R = I4R, I5R = I5R, I6R = I6R,
                I7R = I7R, I8R = I8R, I9R = I9R, I10R = I10R, I11R = I11R, I12R = I12R,
                I13R = I13R, I14R = I14R, I15R = I15R, I16R = I16R, I17R = I17R, I18R = I18R,
                D1R = D1R, D2R = D2R, D3R = D3R, D4R = D4R, D5R = D5R, D6R = D6R,
                D7R = D7R, D8R = D8R, D9R = D9R, D10R = D10R, D11R = D11R, D12R = D12R,
                D13R = D13R, D14R = D14R, D15R = D15R, D16R = D16R, D17R = D17R, D18R = D18R,
                forecast = forecast)
  return(result)
}
