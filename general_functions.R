forecast_band <- function(ZR, linecol, shadecol, description){
  Zm = colMeans(ZR)
  tseq = 2:ncol(ZR)
  xseq = c(1, tseq, rev(tseq), 1)
  
  yseq1 = rep(0, ncol(ZR))
  yseq2 = rep(0, ncol(ZR))
  for (t in 1:length(yseq1)){
    yseq1[t] = quantile(ZR[,t], 0.10)
    yseq2[t] = quantile(ZR[,t], 0.90)
  }
  yseq = c(yseq1, rev(yseq2))
  
  if (logscale){
    plot(c(1, tseq*td), Zm, col="black", ylim=c(1, max(ZR)), log='y', typ="l", main=description, xlab="", ylab="Count")
  }else{
    plot(c(1, tseq*td), Zm, col="black", ylim=c(1, max(ZR)), typ="l", main=description, xlab="", ylab="Count")
  }
  polygon(xseq, yseq,
          col = shadecol,
          border = shadecol)
  lines(c(1, tseq), Zm, col=linecol, lwd=2)
  return(0)
}




forecast_band2 <- function(ZR1, linecol1, shadecol1, ZR2, linecol2, shadecol2, description){
  Zm1 = colMeans(ZR1)
  Zm2 = colMeans(ZR2)
  tseq = 2:ncol(ZR1)
  xseq = c(1, tseq, rev(tseq), 1)
  
  yseq1 = rep(0, ncol(ZR1))
  yseq2 = rep(0, ncol(ZR1))
  for (t in 1:length(yseq1)){
    yseq1[t] = quantile(ZR1[,t], 0.10)
    yseq2[t] = quantile(ZR1[,t], 0.90)
  }
  yseq = c(yseq1, rev(yseq2))
  
  if (logscale){
    plot(c(1, tseq*td), Zm1, col="black", ylim=c(1, max(ZR1, ZR2)), log='y', typ="l", main=description, xlab="", ylab="Count")
  }else{
    plot(c(1, tseq*td), Zm1, col="black", ylim=c(1, max(ZR1, ZR2)), typ="l", main=description, xlab="", ylab="Count", cex.main=1.5)
  }
  polygon(xseq, yseq,
          col = shadecol1,
          border = shadecol1)
  lines(c(1, tseq*td), Zm1, col=linecol1, lwd=2)
  
  yseq1 = rep(0, ncol(ZR2))
  yseq2 = rep(0, ncol(ZR2))
  for (t in 1:length(yseq1)){
    yseq1[t] = quantile(ZR2[,t], 0.10)
    yseq2[t] = quantile(ZR2[,t], 0.90)
  }
  yseq = c(yseq1, rev(yseq2))
  polygon(xseq, yseq,
          col = shadecol2,
          border = shadecol2)
  
  lines(c(1, tseq*td), Zm2, col=linecol2, lwd=2)
  return(0)
}
