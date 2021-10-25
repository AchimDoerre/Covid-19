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
    plot(c(1, tseq*td), rep(Inf, length(tseq)+1), col="black", ylim=c(1, max(yseq)), log='y', typ="l", main=description, xlab="Days", ylab="Count")
  }else{
    plot(c(1, tseq*td), rep(-10000, length(tseq)+1), col="black", ylim=c(1, max(yseq)), typ="l", main=description, xlab="Days", ylab="Count")
  }
  polygon(xseq*td, yseq,
          col = shadecol,
          border = shadecol)
  lines(c(1, tseq*td), Zm, col=linecol, lwd=2)
  
  return(0)
}
