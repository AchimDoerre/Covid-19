get_daynr <- function(theday, tage){
  r = 0
  m = length(tage)
  j=0
  while ((j<m) && (r == 0)){
    j=j+1
    if (tage[j] == theday){
      r = j
    }
  }
  return(r)
}




get_ticks_and_labels <- function(day){
  floorday = 100*floor(as.numeric(day)/100)+1
  r0 = get_daynr(floorday, tage)
  rt = get_daynr(day, tage)
  theticks = cumsum(c(1, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,
                      31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30))
  themonths = rep(c("Jan", "Feb", "Mrz", "Apr", "Mai", "Jun", "Jul", "Aug", "Sep", "Okt", "Nov", "Dez"), 2)
  
  j = 1
  while (theticks[j] < r0){
    j=j+1
  }
  ticks = (theticks[j:length(theticks)]-theticks[4]-(rt-r0))/td
  monlabels = themonths[j:length(themonths)]
  result = list(ticks = ticks, monlabels = monlabels)
  return(result)
}


convert_day <- function(day){
  YY = floor(day/10000)
  YYYY = 2000 + YY
  MM = floor((day - 10000*YY)/100)
  DD = day - 10000*YY - 100*MM
  r = paste(DD, ".", MM, ".", YYYY, sep="")
  return(r)
}


convert_days <- function(days){
  r = c()
  for (j in 1:length(days)){
    r = c(r, convert_day(days[j]))
  }
  return(r)
}




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
    plot(c(1, tseq), Zm, col="black", ylim=c(1, max(ZR)), log='y', xaxt='n', typ="l", main=description, xlab="", ylab="Count")
  }else{
    plot(c(1, tseq), Zm, col="black", ylim=c(1, max(ZR)), xaxt='n', typ="l", main=description, xlab="", ylab="Count")
  }
  polygon(xseq, yseq,
          col = shadecol,
          border = shadecol)
  lines(c(1, tseq), Zm, col=linecol, lwd=2)
  axis(1, at=ticks-tick0, labels=monlabels)
  return(0)
}




forecast_band2 <- function(ZR1, linecol1, shadecol1, ZR2, linecol2, shadecol2, description="", usemonticks=TRUE){
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
    if (usemonticks){
      plot(c(1, tseq), Zm1, col="black", ylim=c(1, max(ZR1, ZR2)), log='y', xaxt='n', typ="l", main=description, xlab="", ylab="Count")
    }else{
      plot(c(1, tseq), Zm1, col="black", ylim=c(1, max(ZR1, ZR2)), log='y', typ="l", main=description, xlab="", ylab="Count")
    }
  }else{
    if (usemonticks){
      plot(c(1, tseq), Zm1, col="black", ylim=c(1, max(ZR1, ZR2)), xaxt='n', typ="l", main=description, xlab="", ylab="Count", cex.main=1.5)
    }else{
      plot(c(1, tseq), Zm1, col="black", ylim=c(1, max(ZR1, ZR2)), typ="l", main=description, xlab="", ylab="Count", cex.main=1.5)
    }
  }
  polygon(xseq, yseq,
          col = shadecol1,
          border = shadecol1)
  lines(c(1, tseq), Zm1, col=linecol1, lwd=2)
  
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
  
  lines(c(1, tseq), Zm2, col=linecol2, lwd=2)
  if (usemonticks){
    axis(1, at=ticks-tick0, labels=monlabels, cex.axis=1.5)
  }
  return(0)
}
