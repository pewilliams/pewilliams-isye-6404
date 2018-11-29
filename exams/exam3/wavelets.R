library(data.table)
library(waveslim)


aedat <- fread('/Users/peterwilliams/Projects/pewilliams-isye-6404/exams/exam3/data/AEP_hourly.csv', 
               header = T, stringsAsFactors = F)
aep <- aedat[(nrow(aedat) - 1024 + 1):nrow(aedat)]
aep <- transform(aep, dt = as.POSIXct(Datetime, format = '%m/%d/%y %H:%M', tz = 'EST5EDT') )
x <- as.numeric(aep$AEP_MW)
#mra.out <- mra(X =  x, n.levels = 4, filter = 'haar', method = 'dwt')

#MRA plots
n.levels <- 2
aep.haar <- mra(x = x, "haar", n.levels, "dwt")
#names(aep.haar) <- c("d1", "d2", "d3", "d4", "s4")
nplots <- n.levels + 2
par(mfcol=c(nplots,2), pty="m", mar=c(5-2.5,4,4-2,1.5))
plot.ts(x, axes=F, main="Haar Filter: MRA Resolution Plot", ylab = 'original', bty='n')
axis(side=1, at=seq(0,1024,by=128), 
     labels=seq(0,1024,by=128))
axis(side=2, at=c(round(min(x),digits = -4),round(max(x),digits = -4)), 
     labels=c(round(min(x),digits = -4),round(max(x),digits = -4)))
for(i in 1:(n.levels +1)){
  plot.ts(aep.haar[[i]], axes=F, ylab=names(aep.haar)[i],bty='n')
  axis(side=1, at=seq(0,1024,by=128), 
     labels=seq(0,1024,by=128))
  axis(side=2, at=c(round(min(aep.haar[[i]]),digits = -2),round(max(aep.haar[[i]]),digits = -2)), 
       labels=c(round(min(aep.haar[[i]]),digits = -2),round(max(aep.haar[[i]]),digits = -2)))
}



n.levels <- 2
aep.haar <- mra(x = x, "la8", n.levels, "dwt")
#names(aep.haar) <- c("d1", "d2", "d3", "d4", "s4")
nplots <- n.levels + 2
#par(mfcol=c(nplots,1), pty="m", mar=c(5-3.5,4,4-2,1.5))
plot.ts(x, axes=F, main="Daubechies Filter (L=8): MRA Resolution Plot", ylab = 'Original', bty='n')
axis(side=1, at=seq(0,1024,by=128), 
     labels=seq(0,1024,by=128))
axis(side=2, at=c(round(min(x),digits = -4),round(max(x),digits = -4)), 
     labels=c(round(min(x),digits = -4),round(max(x),digits = -4)))
for(i in 1:(n.levels +1)){
  plot.ts(aep.haar[[i]], axes=F, ylab=names(aep.haar)[i],bty='n')
  axis(side=1, at=seq(0,1024,by=128), 
       labels=seq(0,1024,by=128))
  axis(side=2, at=c(round(min(aep.haar[[i]]),digits = -2),round(max(aep.haar[[i]]),digits = -2)), 
       labels=c(round(min(aep.haar[[i]]),digits = -2),round(max(aep.haar[[i]]),digits = -2)))
}










da.thresh(wc, alpha = .05, max.level = 4, verbose = FALSE, return.thresh = FALSE)
hybrid.thresh(wc, max.level = 4, verbose = FALSE, seed = 0)
manual.thresh(wc, max.level = 4, value, hard = TRUE)
sure.thresh(wc, max.level = 4, hard = TRUE)
universal.thresh(wc, max.level = 4, hard = TRUE)
universal.thresh.modwt(wc, max.level = 4, hard = TRUE)



