##Yannik Friedli & Thierry Müller Climatology III --> Cosigüina eruption 1835
##16.12.2020

library(ncdf4)
library(maps)

mycol <- c("#4B7DB8","#70A9CF","#99CAE1","#C3E5F0","#E7F6EB","#FFFEBE","#FFE89C","#FDC576","#FA9857","#F26841")
mycol2 <- c("#8e0152", "#c51b7d", "#de77ae", "#f1b6da", "#fde0ef", "#e6f5d0", "#b8e186", "#7fbc41", "#4d9221", "#276419")
mycol3 <- c("#543005", "#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#c7eae5", "#80cdc1", "#35978f", "#01665e", "#003c30")

#reading in EKF400_v2 data
f1 <- nc_open("EKF400v2.0_ensmean.nc")
lon <- f1$dim[[1]]$vals
lat <- f1$dim[[2]]$vals
dlon <- lon[2]-lon[1]
dlat <- lat[1]-lat[2]

## coordinates Cosigüina volcano
clon <- 267.567
clat <- 12.983
sellon <- which((lon<(clon+(dlon/2)))&(lon>(clon-(dlon/2))))
sellat <- which((lat<(clat+(dlat/2)))&(lat>(clat-(dlat/2))))

## coordinates Lewes avalanche
Llon <- 0.022
Llat <- 51.872
selLlon <- which((lon<(Llon+(dlon/2)))&(lon>(Llon-(dlon/2))))
selLlat <- which((lat<(Llat+(dlat/2)))&(lat>(Llat-(dlat/2))))

#Europe
ElonW <- 0
ElonE <- 30
ElatN <- 55
ElatS <- 40

selElon <- which((lon<(ElonE+(dlon/2)))&(lon>(ElonW-(dlon/2))))
selElat <- which((lat<(ElatN+(dlat/2)))&(lat>(ElatS-(dlat/2))))

#aggregate monthly data to seasonal composites (sommer and winter)
yr <- c(1602:2002) 

temp.sum.mean <- array(NA,dim=c(192,96,401)) #array filled with NA; dimensions of 192 x 96 x 401 (number of years)
temp.win.mean <- array(NA,dim=c(192,96,401))
#temp.yr.mean <- array(NA, dim=c(192,96,401))
prec.sum.mean <- array(NA, dim=c(192,96,401))
prec.win.mean <- array(NA, dim=c(192,96,401))
slp.sum.mean <- array(NA, dim=c(192,96,401))
slp.win.mean <- array(NA, dim=c(192,96,401))
hgt.sum.mean <- array(NA, dim=c(192,96,401))
hgt.win.mean <- array(NA, dim=c(192,96,401))

# years following volcanic eruptions
volc <- c(1620, 1642, 1674, 1694, 1712, 1720, 1730, 1739, 1756, 1762, 1784, 1795, 1797, 1810, 1816, 1817, 1823, 1832, 1836, 1862, 1884, 1887, 1903, 1913, 1926, 1944, 1964, 1977, 1983, 1992)
selvolc <- match(volc,yr) 

#________________________________________________________________________________________________________
#calculating mean values

#temperature winter & summer months
for (i in 1602:2002){
  temp.win.mean[,,i-1601] <- apply(ncvar_get(f1,varid="t2m",start=c(1,1,((i-1602)*12)+12),count=c(-1,-1,3)),c(1,2),mean) #for each longitude and latitude the mean temperature over 3 months
  temp.sum.mean[,,i-1601] <- apply(ncvar_get(f1,varid="t2m",start=c(1,1,((i-1602)*12)+18),count=c(-1,-1,3)),c(1,2),mean)
} 

#whole year
#for (i in 1602:2003){
 # temp.yr.mean[,,i-1601] <- apply(ncvar_get(f1,varid="t2m", #first time stamp is done
  #                                            start=c(1,1,((i-1602)*12)+1), #if we want january as start of year it's +1  /1602 has to stay the same = beginning of the dataset # last part says when we are reading, at what timestamp
   #                                           count=c(-1,-1,12)),c(1,2),mean)
#}

#precipitation winter & summer months
for (i in 1602:2002){
  prec.win.mean[,,i-1601] <- apply(ncvar_get(f1,varid="prec",start=c(1,1,((i-1602)*12)+12),count=c(-1,-1,3)),c(1,2),mean)
  prec.sum.mean[,,i-1601] <- apply(ncvar_get(f1,varid="prec",start=c(1,1,((i-1602)*12)+18),count=c(-1,-1,3)),c(1,2),mean)
}

#sea level pressure winter & summer months
for (i in 1602:2002){
  slp.win.mean[,,i-1601] <- apply(ncvar_get(f1,varid="slp",start=c(1,1,((i-1602)*12)+12),count=c(-1,-1,3)),c(1,2),mean)
  slp.sum.mean[,,i-1601] <- apply(ncvar_get(f1,varid="slp",start=c(1,1,((i-1602)*12)+18),count=c(-1,-1,3)),c(1,2),mean)
}

#geopotential height winter & summer months
for (i in 1602:2002){
  hgt.win.mean[,,i-1601] <- apply(ncvar_get(f1,varid="hgt500",start=c(1,1,((i-1602)*12)+12),count=c(-1,-1,3)),c(1,2),mean)
  hgt.sum.mean[,,i-1601] <- apply(ncvar_get(f1,varid="hgt500",start=c(1,1,((i-1602)*12)+18),count=c(-1,-1,3)),c(1,2),mean)
}


#________________________________________________________________________________________________________

#Definiton of years we want to analyze
cosiW1 <- 1835 + 1 #the first winter after the eruption (which is not the winter in which the eruption took place)
cosiW2 <- 1835 + 2 #the second winter after the erupion
cosiComS <- c(1835, 1835 + 1) #the summer following the eruption and the summer a year after the eruption
#selyr <- cosi-1602 
selcosiW1 <- match(cosiW1,yr)
selcosiW2 <- match(cosiW2,yr)
selcosiComS <- match(cosiComS, yr)
#________________________________________________________________________________________________________
## Case Study Event

CosiguinaW1 <- apply(temp.win.mean[,,selcosiW1],c(1,2), mean)
CosiguinaW2 <- apply(temp.win.mean[,,selcosiW2],c(1,2), mean)
CosiguinaS <- apply(temp.sum.mean[,,selcosiComS],c(1,2), mean)
#CosiguinaYr <- apply(temp.yr.mean[,,selcosiCom],c(1,2), mean)

CosiguinaW1Prec <- apply(prec.win.mean[,,selcosiW1],c(1,2), mean)
CosiguinaW2Prec <- apply(prec.win.mean[,,selcosiW2],c(1,2), mean)
CosiguinaSPrec <- apply(prec.sum.mean[,,selcosiComS],c(1,2), mean)

CosiguinaW1SLP <- apply(slp.win.mean[,,selcosiW1],c(1,2), mean)
CosiguinaW2SLP <- apply(slp.win.mean[,,selcosiW2],c(1,2), mean)
CosiguinaSSLP <- apply(slp.sum.mean[,,selcosiComS],c(1,2), mean)


CosiguinaW1hgt <- apply(hgt.win.mean[,,selcosiW1],c(1,2), mean)
CosiguinaW2hgt <- apply(hgt.win.mean[,,selcosiW2],c(1,2), mean)
CosiguinaShgt <- apply(hgt.sum.mean[,,selcosiComS],c(1,2), mean)
#________________________________________________________________________________________________________
##Reference Period

select <- c(154:234) #year number 1835 = 234; here all years from 1602 until 1900 (c(1:298)); 1750-1800 (c(149:199)); 1755-1835 (c(154:234))
select <- setdiff(select,selvolc)


refWin <- apply(temp.win.mean[,,select],c(1,2),mean) 
refSum <- apply(temp.sum.mean[,,select],c(1,2),mean)
#refyr <- apply(temp.yr.mean, c(1,2), mean)

refWinPrec <- apply(prec.win.mean[,,select],c(1,2),mean) 
refSumPrec <- apply(prec.sum.mean[,,select],c(1,2),mean)
#refyrPrec <- apply(prec.yr.mean, c(1,2), mean)

refWinSLP <- apply(slp.win.mean[,,select],c(1,2),mean) 
refSumSLP <- apply(slp.sum.mean[,,select],c(1,2),mean)
#refyrPrec <- apply(prec.yr.mean, c(1,2), mean)

refWinhgt <- apply(hgt.win.mean[,,select],c(1,2),mean) 
refSumhgt <- apply(hgt.sum.mean[,,select],c(1,2),mean)
#refyrPrec <- apply(prec.yr.mean, c(1,2), mean)
#________________________________________________________________________________________________________

diffWin1 <- CosiguinaW1-refWin
diffWin2 <- CosiguinaW2-refWin
diffSum <- CosiguinaS-refSum
#diffyr <- CosiguinaYr-refyr

diffWin1Prec <- CosiguinaW1Prec-refWinPrec
diffWin2Prec <- CosiguinaW2Prec-refWinPrec
diffSumPrec <- CosiguinaSPrec-refSumPrec
#diffyr <- CosiguinaYr-refyr

diffWin1SLP <- CosiguinaW1SLP-refWinSLP
diffWin2SLP <- CosiguinaW2SLP-refWinSLP
diffSumSLP <- CosiguinaSSLP-refSumSLP
#diffyr <- CosiguinaYr-refyr

diffWin1hgt <- CosiguinaW1hgt-refWinhgt
diffWin2hgt <- CosiguinaW2hgt-refWinhgt
diffSumhgt <- CosiguinaShgt-refSumhgt
#diffyr <- CosiguinaYr-refyr
#________________________________________________________________________________________________________

### then prepare plotting
### change lon is from 0:360
x <- c(lon[lon>180]-360,lon[lon<=180])
y <- rev(lat)
mycol0 <- c("blue4","blue2","cornflowerblue","cadetblue2","azure2","bisque1","burlywood1","brown1","brown3","darkred")
min.lon <- -38 #ratio 5:3 looks good
max.lon <- 45
min.lat <- 28
max.lat <- 76

zW1 <- rbind(diffWin1[lon>180,rev(1:length(lat))],diffWin1[lon<=180,rev(1:length(lat))])
zW2 <- rbind(diffWin2[lon>180,rev(1:length(lat))],diffWin2[lon<=180,rev(1:length(lat))])
zS <- rbind(diffSum[lon>180,rev(1:length(lat))],diffSum[lon<=180,rev(1:length(lat))])
#zY <- rbind(diffyr[lon>180,rev(1:length(lat))],diffyr[lon<=180,rev(1:length(lat))])

zW1P <- rbind(diffWin1Prec[lon>180,rev(1:length(lat))],diffWin1Prec[lon<=180,rev(1:length(lat))])
zW2P <- rbind(diffWin2Prec[lon>180,rev(1:length(lat))],diffWin2Prec[lon<=180,rev(1:length(lat))])
zSP <- rbind(diffSumPrec[lon>180,rev(1:length(lat))],diffSumPrec[lon<=180,rev(1:length(lat))])
#zY <- rbind(diffyr[lon>180,rev(1:length(lat))],diffyr[lon<=180,rev(1:length(lat))])

zW1S <- rbind(diffWin1SLP[lon>180,rev(1:length(lat))],diffWin1SLP[lon<=180,rev(1:length(lat))])
zW2S <- rbind(diffWin2SLP[lon>180,rev(1:length(lat))],diffWin2SLP[lon<=180,rev(1:length(lat))])
zSS <- rbind(diffSumSLP[lon>180,rev(1:length(lat))],diffSumSLP[lon<=180,rev(1:length(lat))])

zW1h <- rbind(diffWin1hgt[lon>180,rev(1:length(lat))],diffWin1hgt[lon<=180,rev(1:length(lat))])
zW2h <- rbind(diffWin2hgt[lon>180,rev(1:length(lat))],diffWin2hgt[lon<=180,rev(1:length(lat))])
zSh <- rbind(diffSumhgt[lon>180,rev(1:length(lat))],diffSumhgt[lon<=180,rev(1:length(lat))])
#________________________________________________________________________________________________________

### plots

#----Temperature---------------------------------------------------------------------------------------------------------
minmax <- max(abs(c(min(zW1),max(zW1))))
mylevs <- (c(0:10)-5)
filled.contour(x,y,zW1,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)})
filled.contour(x[x>min.lon&x<max.lon],y[y>min.lat&y<max.lat],zW1[x>min.lon&x<max.lon,y>min.lat&y<max.lat],
               levels=mylevs,col=mycol,
               plot.title = {par(cex.main=1); title(main = "Temperature Anomalies Winter 1835/1836")},
               plot.axes={map("world",interior=F,add=T)},
               #key.title = {par(cex.main=0.5); title(main = "\n\n\n\n\n Temperature\n[°C]")},
               key.axes = axis(4))

minmax <- max(abs(c(min(zW2),max(zW2))))
mylevs <-  0.2*minmax*(c(0:10)-5)
filled.contour(x,y,zW2,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)})
filled.contour(x[x>min.lon&x<max.lon],y[y>min.lat&y<max.lat],zW2[x>min.lon&x<max.lon,y>min.lat&y<max.lat],levels=mylevs,col=mycol,plot.axes={map("worldHires",interior=F,add=T)})
title("Temperature Anomalies Winter 1836/1837")

minmax <- max(abs(c(min(zS),max(zS))))
mylevs <-  (c(0:10)-5)
filled.contour(x,y,zS,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)})
filled.contour(x[x>min.lon&x<max.lon],y[y>min.lat&y<max.lat],zS[x>min.lon&x<max.lon,y>min.lat&y<max.lat],levels=mylevs,col=mycol,plot.axes={map("world",add=T)})
title("Temperatur Anomalies Summer 1835 & 1836")

#minmax <- max(abs(c(min(zY),max(zY))))
#mylevs <-  0.2*minmax*(c(0:10)-5)
#filled.contour(x,y,zY,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)})
#title("whole year")


#-----Precipitation--------------------------------------------------------------------------------------------------------
minmax <- max(abs(c(min(zW1P),max(zW1P))))
mylevs <-  0.2*minmax*(c(0:10)-5)
filled.contour(x,y,zW1P,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)})
filled.contour(x[x>min.lon&x<max.lon],y[y>min.lat&y<max.lat],zW1P[x>min.lon&x<max.lon,y>min.lat&y<max.lat],levels=mylevs,col=mycol2,plot.axes={map("world",interior=F,add=T)})
title("Precipitation Anomalies Winter 1835/1836")

minmax <- max(abs(c(min(zW2P),max(zW2P))))
mylevs <-  0.2*minmax*(c(0:10)-5)
filled.contour(x,y,zW2P,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)})
filled.contour(x[x>min.lon&x<max.lon],y[y>min.lat&y<max.lat],zW2P[x>min.lon&x<max.lon,y>min.lat&y<max.lat],levels=mylevs,col=mycol2,plot.axes={map("world",interior=F,add=T)})
title("Precipitation Anomalies Winter 1836/1837")

minmax <- max(abs(c(min(zSP),max(zSP))))
mylevs <-  0.2*minmax*(c(0:10)-5)
filled.contour(x,y,zSP,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)})
filled.contour(x[x>min.lon&x<max.lon],y[y>min.lat&y<max.lat],zSP[x>min.lon&x<max.lon,y>min.lat&y<max.lat],levels=mylevs,col=mycol2,plot.axes={map("world",interior=F,add=T)})
title("Precipitation Anomalies Summer 1835 & 1836")



#-----sea level pressure--------------------------------------------------------------------------------------------------------
min.lonSLP <- -179
max.lonSLP <- 179
min.latSLP <- 10
max.latSLP <- 88


minmax <- max(abs(c(min(zW1S),max(zW1S))))
mylevs <-  0.2*minmax*(c(0:10)-5)
filled.contour(x,y,zW1S,levels=mylevs,col=mycol3,plot.axes={map("world",interior=F,add=T)})
filled.contour(x[x>min.lon&x<max.lon],y[y>min.lat&y<max.lat],zW1S[x>min.lon&x<max.lon,y>min.lat&y<max.lat],levels=mylevs,col=mycol3,plot.axes={map("world",interior=F,add=T)})
filled.contour(x[x>min.lonSLP&x<max.lonSLP],y[y>min.latSLP&y<max.latSLP],zW1S[x>min.lonSLP&x<max.lonSLP, y>min.latSLP&y<max.latSLP],levels=mylevs,col=rev(mycol3),plot.axes={map("world",interior=F,add=T)})
title("Sea Level Pressure Anomalies Winter 1835/1836")

minmax <- max(abs(c(min(zW2S),max(zW2S))))
mylevs <-  0.2*minmax*(c(0:10)-5)
filled.contour(x,y,zW2S,levels=mylevs,col=mycol3,plot.axes={map("world",interior=F,add=T)})
filled.contour(x[x>min.lon&x<max.lon],y[y>min.lat&y<max.lat],zW2S[x>min.lon&x<max.lon,y>min.lat&y<max.lat],levels=mylevs,col=mycol3,plot.axes={map("world",interior=F,add=T)})
filled.contour(x[x>min.lonSLP&x<max.lonSLP],y[y>min.latSLP&y<max.latSLP],zW2S[x>min.lonSLP&x<max.lonSLP, y>min.latSLP&y<max.latSLP],levels=mylevs,col=mycol3,plot.axes={map("world",interior=F,add=T)})
title("Sea Level Pressure Anomalies Winter 1836/1837")

minmax <- max(abs(c(min(zSS),max(zSS))))
mylevs <-  0.2*minmax*(c(0:10)-5)
filled.contour(x,y,zSS,levels=mylevs,col=mycol,plot.axes={map("world",interior=F,add=T)})
filled.contour(x[x>min.lon&x<max.lon],y[y>min.lat&y<max.lat],zSS[x>min.lon&x<max.lon,y>min.lat&y<max.lat],levels=mylevs,col=mycol3,plot.axes={map("world",interior=F,add=T)})
title("Seal Level Pressure Anomalies Summer 1835 & 1836")

#------geopotential height-------------------------------------------------------------------------------------------------------

minmax <- max(abs(c(min(zW1h),max(zW1h))))
mylevs <-  0.2*minmax*(c(0:10)-5)
filled.contour(x,y,zW1h,levels=mylevs,col=mycol3,plot.axes={map("world",interior=F,add=T)})
filled.contour(x[x>min.lon&x<max.lon],y[y>min.lat&y<max.lat],zW1h[x>min.lon&x<max.lon,y>min.lat&y<max.lat],levels=mylevs,col=mycol3,plot.axes={map("world",interior=F,add=T)})
filled.contour(x[x>min.lonSLP&x<max.lonSLP],y[y>min.latSLP&y<max.latSLP],zW1h[x>min.lonSLP&x<max.lonSLP, y>min.latSLP&y<max.latSLP],levels=mylevs,col=rev(mycol3),plot.axes={map("world",interior=F,add=T)})
title("Geopotential hight at 500 Pa  Anomalies Winter 1835/1836")

minmax <- max(abs(c(min(zW2h),max(zW2h))))
mylevs <-  0.2*minmax*(c(0:10)-5)
filled.contour(x,y,zW2h,levels=mylevs,col=mycol3,plot.axes={map("world",interior=F,add=T)})
filled.contour(x[x>min.lon&x<max.lon],y[y>min.lat&y<max.lat],zW2h[x>min.lon&x<max.lon,y>min.lat&y<max.lat],levels=mylevs,col=mycol3,plot.axes={map("world",interior=F,add=T)})
filled.contour(x[x>min.lonSLP&x<max.lonSLP],y[y>min.latSLP&y<max.latSLP],zW2h[x>min.lonSLP&x<max.lonSLP, y>min.latSLP&y<max.latSLP],levels=mylevs,col=mycol3,plot.axes={map("world",interior=F,add=T)})
title("Geopotential hight at 500 Pa  Anomalies Winter 1836/1837")

minmax <- max(abs(c(min(zSh),max(zSh))))
mylevs <-  0.2*minmax*(c(0:10)-5)
filled.contour(x,y,zSh,levels=mylevs,col=mycol3,plot.axes={map("world",interior=F,add=T)})
filled.contour(x[x>min.lon&x<max.lon],y[y>min.lat&y<max.lat],zSh[x>min.lon&x<max.lon,y>min.lat&y<max.lat],levels=mylevs,col=mycol3,plot.axes={map("world",interior=F,add=T)})
title("Geopotential hight at 500 Pa Anomalies Summer 1835 & 1836")


## some random statistics

temp3 <- ncvar_get(f1,varid="t2m",start=c(selLlon,selLlat,1),count=c(1,1,-1)) #value 2796 corresponds to Januar 1835
time <- f1$dim[[3]]$vals

plot(time,filter(temp3,rep(1/12,12)),type="l")

ttE <- temp.win.mean[selElon, selElat,c(199:234)]
  hist(ttE, freq = FALSE)
  lines(density(ttE))
  ttestE36 <- t.test(ttE, mu=mean(temp.win.mean[,,235]))
  plot(density(ttE))
  CIE <- (mean(ttE) + c(-1,1)*sd(ttE)*qt(1-0.05/2, length(ttE)-1)/sqrt(length(ttE)))
  abline(v=mean(ttE), col = "blue")
  abline(v=CIE)
  abline(v=mean(temp.win.mean[,,235]))

ttt <- temp.win.mean[selLlon,selLlat,select] #all year until 1900 
  ttC <- temp.win.mean[sellon,sellat,c(1:234)]
  t36 <- temp.win.mean[selLlon,selLlat,235]#year 1 after eruption 
  t37 <- temp.win.mean[selLlon,selLlat,236] #year 2 after eruption
  hist(ttt, freq = FALSE)
  lines(density(ttt))
  ttest36 <- t.test(ttt, mu=mean(t36))
  ttest37 <- t.test(ttt, mu=t37)
  plot(density(ttt))
  abline(v=mean(ttt), col = "blue")
  CI <- (mean(ttt) + c(-1,1)*sd(ttt)*qt(1-0.05/2, length(ttt)-1)/sqrt(length(ttt)))
  abline(v=CI)
  abline(v=mean(t36))
  abline(v=t37)
  quantile(ttt, c(0.05, 0.95))

timeseries <- temp.win.mean[selLlon,selLlat,c(200:249)]
plot(timeseries, type = "b", xaxt = "none", xlab = NA, ylab = "[°C]")
abline(v=235, col="blue")
axis(1, at=1:31, labels = seq(1820, 1850))

temp3 <- ncvar_get(f1,varid="t2m",start=c(selLlon,selLlat,1),count=c(1,1,-1)) #value 2796 corresponds to Januar 1835
  ts <- temp3[2736:2856] #eruption +/- 5 years
  ts <- temp3[2796:2844]
  plot(ts, type = "l")
time <- f1$dim[[3]]$vals

plot(time,filter(temp3,rep(1/12,12)),type="l")


