## This is code that takes the 11 ensemble members of a POAMA run
## Started on a particular day
## And outputs some interested results for further analysis

## I care about ECLs, so main result is a daily timeseries of # of members that
## have an ECL on that day

## Also outputs a figure for each of the months in the file
## That shows the number of cyclones centred in that month
## Option to set a higher intensity threshold

library(maps)
library(ncdf4)

color.palette <- function(steps, n.steps.between=NULL, ...){
  
  if(is.null(n.steps.between)) n.steps.between <- rep(0, (length(steps)-1))
  if(length(n.steps.between) != length(steps)-1) stop("Must have one less n.steps.between value than steps")
  
  fill.steps <- cumsum(rep(1, length(steps))+c(0,n.steps.between))
  RGB <- matrix(NA, nrow=3, ncol=fill.steps[length(fill.steps)])
  RGB[,fill.steps] <- col2rgb(steps)
  
  for(i in which(n.steps.between>0)){
    col.start=RGB[,fill.steps[i]]
    col.end=RGB[,fill.steps[i+1]]
    for(j in seq(3)){
      vals <- seq(col.start[j], col.end[j], length.out=n.steps.between[i]+2)[2:(2+n.steps.between[i]-1)]  
      RGB[j,(fill.steps[i]+1):(fill.steps[i+1]-1)] <- vals
    }
  }
  
  new.steps <- rgb(RGB[1,], RGB[2,], RGB[3,], maxColorValue = 255)
  pal <- colorRampPalette(new.steps, ...)
  return(pal)
}

col_anom <- color.palette(c("darkblue","blue","white","red","darkred"),c(10,20,20,10))
col_val <- color.palette(c("white","blue","darkblue","black"),c(20,10,5))

ColorBar <- function(brks,cols,vert=T,subsampleg=1)
{
  if(vert) {
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(1.5, length(brks) - 1.5, subsampleg), tick = TRUE, 
       labels = brks[seq(2, length(brks)-1, subsampleg)])
  } else {
    par(mar = c(1.5, 1, 1, 1), mgp = c(1.5, 0.3, 0), las = 1, cex = 1)
    image(1:length(cols), 1, t(t(1:length(cols))), axes = FALSE, col = cols,
          xlab = '', ylab = '')
    box()
    axis(1, at = seq(1.5, length(brks) - 1.5, subsampleg), 
         labels = brks[seq(2, length(brks)-1, subsampleg)])
  }
}

plot_corrs_slp<-function(year1,year2,dir="gcyc_out",anti="proj100_highs_rad10cv0.075",cyc="proj100_lows_rad5cv0.15",hcv=NA,lcv=NA,dur=NA,hslp=NA,lslp=NA,month1=1,month2=12,fout=NA)
{
years=seq(year1,year2,1)

### Next, get the cyclones for that grid point +- 5
lat=seq(-85,85,10)
lon=seq(5,355,10)  ### Can always combine into bigger cells later
if(month2>=month1) systems<-array(0,c(length(lon),length(lat),length(years),2,4)) else systems<-array(0,c(length(lon),length(lat),length(years)-1,2,4))
dimnames(systems)[[5]]=c("Freq","MSLP","CV","Depth")
dimnames(systems)[[4]]=c("Cyclones","Anticyclones")

name=c(cyc,anti)
cv=c(lcv,hcv)

for(n in 1:2)
for(y in 1:length(years))
{
print(years[y])
fname=paste(dir,name[n],"/tracks_",years[y],".dat",sep="")
read.table(fname, sep="",skip=1)->fixes
colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Depth","Radius","Up","Vp")
#fixes=fixes[fixes$Lat<0,]
fixes$Year=floor(fixes$Date/10000)
fixes$CV=abs(fixes$CV)
fixes$Depth=abs(fixes$Depth)
if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
fixes$Month=floor(fixes$Date/100)%%100
fixes$Lat2=10*floor(fixes$Lat/10) ## Let's skip the first step, and just make it the 10 degrees
fixes$Lon2=10*floor((fixes$Lon%%360)/10)

### Make table of events to combine with DJF for exclusion
 if(!is.na(dur))
  {
  x<-rle(fixes$ID)
  events<-cbind(x$values,x$lengths,matrix(data=0,nrow=length(x$values),ncol=1))
  events=events[events[,2]>=dur,]
  include<-match(fixes[,1],events[,1])
  J<-which(is.na(include)==0)
  fixes=fixes[J,]
  }
if(!is.na(cv[n])) fixes=fixes[fixes$CV>=cv[n],]

if(n==1 & !is.na(lslp)) fixes=fixes[fixes$MSLP<=lslp,]
if(n==2 & !is.na(hslp)) fixes=fixes[fixes$MSLP>=hslp,]

if(month2>=month1)
 {
  I=which(fixes$Month>=month1 & fixes$Month<=month2)
  systems[,,y,n,1]=table(factor(fixes$Lon2[I],levels=seq(0,350,10)),factor(fixes$Lat2[I],levels=seq(-90,80,10)))
  tmp2=aggregate(fixes[I,8:10],by=list(fixes$Lon2[I],fixes$Lat2[I]),FUN=mean)
  jlat=match(tmp2[,2],seq(-90,80,10))
  ilon=match(tmp2[,1],seq(0,350,10))
  for(nn in which(!is.na(jlat))) systems[ilon[nn],jlat[nn],y,n,2:4]=as.numeric(tmp2[nn,3:5])
 } else {
  if(y==1)
  {
   I=which(fixes$Month>=month1)
   store=fixes[I,]
  } else {
   I=which(fixes$Month<=month2)
   fixes2=rbind(store,fixes[I,])

   systems[,,y-1,n,1]=table(factor(fixes2$Lon2,levels=seq(0,350,10)),factor(fixes2$Lat2,levels=seq(-90,80,10)))
   tmp2=aggregate(fixes2[,8:10],by=list(fixes2$Lon2,fixes2$Lat2),FUN=mean)
   jlat=match(tmp2[,2],seq(-90,80,10))
   ilon=match(tmp2[,1],seq(0,350,10))
   for(nn in which(!is.na(jlat))) systems[ilon[nn],jlat[nn],y-1,n,2:4]=as.numeric(tmp2[nn,3:5])

   I=which(fixes$Month>=month1)
   store=fixes[I,]
  }
}
} # End year loop

print('Cyclones done')

if(month2<month1) years=seq(year1,year2-1,1)

### Correlations
meanfreq=apply(systems,c(1,2,4,5),mean)
cyccorr<-array(NaN,c(length(lon),length(lat),4,4))

for(i in 1:length(lon))
  for(j in 1:length(lat))
    for(k in 1:4)
     if(meanfreq[i,j,1,1]>=1 & meanfreq[i,j,2,1]>=1)
     {
      a=cor.test(systems[i,j,,1,k],systems[i,j,,2,k])
      cyccorr[i,j,1,k]=a$estimate
      cyccorr[i,j,2,k]=a$p.value      
      a=cor.test(diff(systems[i,j,,1,k]),diff(systems[i,j,,2,k]))
      cyccorr[i,j,3,k]=a$estimate
      cyccorr[i,j,4,k]=a$p.value      
     }

print('Corrs done')

### Plot
breaks=seq(-1,1,0.1)
cols=col_anom(length(breaks)-1)
type=c("freq","MSLP","CV","Depth")
print("Plotting")

for(k in 1:4)
{
if(is.na(fout)) fname=paste("NCEP1_cycanticorr_",year1,"_",year2,"_",type[k],"_global.pdf",sep="") else 
  fname=paste(fout,"_",type[k],"_global.pdf",sep="")

pdf(file=fname,width=7,height=10)
layout(c(1,2,3),height=c(1,1,0.15))
par(mar=c(2,2,4,1))

image(lon,lat,cyccorr[,,1,k],breaks=breaks,col=cols,xlab="",ylab="",
          main=paste("Corr:",year1,"-",year2))
map('world2',add=T)
sigmask=which(cyccorr[,,2,k]<0.05 & meanfreq[,,1,k]>=1 & meanfreq[,,2,k]>=1,arr.ind=T)
points(lon[sigmask[,1]],lat[sigmask[,2]],col="black",pch=19,cex=0.6)

image(lon,lat,cyccorr[,,3,k],breaks=breaks,col=cols,xlab="",ylab="",
      main=paste("DifferenceCorr:",year1,"-",year2))
map('world2',add=T)
sigmask=which(cyccorr[,,4,k]<0.05 & meanfreq[,,1,k]>=1 & meanfreq[,,2,k]>=1,arr.ind=T)
points(lon[sigmask[,1]],lat[sigmask[,2]],col="black",pch=19,cex=0.6)

ColorBar(breaks,cols,subsampleg=1,vert=F)
dev.off()
}
}

seas=c("","_MAM","_JJA","_SON","_DJF")
m1=c(1,3,6,9,12)
m2=c(12,5,8,11,2)

for(s in 1:5)
{
  plot_corrs_slp(1950,2016,dir="/short/eg3/asp561/cts.dir/gcyc_out/NCEP1/",anti="proj100_highs_rad10cv0.075_v2",cyc="proj100_lows_rad5cv0.15_v2",
                 fout=paste("UM_NCEP1_hilocorr_1950-2016_D2",seas[s],sep=""),
                 month1=m1[s],month2=m2[s],dur=2,hslp=NaN,lslp=NaN)
}
