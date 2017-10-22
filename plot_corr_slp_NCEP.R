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

plot_corrs_slp<-function(year1,year2,dir="gcyc_out",name="proj100_highs_rad10cv0.075",cv=NA,dur=NA,slpthresh=NA,month1=1,month2=12,fout=NA,hilo="low",members=56,
                          iname="SOI",ifile=NaN)
{
years=seq(year1,year2,1)

## First, get the slp and make an slp grid

a=nc_open("/g/data/eg3/asp561/NCEP1/slp.mon.mean.nc")
lat=ncvar_get(a,"lat")
lon=ncvar_get(a,"lon")
slp=ncvar_get(a,"slp")

dd=dim(slp)
if(month2>=month1) dd2=floor(dd[3]/12)-1 else dd2=floor(dd[3]/12)-2
if(month2>=month1) month2a=month2 else month2a=month2+12
years2=seq(1948,by=1,length.out=dd2+1)
  
slp2=array(NaN,c(dd[1],dd[2],dd2+1))

for(n in 0:dd2) slp2[,,n+1]=apply(slp[,,(n*12+month1):(n*12+month2a)],c(1,2),mean)

if(lat[2]<lat[1])
{
  lat=lat[length(lat):1]
  slp2=slp2[,length(lat):1,]
}

lat2=seq(-89.5,89.5)
lon2=seq(0.5,359.5)

print('SLP done')

### Next, get the cyclones for that grid point +- 5

if(month2>=month1) systems<-array(0,c(length(lon),length(lat),length(years),2)) else systems<-array(0,c(length(lon),length(lat),length(years)-1,2))
dimnames(systems)[[4]]=c("Freq","MSLP")

for(y in 1:length(years))
{
print(years[y])
  fname=paste(dir,"/NCEP1/",name,"/tracks_",years[y],".dat",sep="")
read.table(fname, sep="",skip=1)->fixes
colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Depth","Radius","Up","Vp")
#fixes=fixes[fixes$Lat<0,]
fixes$Year=floor(fixes$Date/10000)
fixes$CV=abs(fixes$CV)
fixes$Depth=abs(fixes$Depth)
if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
fixes$Month=floor(fixes$Date/100)%%100
fixes$Lat2=floor(fixes$Lat) ## Let's skip the first step, and just make it the 10 degrees
fixes$Lon2=floor((fixes$Lon%%360))
print(paste(dur,cv,slpthresh))
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
if(!is.na(cv)) fixes=fixes[fixes$CV>=cv,]
if(!is.na(slpthresh)) fixes=fixes[fixes$MSLP>=slpthresh,]

if(month2>=month1)
 {
  I=which(fixes$Month>=month1 & fixes$Month<=month2)
  tmp=table(factor(fixes$Lon2[I],levels=seq(0,359)),factor(fixes$Lat2[I],levels=seq(-90,89)))
  tmp2=aggregate(fixes$MSLP[I],by=list(fixes$Lon2[I],fixes$Lat2[I]),FUN=mean)
  tmp3=array(NaN,dim(tmp))
  jlat=match(tmp2[,2],seq(-90,89))
  ilon=match(tmp2[,1],seq(0,359))
  for(nn in 1:length(tmp2[,1])) tmp3[ilon[nn],jlat[nn]]=as.numeric(tmp2[nn,3])
  
    for(i in 1:length(lon))
    for(j in 1:length(lat))
    {
      I=which((lon2>=lon[i]-5 & lon2<lon[i]+5) | lon2>=lon[i]+355 | lon2<lon[i]-355)
      J=which(lat2>=lat[j]-5 & lat2<lat[j]+5)
      
      systems[i,j,y,1]=sum(tmp[I,J],na.rm=T)
      systems[i,j,y,2]=sum(tmp[I,J]*tmp3[I,J],na.rm=T)/sum(tmp[I,J],na.rm=T)
    }
 } else {
  if(y==1)
  {
   I=which(fixes$Month>=month1)
   store=fixes[I,]
  } else {
   I=which(fixes$Month<=month2)
   fixes2=rbind(store,fixes[I,])

   tmp=table(factor(fixes2$Lon2[I],levels=seq(0,359)),factor(fixes2$Lat2[I],levels=seq(-90,89)))
   tmp2=aggregate(fixes2$MSLP[I],by=list(fixes2$Lon2[I],fixes2$Lat2[I]),FUN=mean)
   tmp3=array(NaN,dim(tmp))
   jlat=match(tmp2[,2],seq(-90,89))
   ilon=match(tmp2[,1],seq(0,359))
   for(nn in 1:length(tmp2[,1])) tmp3[ilon[nn],jlat[nn]]=as.numeric(tmp2[nn,3])
   
   for(i in 1:length(lon))
     for(j in 1:length(lat))
     {
       I=which((lon2>=lon[i]-5 & lon2<lon[i]+5) | lon2>=lon[i]+355 | lon2<lon[i]-355)
       J=which(lat2>=lat[j]-5 & lat2<lat[j]+5)
       
       systems[i,j,y-1,1]=sum(tmp[I,J],na.rm=T)
       systems[i,j,y-1,2]=sum(tmp[I,J]*tmp3[I,J],na.rm=T)/sum(tmp[I,J],na.rm=T)
     }

   I=which(fixes$Month>=month1)
   store=fixes[I,]
  }
}
} # End year loop

print('Cyclones done')

if(month2<month1) years=seq(year1,year2-1,1)
slp2=slp2[,,years2%in%years]

### Correlations
meanfreq=apply(systems,c(1,2,4),mean)
cyccorr<-array(NaN,c(length(lon),length(lat),3,4))

for(i in 1:length(lon))
  for(j in 1:length(lat))
    for(k in 1:2)
     if(meanfreq[i,j,1]>=1 & !is.na(meanfreq[i,j,k]))
     {
      a=cor.test(systems[i,j,,k],slp2[i,j,])
      cyccorr[i,j,k,1]=a$estimate
      cyccorr[i,j,k,2]=a$p.value      
      a=cor.test(diff(systems[i,j,,k]),diff(slp2[i,j,]))
      cyccorr[i,j,k,3]=a$estimate
      cyccorr[i,j,k,4]=a$p.value      
     }

print('Corrs done')

### Plot

type=c("freq","MSLP","CV")
breaks=seq(-1,1,0.1)
cols=col_anom(length(breaks)-1)

print("Plotting")

if(is.na(fout)) fname=paste("NCEP1_slpcorr_",year1,"_",year2,"_global.pdf",sep="") else 
  fname=paste(fout,"_slpcorr_global.pdf",sep="")

pdf(file=fname,width=14,height=10)
layout(cbind(c(1,2,5),c(3,4,5)),height=c(1,1,0.15))
par(mar=c(2,2,4,1))

for(k in 1:2)
{
image(lon,lat,cyccorr[,,k,1],breaks=breaks,col=cols,xlab="",ylab="",
          main=paste(type[k],"Corr:",year1,"-",year2))
map('world2',add=T)
sigmask=which(cyccorr[,,k,2]<0.05 & meanfreq[,,k]>=5,arr.ind=T)
points(lon[sigmask[,1]],lat[sigmask[,2]],col="black",pch=19,cex=0.2)

image(lon,lat,cyccorr[,,k,3],breaks=breaks,col=cols,xlab="",ylab="",
      main=paste(type[k],"DifferenceCorr:",year1,"-",year2))
map('world2',add=T)
sigmask=which(cyccorr[,,k,4]<0.05 & meanfreq[,,k]>=5,arr.ind=T)
points(lon[sigmask[,1]],lat[sigmask[,2]],col="black",pch=19,cex=0.2)
}

ColorBar(breaks,cols,subsampleg=1,vert=F)
dev.off()
}

seas=c("","_MAM","_JJA","_SON","_DJF")
m1=c(1,3,6,9,12)
m2=c(12,5,8,11,2)

for(s in 1:5)
{
  plot_corrs_slp(1950,2016,dir="/short/eg3/asp561/cts.dir/gcyc_out/",name="proj100_highs_rad10cv0.075_v2",
                 fout=paste("UM_NCEP1_anticyclonecorr_1950-2016_rad10cv075_D2",seas[s],sep=""),
                 hilo="high",month1=m1[s],month2=m2[s],dur=2,iname=ind)

  plot_corrs_slp(1950,2016,dir="/short/eg3/asp561/cts.dir/gcyc_out/",name="proj100_lows_rad5cv0.15_v2",
                 fout=paste("UM_NCEP1_cyclonecorr_1950-2016_rad5cv15_D2",seas[s],sep=""),
                 hilo="low",month1=m1[s],month2=m2[s],dur=2,iname=ind)



}
