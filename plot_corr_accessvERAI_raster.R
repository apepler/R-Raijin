## This is code that takes the 11 ensemble members of a POAMA run
## Started on a particular day
## And outputs some interested results for further analysis

## I care about ECLs, so main result is a daily timeseries of # of members that
## have an ECL on that day

## Also outputs a figure for each of the months in the file
## That shows the number of cyclones centred in that month
## Option to set a higher intensity threshold
library(ncdf4)
library(raster)
library(maps)
library(abind)
library(sp)
library(abind)

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
  par(mar = c(1, 1, 1, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
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

lat=seq(-89.5,89.5)
lon=seq(0,359.5)
makesmooth<-function(data,winwid=5)
{
  a=dim(data)
  m1=abind(data[(a[1]-winwid):a[1],a[2]:1],data[,a[2]:1],data[1:winwid,a[2]:1],along=1)
  m2=raster(t(m1),xmn=min(lon)-winwid,xmx=max(lon)+winwid,ymn=min(lat),ymx=max(lat))
  w2=focalWeight(m2,winwid,type="circle")
  m3=focal(m2,w2)
  tmp=t(as.matrix(m3)[a[2]:1,(winwid+1):(a[1]+winwid)])*pi*(winwid^2)
  return(tmp)
}


plot_freq_panel<-function(rad=5,lead=c(1,7),ethresh=0.15,closed=T,
       lonlim=c(0,360),latlim=c(-90,90),
       dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/",
       edir=NaN,efile=NaN,afile=NaN,
       fout="output")
{

lat=seq(-89.5,89.5)
lon=seq(0,359.5)  ### Can always combine into bigger cells later
mdays=c(31,28.25,31,30,31,30,31,31,30,31,30,31)
llist=seq(lead[1],lead[2])
mapname="world2"

breaks=c(-1,seq(-0.7,0.7,0.1),1)
breaks[9]=0
cols=col_anom(length(breaks)-1)

pdf(file=paste0(fout,".pdf"),width=16,height=10)
layout(rbind(1:4,5:8,9:12,rep(13,4)),height=c(1,1,1,0.25))

par(mar=c(2,2,4,1))

# First, get files
## Need to get ERAI for only the period of interest

years=seq(1990,2012)
erai=array(0,c(length(lon),length(lat),length(years),12))

for(y in 1:length(years))
{
fname=paste(edir,"/tracks_",years[y],".dat",sep="")
read.table(fname, sep="",skip=1)->fixes
colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Meh")
fixes$Year=floor(fixes$Date/10000)
if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
fixes$Date=10000*years[y]+fixes$Date%%10000
fixes$Month=floor(fixes$Date/100)%%100
fixes$Lat2=floor(fixes$Lat)
fixes$Lon2=floor(fixes$Lon)%%360
fixes$CV=abs(fixes$CV)
fixes=fixes[fixes$CV>=ethresh,]
if(closed) fixes=fixes[(fixes$Open==0 | fixes$Open==10),]
fixes$Date2=as.POSIXct(paste(fixes$Date,fixes$Time,sep=""),
                       format="%Y%m%d %H:%M",tz="GMT")

for(m in 1:12)
{
day1=as.POSIXct(paste0(years[y],sprintf("%2.2i",m),"01 00:00"),format="%Y%m%d %H:%M",tz="GMT")
I=which(fixes$Date2>=(day1+(lead[1]-1)*60*60*24) & fixes$Date2<(day1+lead[2]*60*60*24))
tmp=table(factor(fixes$Lon2[I],levels=0:359),factor(fixes$Lat2[I],levels=-90:89))
erai[,,y,m]=tmp
}
}
print(paste(mean(erai),max(erai)))

a=nc_open(paste0(dir,afile)) ## All, including open
access=ncvar_get(a,"cyclones")

## Reshape if necessary
if(lonlim[1]<0) {
   lon=seq(-179.5,179.5,1)
   erai=abind(erai[181:360,,,],erai[1:180,,,],along=1)
   access=abind(access[181:360,,,,],access[1:180,,,,],along=1)
   mapname="world"
   }

print(paste(mean(access),max(access)))

for(m in 1:12)
{
  accessM<-eraiM<-array(0,c(length(lon),length(lat),length(years)))
  for(y in 1:23) accessM[,,y]=makesmooth(apply(access[,,y,m,llist],c(1,2),sum),rad)
  for(y in 1:23) eraiM[,,y]=makesmooth(erai[,,y,m],rad)

  cyccorr<-array(NaN,c(length(lon),length(lat),2))

  meanfreq=abind(apply(eraiM,c(1,2),mean),apply(accessM,c(1,2),mean),along=3)
  print(apply(meanfreq,3,max,na.rm=T))
   for(i in 1:length(lon))
    for(j in 1:length(lat))
     if(!is.na(min(meanfreq[i,j,])))
     if(min(meanfreq[i,j,]>=0.5))
     {
      a=cor.test(accessM[i,j,],eraiM[i,j,],na.rm=T)
      cyccorr[i,j,1]=a$estimate
      cyccorr[i,j,2]=a$p.value
      }


  image(lon,lat,cyccorr[,,1],breaks=breaks,col=cols,
      xlab="",ylab="",xlim=lonlim,ylim=latlim,
      main=paste0(letters[m],") ACCESS-S anomcorr: ",month.name[m]," Days ",lead[1],"-",lead[2]))
  map(mapname,add=T)
   contour(lon,lat,cyccorr[,,2],levels=c(-100,0.05,100),add=T,lwd=2,col="black",drawlabels=F)
}

ColorBar(breaks,cols,subsampleg=1,vert=F)
dev.off()
}

plot_freq_panel(rad=5,lead=c(1,7),ethresh=0.25,closed=T,
        dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/",
        edir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/daily_proj240_lows_rad5cv0.15/",
        afile="ACCESS_globalcyclones_proj240_rad5cv0.25_40daylead.nc",
        fout=paste0("cyccorr_ACCESSvERAI_proj240_rad5cv0.25_globalrad5_lead1-7"))

plot_freq_panel(rad=5,lead=c(1,30),ethresh=0.25,closed=T,
        dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/",
        edir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/daily_proj240_lows_rad5cv0.15/",
        afile="ACCESS_globalcyclones_proj240_rad5cv0.25_40daylead.nc",
        fout=paste0("cyccorr_ACCESSvERAI_proj240_rad5cv0.25_globalrad5_lead1-30"))

plot_freq_panel(rad=10,lead=c(1,7),ethresh=0.25,closed=T,
        dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/",
        edir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/daily_proj240_lows_rad5cv0.15/",
        afile="ACCESS_globalcyclones_proj240_rad5cv0.25_40daylead.nc",
        fout=paste0("cyccorr_ACCESSvERAI_proj240_rad5cv0.25_globalrad10_lead1-7"))

plot_freq_panel(rad=10,lead=c(1,30),ethresh=0.25,closed=T,
        dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/",
        edir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/daily_proj240_lows_rad5cv0.15/",
        afile="ACCESS_globalcyclones_proj240_rad5cv0.25_40daylead.nc",
        fout=paste0("cyccorr_ACCESSvERAI_proj240_rad5cv0.25_globalrad10_lead1-30"))

