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

ColorBar <- function(brks,cols,vert=T,subsampleg=1,nstart=2)
{
  if(vert) {
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols,
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(nstart-0.5, length(brks) - 1.5, subsampleg), tick = TRUE,
       labels = brks[seq(nstart, length(brks)-1, subsampleg)])
  } else {
    par(mar = c(1.5, 3, 1, 3), mgp = c(1.5, 0.3, 0), las = 1, cex = 0.9)
    image(1:length(cols), 1, t(t(1:length(cols))), axes = FALSE, col = cols,
          xlab = '', ylab = '')
    box()
    axis(1, at = seq(nstart-0.5, length(brks) - 1.5, subsampleg),
         labels = brks[seq(nstart, length(brks)-1, subsampleg)])
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
  tmp=t(as.matrix(m3)[a[2]:1,(winwid+1):(a[1]+winwid)])
  return(tmp)
}

plot_freq_panel<-function(year1,year2,season=c(1,12),
       dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/",reanal="ERAI",
       latlim=c(-90,90),lonlim=c(0,360),breaks=c(0,0.05,seq(0.1,1,0.1),1000),
       closed=F,cv=NaN,dur=NA,move=NA,fout="output",type="high")
{

years=seq(year1,year2,1)

vars=c(expression("a) Mean Laplacian (hPa (deg.lat)"^-2 * ")"),"b) Mean Central Pressure (hPa)","c) Mean Radius (degrees)","d) Mean Movement Speed (km/hr)")
vnames=c("CV","MSLP","Radius","Move")

if(type=="high" | type=="anticyclone") {
#vbreaks=list(c(0,seq(0.1,0.4,0.05),100),
#             c(0,seq(1020,1028,2),9999),
#             c(0,seq(5,10,1),1000),
#             c(0,seq(20,60,10),1000))
vbreaks=list(c(0,seq(0.1,0.3,0.05),100),
             c(0,seq(1020,1036,4),9999),
             c(0,seq(4,10,1),1000),
             c(0,seq(20,50,10),1000))
} else {
vbreaks=list(c(seq(0.2,0.9,0.1),100),
             c(seq(1020,970,-10),0),
             c(seq(3,7,0.5),1000),
             c(seq(0,50,10),1000))
}

lat=seq(-89.5,89.5)
lon=seq(0,359.5)  ### Can always combine into bigger cells later
if(season[2]>=season[1]) mlist=seq(season[1],season[2]) else mlist=c(seq(season[1],12),seq(1,season[2]))

#ss=length(mlist)
#if(ss<=2) vbreaks[[1]]=c(0,0.05,0.25,0.5,0.75,1,1.5,2,100)
vv=length(vars)

varfreq=array(NaN,c(length(lon),length(lat),length(years),vv,2))
varfreq[,,,,1]=0
for(y in 1:length(years))
{
  print(years[y])
  fname=paste(dir,"/tracks_",years[y],".dat",sep="")
  read.table(fname, sep="",skip=1)->fixes
  colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Depth","Radius","Up","Vp")
  fixes$Year=floor(fixes$Date/10000)
  if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
  fixes$Month=floor(fixes$Date/100)%%100
  fixes$Lat2=floor(fixes$Lat)
  fixes$Lon2=floor(fixes$Lon)%%360
  fixes$CV=abs(fixes$CV)
    fixes$Radius[fixes$Radius==0 | fixes$Radius > 100] = NaN
    fixes$Depth[fixes$Depth==0 | fixes$Depth > 100] = NaN

  fixes$Move<-NaN
  I=which(fixes$Fix>1)
  if(I[1]==1) I=I[-1]
  for(i in 1:length(I)) fixes$Move[I[i]]=spDistsN1(cbind(fixes$Lon[I[i]],fixes$Lat[I[i]]),cbind(fixes$Lon[I[i]-1],fixes$Lat[I[i]-1]),longlat=T)

  if(!is.na(move))
  {
    x<-rle(fixes$ID)
    events<-data.frame(ID=x$values,Length=x$lengths,Date1=rep(0,length(x$values)),Move=rep(0,length(x$values)))
    for(i in 1:length(events$ID))
    {
    I=which(fixes$ID==events[i,1])
    events$Move[i]=spDistsN1(cbind(fixes$Lon[min(I)],fixes$Lat[min(I)]),
                           cbind(fixes$Lon[max(I)],fixes$Lat[max(I)]),longlat=T)
    }
  events=events[events$Move>=move,]
  include<-match(fixes[,1],events[,1])
  J<-which(is.na(include)==0)
  fixes=fixes[J,]
  }

  if(!is.na(dur))
  {
  x<-rle(fixes$ID)
  events<-cbind(x$values,x$lengths,matrix(data=0,nrow=length(x$values),ncol=1))
  events=events[events[,2]>=dur,]
  include<-match(fixes[,1],events[,1])
  J<-which(is.na(include)==0)
  fixes=fixes[J,]
  }

  fixes$Move=fixes$Move/6 # Speed per hour
#  fixes=fixes[!is.na(fixes$Move),]
  if(!is.na(cv)) fixes=fixes[fixes$CV>=cv,]
  if(closed) fixes=fixes[(fixes$Open==0 | fixes$Open==10),] # Remove all open systems 

   for(v in 1:vv)
   {
    I=which(fixes$Month%in%mlist & !is.na(fixes[,vnames[v]]))
    if(length(I)>0) tmp=table(factor(fixes$Lon2[I],levels=0:359),factor(fixes$Lat2[I],levels=-90:89))
    varfreq[,,y,v,1]=tmp

    tmp2=aggregate(fixes[I,vnames[v]],by=list(fixes$Lon2[I],fixes$Lat2[I]),FUN=mean,na.rm=T)
    jlat=match(tmp2[,2],seq(-90,89))
    ilon=match(tmp2[,1],seq(0,359))
    for(nn in which(!is.na(jlat))) varfreq[ilon[nn],jlat[nn],y,v,2]=as.numeric(tmp2[nn,3])
  }
}

   if(lonlim[1]<0) {
   lon=seq(-179.5,179.5,1)
   library(abind)
   freq=abind(freq[181:360,,,],freq[1:180,,,],along=1)
   varfreq=abind(varfreq[181:360,,,,,],varfreq[1:180,,,,,],along=1)
  }

  pnum=1

  pdf(file=paste0(fout,".pdf"),width=10,height=6)
  layout(rbind(c(1,5,2,6),c(3,7,4,8)),width=c(1,0.25,1,0.25))

  par(mar=c(2,2,4,1))

   for(v in 1:vv)
   {
   breaks=vbreaks[[v]]
   #col1=col_val(length(breaks))
   library(RColorBrewer)
   col1=colorRampPalette(brewer.pal(9,"Blues"))(length(breaks))
   col1=col1[-1]

   meanvar=apply(varfreq[,,,v,2]*varfreq[,,,v,1],c(1,2),sum,na.rm=T)/apply(varfreq[,,,v,1],c(1,2),sum,na.rm=T)
   meanvar[apply(varfreq[,,,v,1],c(1,2),mean,na.rm=T)<0.05]=NaN
   meanvar2=makesmooth(apply(varfreq[,,,v,2]*varfreq[,,,v,1],c(1,2),sum,na.rm=T))/makesmooth(apply(varfreq[,,,v,1],c(1,2),sum,na.rm=T))
   meanvar2[makesmooth(apply(varfreq[,,,v,1],c(1,2),mean,na.rm=T))<0.05]=NaN

   print(range(meanvar,na.rm=T))

   print(mean(meanvar,na.rm=T))
   if(breaks[2]>breaks[1]) image(lon,lat,meanvar,breaks=breaks,col=col1,xlab="",ylab="",xlim=lonlim,ylim=latlim,main=vars[v],cex.main=1.5,cex.axis=1) else
     image(lon,lat,meanvar,breaks=rev(breaks),col=rev(col1),xlab="",ylab="",xlim=lonlim,ylim=latlim,main=vars[v],cex.main=1.5,cex.axis=1)
   if(lonlim[1]<0) map('world',add=T) else map('world2',add=T)   
   contour(lon,lat,meanvar2,levels=breaks,add=T,lwd=2,col="black")
   pnum=pnum+1
   }


for(v in 1:vv)
{
breaks=vbreaks[[v]]
col1=col_val(length(breaks))
col1=col1[-1]
ColorBar(breaks,col1,subsampleg=1,vert=T,nstart=2)
}
dev.off()

}


slist2=c("MAM","JJA","SON","DJF")
smons2=rbind(c(3,5),c(6,8),c(9,11),c(12,2))
slist=c("MJJASO","NDJFMA")
smons=rbind(c(5,10),c(11,4))

plot_freq_panel(1980,2016,season=c(1,12),breaks=breaks,closed=T,move=500,
        dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_highs_rad10cv0.075/",#latlim=c(-50,-10),lonlim=c(100,180),
        reanal="ERAI",fout="paperfig_anticycfreq_ERAI_rad10cv0.075_500stats_global2")

#plot_freq_panel(1980,2016,season=c(1,12),breaks=breaks,closed=T,move=500,type="cyclone",
#        dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_lows_rad5cv0.15/",
#        reanal="ERAI",fout="paperfig_cycfreq_ERAI_rad5cv0.15_500kmstats")

