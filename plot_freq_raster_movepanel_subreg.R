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


plot_freq_panel<-function(year1,year2,season=c(1,12),mthresh=c(0,100,200,300,400,600),
       dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/",reanal="ERAI",
       latlim=c(-90,90),lonlim=c(0,360),breaks=c(0,0.05,seq(0.1,1,0.1),1000),
       closed=F,cv=NaN,fout="output")
{

years=seq(year1,year2,1)
if(season[2]>=season[1]) mlist=seq(season[1],season[2]) else mlist=c(seq(season[1],12),seq(1,season[2]))

if(reanal=="20CR") {
varname="EnsMean"
fend="_EnsMean.nc" 
} else {
varname="systems"
fend=".nc"
}

lat=seq(-89.5,89.5)
lon=seq(0,359.5)  ### Can always combine into bigger cells later

ss=length(mthresh)
ss2=ceiling(ss/2)

freq=array(0,c(length(lon),length(lat),length(years),ss))
for(y in 1:length(years))
{
  print(years[y])
  fname=paste(dir,"/tracks_",years[y],".dat",sep="")
  read.table(fname, sep="",skip=1)->fixes
  colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Meh")
  fixes$Year=floor(fixes$Date/10000)
  if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
  fixes$Month=floor(fixes$Date/100)%%100
  fixes$Lat2=floor(fixes$Lat)
  fixes$Lon2=floor(fixes$Lon)%%360
  fixes$CV=abs(fixes$CV)

  fixes$Move<-NaN
  I=which(fixes$Fix>1)
  if(I[1]==1) I=I[-1]
  for(i in 1:length(I)) fixes$Move[I[i]]=spDistsN1(cbind(fixes$Lon[I[i]],fixes$Lat[I[i]]),cbind(fixes$Lon[I[i]-1],fixes$Lat[I[i]-1]),longlat=T)
  fixes$Move=fixes$Move/6 # Speed per hour

  fixes=fixes[!is.na(fixes$Move),]
  if(!is.na(cv)) fixes=fixes[fixes$CV>=cv,]
  if(closed) fixes=fixes[(fixes$Open==0 | fixes$Open==10),] # Remove all open systems 

  for(m in 1:(ss-1))
  {
  I=which(fixes$Move>=mthresh[m] & fixes$Move<mthresh[m+1] & fixes$Month%in%mlist)
  if(length(I)>0) tmp=table(factor(fixes$Lon2[I],levels=0:359),factor(fixes$Lat2[I],levels=-90:89))
  freq[,,y,m]=tmp
  }
  I=which(fixes$Move>=mthresh[ss] & fixes$Month%in%mlist)
  if(length(I)>0) tmp=table(factor(fixes$Lon2[I],levels=0:359),factor(fixes$Lat2[I],levels=-90:89))
  freq[,,y,ss]=tmp

}
   meanfreq=apply(freq,c(1,2,4),mean,na.rm=T)

  if(lonlim[1]<0) {
   lon=seq(-179.5,179.5,1)
   library(abind)
   meanfreq=abind(meanfreq[181:360,,],meanfreq[1:180,,],along=1)
  }

  col1=col_val(length(breaks)-1)
  pnum=1

  pdf(file=paste0(fout,".pdf"),width=2*ss,height=6.3)
  layout(rbind(seq(1,ss2),seq(ss2+1,ss),rep(ss+1,ss2)),height=c(1,1,0.3))
  par(mar=c(2,2,4,1))

   for(m in 1:ss)
   {
   meanfreq2=makesmooth(meanfreq[,,m])
   if(m==ss) tit=paste("Move >=",max(mthresh),"km/hr") else 
      tit=paste0("Move ",mthresh[m],"-",mthresh[m+1]," km/hr") 

   image(lon,lat,meanfreq[,,m],breaks=breaks,col=col1,xlab="",ylab="",xlim=lonlim,ylim=latlim,
          main=paste0(letters[pnum],") ",tit))
   if(lonlim[1]<0) map('world',add=T) else map('world2',add=T)   
   contour(lon,lat,meanfreq2,levels=breaks[seq(3,length(breaks),2)],add=T,lwd=2,col="black",drawlabels=F)
   pnum=pnum+1
  }

ColorBar(breaks,col1,subsampleg=1,vert=F)
dev.off()

  breaks=c(seq(0,40,5),100)
  col1=col_val(length(breaks)-1)
  pnum=1

  pdf(file=paste0(fout,"_percent.pdf"),width=2*ss,height=6.3)
  layout(rbind(seq(1,ss2),seq(ss2+1,ss),rep(ss+1,ss2)),height=c(1,1,0.3))
  par(mar=c(2,2,4,1))

   for(m in 1:ss)
   {
   meanfreq2=100*makesmooth(meanfreq[,,m])/makesmooth(apply(meanfreq,c(1,2),sum))
   I=which(apply(meanfreq,c(1,2),sum)<0.2)
   meanfreq2[I]=NaN

   if(m==ss) tit=paste("Move >=",max(mthresh),"km/hr") else
      tit=paste0("Move ",mthresh[m],"-",mthresh[m+1]," km/hr")

   image(lon,lat,meanfreq2,breaks=breaks,col=col1,xlab="",ylab="",xlim=lonlim,ylim=latlim,
          main=paste0(letters[pnum],") ",tit))
   if(lonlim[1]<0) map('world',add=T) else map('world2',add=T)
   contour(lon,lat,meanfreq2,levels=breaks[seq(3,length(breaks),2)],add=T,lwd=2,col="black",drawlabels=F)


   pnum=pnum+1
  }

ColorBar(breaks,col1,subsampleg=1,vert=F)
dev.off()

}

slist=c("","_MJJASO","_NDJFMA","_MAM","_JJA","_SON","_DJF")
smons=rbind(c(1,12),c(5,10),c(11,4),c(3,5),c(6,8),c(9,11),c(12,2))

for(s in 1:length(slist))
{
if(s==1) breaks=c(0,0.05,seq(0.1,1,0.1),1000) else breaks=c(seq(0,0.5,0.05),1000)

plot_freq_panel(1980,2016,mthresh=seq(0,50,10),season=smons[s,],breaks=breaks,
        dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_highs_rad10cv0.075/",
        reanal="ERAI",fout=paste0("paperfig_anticycfreq_ERAI_rad10cv0.075_movement",slist[s]))

plot_freq_panel(1980,2016,mthresh=seq(0,50,10),season=smons[s,],breaks=breaks,
        dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_lows_rad5cv0.15/",
        reanal="ERAI",fout=paste0("paperfig_cycfreq_ERAI_rad5cv0.15_movement",slist[s]))
}

