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
  par(mar = c(1, 1, 1, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols,
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(nstart-0.5, length(brks) - 1.5, subsampleg), tick = TRUE,
       labels = brks[seq(nstart, length(brks)-1, subsampleg)])
  } else {
    par(mar = c(1.5, 1, 1, 1), mgp = c(1.5, 0.3, 0), las = 1, cex = 1)
    image(1:length(cols), 1, t(t(1:length(cols))), axes = FALSE, col = cols,
          xlab = '', ylab = '')
    box()
    axis(1, at = seq(nstart-0.5, length(brks) - 1.5, subsampleg),
         labels = brks[seq(nstart=2, length(brks)-1, subsampleg)])
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

plot_freq_panel<-function(year1,year2,mthresh=c(0,100,200,300,400,600),
       seasons=rbind(c(5,10),c(11,4)),snames=c("MJJASO","NDJFMA"),
       dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/",reanal="ERAI",
       latlim=c(-90,90),lonlim=c(0,360),breaks=c(0,0.05,seq(0.1,1,0.1),1000),
       closed=F,cv=NaN,dur=NA,move=NA,fout="output")
{

years=seq(year1,year2,1)

vars=c("Frequency","Laplacian","Central Pressure (hPa)","Movement Speed (km/hr)")
vnames=c(NaN,"CV","MSLP","Move")
vbreaks=list(c(0,0.05,seq(0.1,1,0.1),1000),
             c(seq(0,0.4,0.05),100),
             c(seq(1000,1040,5),9999),
             c(seq(0,50,10),1000))

lat=seq(-89.5,89.5)
lon=seq(0,359.5)  ### Can always combine into bigger cells later

ss=length(snames)
if(ss<=2) vbreaks[[1]]=c(0,0.05,0.25,0.5,0.75,1,1.5,2,100)
vv=length(vars)

freq=array(0,c(length(lon),length(lat),length(years),ss))
varfreq=array(NaN,c(length(lon),length(lat),length(years),ss,vv,2))
varfreq[,,,,,1]=0
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

  for(s in 1:ss)
  {
   if(seasons[s,2]>=seasons[s,1]) mlist=seq(seasons[s,1],seasons[s,2]) else mlist=c(seq(seasons[s,1],12),seq(1,seasons[s,2]))

   I=which(fixes$Month%in%mlist)
   if(length(I)>0) tmp=table(factor(fixes$Lon2[I],levels=0:359),factor(fixes$Lat2[I],levels=-90:89))
   freq[,,y,s]=tmp

   for(v in 2:vv)
   {
    I=which(fixes$Month%in%mlist & !is.na(fixes[,vnames[v]]))
    if(length(I)>0) tmp=table(factor(fixes$Lon2[I],levels=0:359),factor(fixes$Lat2[I],levels=-90:89))
    varfreq[,,y,s,v,1]=tmp

    tmp2=aggregate(fixes[I,vnames[v]],by=list(fixes$Lon2[I],fixes$Lat2[I]),FUN=mean)
    jlat=match(tmp2[,2],seq(-90,89))
    ilon=match(tmp2[,1],seq(0,359))
    for(nn in which(!is.na(jlat))) varfreq[ilon[nn],jlat[nn],y,s,v,2]=as.numeric(tmp2[nn,3])
   }
  }
}

   if(lonlim[1]<0) {
   lon=seq(-179.5,179.5,1)
   library(abind)
   freq=abind(freq[181:360,,,],freq[1:180,,,],along=1)
   varfreq=abind(varfreq[181:360,,,,,],varfreq[1:180,,,,,],along=1)
  }

  pnum=1

  pdf(file=paste0(fout,".pdf"),width=(4*ss)+1.2,height=(2.7*vv))
  tmp=matrix(0,vv,ss+1)
  n=1
  for(s in 1:ss)
  for(i in 1:vv)
  {
  tmp[i,s]=n
  n=n+1
  }
  tmp[,ss+1]=seq(n,n+vv-1)

  layout(tmp,width=c(rep(1,ss),0.3))

  par(mar=c(2,2,4,1))

   for(s in 1:ss)
   {
   breaks=vbreaks[[1]]
   col1=col_val(length(breaks)-1)
   meanfreq=apply(freq[,,,s],c(1,2),mean,na.rm=T)
   meanfreq2=makesmooth(meanfreq)
   tit=paste0(snames[s],": Mean ",vars[1])
   image(lon,lat,meanfreq,breaks=breaks,col=col1,xlab="",ylab="",xlim=lonlim,ylim=latlim,
          main=paste0(letters[pnum],") ",tit))
   if(lonlim[1]<0) map('world',add=T) else map('world2',add=T)
   contour(lon,lat,meanfreq2,levels=breaks[seq(3,length(breaks),2)],add=T,lwd=2,col="black",drawlabels=F)
   pnum=pnum+1
   
   for(v in 2:vv)
   {
   print(mean(varfreq[,,,s,v,2],na.rm=T))
   breaks=vbreaks[[v]]
   col1=col_val(length(breaks))
   col1=col1[-1]

   meanvar=apply(varfreq[,,,s,v,2]*varfreq[,,,s,v,1],c(1,2),sum,na.rm=T)/apply(varfreq[,,,s,v,1],c(1,2),sum,na.rm=T)
   meanvar[apply(freq[,,,s],c(1,2),mean,na.rm=T)<0.05]=NaN
   meanvar2=makesmooth(apply(varfreq[,,,s,v,2]*varfreq[,,,s,v,1],c(1,2),sum,na.rm=T))/makesmooth(apply(varfreq[,,,s,v,1],c(1,2),sum,na.rm=T))
   meanvar2[makesmooth(apply(freq[,,,s],c(1,2),mean,na.rm=T))<0.05]=NaN

   tit=paste0(snames[s],": Mean ",vars[v])
   print(mean(meanvar,na.rm=T))
   image(lon,lat,meanvar,breaks=breaks,col=col1,xlab="",ylab="",xlim=lonlim,ylim=latlim,
          main=paste0(letters[pnum],") ",tit))
   if(lonlim[1]<0) map('world',add=T) else map('world2',add=T)   
   contour(lon,lat,meanvar2,levels=breaks,add=T,lwd=2,col="black",drawlabels=F)
   pnum=pnum+1
   }
  }

breaks=vbreaks[[1]]
col1=col_val(length(breaks)-1)
ColorBar(breaks,col1,subsampleg=1,vert=T) 

for(v in 2:vv)
{
breaks=vbreaks[[v]]
col1=col_val(length(breaks))
col1=col1[-1]
ColorBar(breaks,col1,subsampleg=2,vert=T,nstart=1)
}
dev.off()

}


slist2=c("MAM","JJA","SON","DJF")
smons2=rbind(c(3,5),c(6,8),c(9,11),c(12,2))
slist=c("MJJASO","NDJFMA")
smons=rbind(c(5,10),c(11,4))

plot_freq_panel(1980,2016,seasons=smons,snames=slist,breaks=breaks,closed=T,move=500,
        dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_highs_rad10cv0.075/",latlim=c(-50,-10),lonlim=c(100,180),
        reanal="ERAI",fout="paperfig_anticycfreq_ERAI_rad10cv0.075_500kmstats_2seasons_australia")

#plot_freq_panel(1980,2016,seasons=smons,snames=slist,breaks=breaks,closed=T,move=500,
#        dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_lows_rad5cv0.15/",
#        reanal="ERAI",fout="paperfig_cycfreq_ERAI_rad5cv0.15_500kmstats_vseason")

plot_freq_panel(1980,2016,seasons=smons2,snames=slist2,breaks=breaks,closed=T,move=500,
        dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_highs_rad10cv0.075/",latlim=c(-50,-10),lonlim=c(100,180),
        reanal="ERAI",fout="paperfig_anticycfreq_ERAI_rad10cv0.075_500kmstats_4seasons_australia")

#plot_freq_panel(1980,2016,seasons=smons2,snames=slist2,breaks=breaks,closed=T,move=500,
#        dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_lows_rad5cv0.15/",
#        reanal="ERAI",fout="paperfig_cycfreq_ERAI_rad5cv0.15_500kmstats_4seasons")
