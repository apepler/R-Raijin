## This is code that takes the 11 ensemble members of a POAMA run
## Started on a particular day
## And outputs some interested results for further analysis

## I care about ECLs, so main result is a daily timeseries of # of members that
## have an ECL on that day

## Also outputs a figure for each of the months in the file
## That shows the number of cyclones centred in that month
## Option to set a higher intensity threshold

rm(list=ls())

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
  par(mar = c(2, 1, 2, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
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

makesmooth<-function(data,winwid=5,lon=seq(1,dimsizes(data)[1]),lat=seq(1,dimsizes(data)[2]))
{
  a=dim(data)
  m1=abind(data[(a[1]-winwid+1):a[1],a[2]:1],data[,a[2]:1],data[1:winwid,a[2]:1],along=1)
  m2=raster(t(m1),xmn=min(lon)-winwid,xmx=max(lon)+winwid,ymn=min(lat),ymx=max(lat))
  w2=focalWeight(m2,winwid,type="circle")
  m3=focal(m2,w2,pad=T)
  tmp=t(as.matrix(m3)[a[2]:1,(winwid+1):(a[1]+winwid)])
  return(tmp)
}

plot_indcorr_panel<-function(year1,year2,seasons=rbind(c(5,10),c(11,4)),snames=c("MJJASO","NDJFMA"),detrend=F,
       dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/",
       type="low",thresh=0,closed=T,move=NA,dur=NA,fout="output.pdf",
       latlim=c(-90,90),lonlim=c(0,360))
{
years=seq(year1,year2,1)

## First, get the slp and make an slp grid

a=nc_open("/g/data/eg3/asp561/ERAI/ERAI.slp.monmean.19792016.nc")
lat=ncvar_get(a,"lat")
lon=ncvar_get(a,"lon")
slp=ncvar_get(a,"psl")/100

a=nc_open("/g/data/eg3/asp561/ERAI/ERAI.str.mon.mean.longitude.nc")
strnames=c("STRI_SH","STRI_NH","STRP_SH","STRP_NH")
str=abind(ncvar_get(a,"STRI_SH"),ncvar_get(a,"STRI_NH"),ncvar_get(a,"STRP_SH"),ncvar_get(a,"STRP_NH"),along=3)
str=abs(str) # So that position means the same thing in both hemispheres - +ve = poleward

dd=dim(slp)
dd2=floor(dd[3]/12)
years2=seq(1979,by=1,length.out=dd2)

slp2=array(NaN,c(dd[1],dd[2],dd2,12))
for(n in 1:dd2) slp2[,,n,]=slp[,,(n*12-11):(n*12)]
slp2=slp2[,,years2%in%years,]

str2=array(NaN,c(dd[1],dd2,12,4))
for(n in 1:dd2) str2[,n,,]=str[,(n*12-11):(n*12),]
str2=str2[,years2%in%years,,]

if(lat[2]<lat[1])
{
  lat=lat[length(lat):1]
  slp2=slp2[,length(lat):1,,]
}
if(min(lon)<0)
{
I=which(lon<0)
lon=c(lon[-I],lon[I]+360)
slp2=abind(slp2[-I,,,],slp2[I,,,],along=1)
str2=abind(str2[-I,,,],str2[I,,,],along=1)
}

dlat=abs(lat[2]-lat[1])
dlon=abs(lon[2]-lon[1])
lat2=lat/dlat
lon2=lon/dlat

### Next, get the cyclones for that grid point

systems<-array(0,c(length(lon),length(lat),length(years),12))

for(y in 1:length(years))
{
print(years[y])
fname=paste(dir,"/tracks_",years[y],".dat",sep="")
read.table(fname, sep="",skip=1)->fixes
colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Depth","Radius","Up","Vp")
fixes$Year=floor(fixes$Date/10000)
fixes$CV=abs(fixes$CV)
fixes$Depth=abs(fixes$Depth)
if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
fixes$Month=floor(fixes$Date/100)%%100
fixes$Lat2=round(fixes$Lat/dlat,0) ## Let's skip the first step, and just make it the 10 degrees
fixes$Lon2=round((fixes$Lon%%360)/dlon,0)

### Make table of events to combine with DJF for exclusion
 if(!is.na(move))
 {
    fixes$Move<-NaN
    I=which(fixes$Fix>1)
    if(I[1]==1) I=I[-1]
    for(i in 1:length(I)) fixes$Move[I[i]]=spDistsN1(cbind(fixes$Lon[I[i]],fixes$Lat[I[i]]),cbind(fixes$Lon[I[i]-1],fixes$Lat[I[i]-1]),longlat=T)

    x<-rle(fixes$ID)
    events<-data.frame(ID=x$values,Length=x$lengths,Date1=rep(0,length(x$values)),Move=rep(0,length(x$values)))
    for(i in 1:length(events$ID))
    {
    events$Date1[i]=min(fixes$Date[fixes$ID==events$ID[i]])
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
 fixes=fixes[fixes$CV>=thresh,]
 if(closed) fixes=fixes[(fixes$Open==0 | fixes$Open==10),] # Remove all open systems
 tmp=table(factor(fixes$Lon2,levels=lon2),factor(fixes$Lat2,levels=lat2),fixes$Month)
 systems[,,y,]=tmp
 print(mean(systems[,,y,],na.rm=T))

} # End year loop

print('Cyclones done')

## Now, need to make pdf & loop through seasons

ss=length(snames)
inames=c("SLP","STRI","STRP")
inum=3
breaks=c(-10,seq(-7,7,1),10)/10
col1=col_anom(length(breaks)-1)
pnum=1

tmp=matrix(0,inum+1,ss)
n=1
for(s in 1:ss)
for(i in 1:inum)
{
tmp[i,s]=n
n=n+1
}
tmp[inum+1,]=n

pdf(file=fout,width=4*ss,height=(2.7*inum)+0.8)
layout(tmp,height=c(rep(1,inum),0.3))

par(mar=c(2,2,4,1))
cyccorr<-array(NaN,c(length(lon),length(lat),3,ss,2))
n=0
for(s in 1:ss)
  {
   if(seasons[s,2]>=seasons[s,1])
   {
    freq<-array(NaN,c(length(lon),length(lat),length(years)))
    for(y in 1:length(years)) freq[,,y]=makesmooth(apply(systems[,,y,seasons[s,1]:seasons[s,2]],c(1,2),sum,na.rm=T),winwid=6,lon,lat)
    str=apply(str2[,,seasons[s,1]:seasons[s,2],],c(1,2,4),mean,na.rm=T)
    slp=apply(slp2[,,,seasons[s,1]:seasons[s,2]],c(1,2,3),mean,na.rm=T)
   } else {
    freq<-array(NaN,c(length(lon),length(lat),length(years)-1))
    tmp=abind(systems[,,-length(years),seasons[s,1]:12],systems[,,-1,1:seasons[s,2]],along=4)
    for(y in 1:(length(years)-1)) freq[,,y]=makesmooth(apply(tmp[,,y,],c(1,2),sum,na.rm=T),winwid=6,lon,lat)
    tmp=abind(str2[,-length(years),seasons[s,1]:12,],str2[,-1,1:seasons[s,2],],along=3)
    str=apply(tmp,c(1,2,4),mean,na.rm=T)
    tmp=abind(slp2[,,-length(years),seasons[s,1]:12],slp2[,,-1,1:seasons[s,2]],along=4)
    slp=apply(tmp,c(1,2,3),mean,na.rm=T)
   }

### Correlations
   meanfreq=apply(freq,c(1,2),mean,na.rm=T)

for(i in 1:length(lon))
  for(j in 1:length(lat))
    if(!is.na(meanfreq[i,j]))
     if(meanfreq[i,j]>=0.1)
     {
      if(detrend) a=cor.test(detrend(freq[i,j,]),detrend(slp[i,j,]),na.rm=T) else a=cor.test(freq[i,j,],slp[i,j,],na.rm=T)
      cyccorr[i,j,1,s,1]=a$estimate
      cyccorr[i,j,1,s,2]=a$p.value

     if(lat[j]<0)
     {
      if(detrend) a=cor.test(detrend(freq[i,j,]),detrend(str[i,,1]),na.rm=T) else a=cor.test(freq[i,j,],str[i,,1],na.rm=T)
      cyccorr[i,j,2,s,1]=a$estimate
      cyccorr[i,j,2,s,2]=a$p.value      
      if(detrend) a=cor.test(detrend(freq[i,j,]),detrend(str[i,,3]),na.rm=T) else a=cor.test(freq[i,j,],str[i,,3],na.rm=T)
      cyccorr[i,j,3,s,1]=a$estimate
      cyccorr[i,j,3,s,2]=a$p.value
     } else {
      if(detrend) a=cor.test(detrend(freq[i,j,]),detrend(str[i,,2]),na.rm=T) else a=cor.test(freq[i,j,],str[i,,2],na.rm=T)
      cyccorr[i,j,2,s,1]=a$estimate
      cyccorr[i,j,2,s,2]=a$p.value
      if(detrend) a=cor.test(detrend(freq[i,j,]),detrend(str[i,,4]),na.rm=T) else a=cor.test(freq[i,j,],str[i,,4],na.rm=T)
      cyccorr[i,j,3,s,1]=a$estimate
      cyccorr[i,j,3,s,2]=a$p.value
     }
    }

 for(ii in 1:3)
 {
  n=n+1
  tit=paste0(letters[n],") ",snames[s]," correlation with ",inames[ii])
  tmp=cyccorr[,,ii,s,1]
  tmp[cyccorr[,,ii,s,2]>=0.05]=NaN
  image(lon,lat,tmp,breaks=breaks,col=col1,xlab="",ylab="",
          main=tit)
  map('world2',add=T)
  pval=array(p.adjust(cyccorr[,,ii,s,2],"fdr"),dim(cyccorr[,,ii,s,2]))
  contour(lon,lat,pval,levels=c(-100,0.05,100),add=T,lwd=2,col="black",drawlabels=F)  
 }
}

print('Corrs done')

ColorBar(breaks,col1,subsampleg=2,vert=F)
dev.off()

}

slist2=c("MAM","JJA","SON","DJF")
smons2=rbind(c(3,5),c(6,8),c(9,11),c(12,2))
slist=c("MJJASO","NDJFMA")
smons=rbind(c(5,10),c(11,4))

plot_indcorr_panel(1980,2016,seasons=smons,snames=slist,
       dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_highs_rad5cv0.075/",
       type="high",closed=T,move=500,
       fout=paste0("paperfig_anticyccorr_strslp_ERAIpanel_proj100_rad10cv0.075_500km_2seasons_fieldsig.pdf"))


plot_indcorr_panel(1980,2016,seasons=smons2,snames=slist2,
       dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_highs_rad5cv0.075/",
       type="high",closed=T,move=500,
       fout=paste0("paperfig_anticyccorr_strslp_ERAIpanel_proj100_rad10cv0.075_500km_4seasons_fieldsig.pdf"))

