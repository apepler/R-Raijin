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

detrend<-function(x)
{
  regCoef=lm(x~seq(1,length(x)))
  return(resid(regCoef)) # Detrended
}

fieldsig<-function(x,alpha=0.05) # x is an array of p values
{
  ## Follows Wilks 2016
  x2=sort(as.vector(x)) ## Sorted vector
  n=alpha*seq(1,length(x2))/length(x2)
  pval=max(x2[x2<=n])
  y=array(0,dim(x))
  y[x<pval]=1
  return(y)
}

getInd<-function(iname)
{
ifile=switch(iname,
             SOI = "/short/eg3/asp561/Timeseries/SOI_18762016.csv",
             SAM = "/short/eg3/asp561/Timeseries/SAM_NOAA_19792016.csv",
             STRI = "/short/eg3/asp561/Timeseries/STRI_18902016.csv",
             STRP = "/short/eg3/asp561/Timeseries/STRP_18902016.csv",
             Hadley.SH = "/short/eg3/asp561/Timeseries/HadleyCell.SH.csv",
             Hadley.AsiaPac = "/short/eg3/asp561/Timeseries/HadleyCell.AsiaPac.csv",
             Hadley.NH = "/short/eg3/asp561/Timeseries/HadleyCell.NH.csv",
             Hadley.SH.Intensity = "/short/eg3/asp561/Timeseries/HadleyCell.SH.Intensity.csv",
              Hadley.NH.Intensity = "/short/eg3/asp561/Timeseries/HadleyCell.NH.Intensity.csv",
             NINO3.4 = "/short/eg3/asp561/Timeseries/NOAA.N34.OISST.monthly.csv",
             NINO3.4_ERSST = "/short/eg3/asp561/Timeseries/NOAA.N34.ERSST5.monthly.csv",
             DMI = "/short/eg3/asp561/Timeseries/NOAA.DMI.monthly.csv",
             AOI = "/short/eg3/asp561/Timeseries/NOAA.AOI.csv",
             NAO = "/short/eg3/asp561/Timeseries/NOAA.NAO.csv",
             PNA = "/short/eg3/asp561/Timeseries/NOAA.PNA.csv",
             IPO = "/short/eg3/asp561/Timeseries/Henley.IPO.ERSST4.csv",
             IPOf = "/short/eg3/asp561/Timeseries/Henley.IPO.ERSST4.filt.csv",
             NaN)
return(ifile)
}


plot_indcorr_panel<-function(year1,year2,seasons=rbind(c(5,10),c(11,4)),snames=c("MJJASO","NDJFMA"),iname="SOI",detrend=F,reanal="ERAI",
       dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/",type2="_500km",latlim=c(-90,90),lonlim=c(0,360),
       type="cyclone",proj="proj100_rad5cv0.15",fout="output.pdf")
{

years=seq(year1,year2,1)
if(year1<1979) startN=2 else startN=1

if(reanal=="20CR")
{
varname="EnsMean"
fend="_EnsMean.nc"
} else {
varname="systems"
fend=".nc"
}
lat=seq(-89.5,89.5)
lon=seq(0,359.5)  ### Can always combine into bigger cells later

ss=length(snames)

breaks=c(-1,seq(-0.7,-0.3,0.1),0,seq(0.3,0.7,0.1),1)
#breaks=c(1,seq(-0.7,0.7,0.1),1)
col1=col_anom(length(breaks)-1)

ss=length(snames)
pnum=1
pdf(file=fout,width=2*ss+1,height=5.5)
layout(rbind(c(seq(1,ss/2),ss+1),c(seq((ss/2)+1,ss),ss+1)),width=c(rep(1,ss/2),0.3))

par(mar=c(2,2,4,1))
n=0
  a=nc_open(paste0(dir,"/",reanal,"_UM_global",type,"s_",proj,type2,fend))
  systems=ncvar_get(a,varname)
  years2=ncvar_get(a,"year")

  I=which(years2%in%years)
  J=which(years%in%years2)
  freq=array(NaN,c(length(lon),length(lat),length(I)))  
  
  for(s in 1:ss)
  {
   if(seasons[s,2]>=seasons[s,1])
   {
    for(y in 1:length(I)) freq[,,y]=makesmooth(apply(systems[,,I[y],seasons[s,1]:seasons[s,2]],c(1,2),sum,na.rm=T),5)    
   } else {
    tmp=abind(systems[,,I[-length(I)],seasons[s,1]:12],systems[,,I[-1],1:seasons[s,2]],along=4)
    for(y in 1:(length(I)-1)) freq[,,y]=makesmooth(apply(tmp[,,y,],c(1,2),sum,na.rm=T),5)
   }

   meanfreq=apply(freq,c(1,2),mean,na.rm=T)
   cyccorr<-array(NaN,c(length(lon),length(lat),2))

    ifile=getInd(iname)
    if(is.na(ifile)) stop(paste("No file provided for",iname))
    index=read.csv(ifile)
    if(min(index$Year)>year1 | max(index$Year)<year2) stop("Index has shorter length of record than selected years")
    K=which(index$Year%in%years2[I])
    if(seasons[s,2]>=seasons[s,1])
    {
      index2=apply(index[K,(seasons[s,1]:seasons[s,2])+1],1,mean)
    } else {
      tmp=cbind(index[K[-length(K)],(seasons[s,1]+1):13],index[K[-1],2:(seasons[s,2]+1)])
      index2=c(apply(tmp,1,mean),NaN)
   }

   cyccorr<-array(NaN,c(length(lon),length(lat),2))
   n=n+1
   tit=paste0(letters[n],") ",snames[s]," correlation with ",iname)

   for(i in 1:length(lon))
    for(j in 1:length(lat))
     if(!is.na(meanfreq[i,j]))
     if(meanfreq[i,j]>=0.1)
     {
     if(detrend) a=cor.test(detrend(freq[i,j,]),detrend(index2),na.rm=T) else a=cor.test(freq[i,j,],index2,na.rm=T)
      cyccorr[i,j,2]=a$p.value
      if(a$p.value<0.05) cyccorr[i,j,1]=a$estimate # Only plot where individually significant
      }

   if(lonlim[1]<0) {
   lon2=seq(-179.5,179.5,1)
   cyccorr=abind(cyccorr[181:360,,],cyccorr[1:180,,],along=1)
   } else lon2=lon 

   image(lon2,lat,cyccorr[,,1],breaks=breaks,col=col1,xlab="",ylab="",
          xlim=lonlim,ylim=latlim,main=tit)
   map('world2',add=T)   
   pval=array(p.adjust(cyccorr[,,2],"fdr"),dim(cyccorr[,,2]))
#   pval=fieldsig(cyccorr[,,2],alpha=0.05)
   contour(lon2,lat,pval<0.05,levels=c(-100,0.05,100),add=T,lwd=1.5,col="black",drawlabels=F)
#   contour(lon,lat,(pval<0.05 & cyccorr[,,1]<0),levels=c(-100,0.05,100),add=T,lwd=1.5,lty=2,drawlabels=F)
  

}

ColorBar(breaks,col1,subsampleg=1,vert=T)
dev.off()
}

slist2=c("MAM","JJA","SON","DJF")
smons2=rbind(c(3,5),c(6,8),c(9,11),c(12,2))
slist=c("MJJASO","NDJFMA")
smons=rbind(c(5,10),c(11,4))

for(index in c("STRI","STRP"))
{

plot_indcorr_panel(1980,2016,seasons=smons,
       snames=slist,iname=index,
       dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf",reanal="ERAI",
       type="anticyclone",proj="proj100_rad10cv0.075",type2="_500km",latlim=c(-50,-10),lonlim=c(100,180),
       fout=paste0("paperfig_anticyccorr_",index,"_ERAIpanel_proj100_rad5cv0.15_500km_fieldsig2.pdf"))
}
