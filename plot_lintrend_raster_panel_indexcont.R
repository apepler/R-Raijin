# This is code that takes the 11 ensemble members of a POAMA run
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

detrend<-function(x,ind=seq(1,length(x)))
{
  regCoef=lm(x~ind)
  return(resid(regCoef)) # Detrended
}

getInd<-function(iname)
{
ifile=switch(iname,
             SOI = "/short/eg3/asp561/Timeseries/SOI_18762016.csv",
             SAM = "/short/eg3/asp561/Timeseries/SAM_NOAA_19792016.csv",
             SAM_M = "/short/eg3/asp561/Timeseries/sam_marshall.csv",
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


plot_indcorr_panel<-function(year1,year2,season=c(5,10),sname="MJJASO",ind="SAM_M",
       dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/",type2="_500km",
       type="cyclone",proj="proj100_rad5cv0.15",fout="output.pdf")
{

years=seq(year1,year2,1)
if(year1<1979) startN=2 else startN=1

reanals=c("ERAI","NCEP1","20CR")
varnames=c("systems","systems","EnsMean")
fend=c(".nc",".nc","_EnsMean.nc")

lat=seq(-89.5,89.5)
lon=seq(0,359.5)  ### Can always combine into bigger cells later

#breaks=c(-1,seq(-0.7,-0.3,0.1),0,seq(0.3,0.7,0.1),1)
#breaks=c(1,seq(-0.7,0.7,0.1),1)
breaks=c(-100,seq(-1.5,1.5,0.25),100)
col1=col_anom(length(breaks)-1)

if(year1<1979)
{
pdf(file=fout,width=8,height=9)
layout(cbind(c(1:3,7),c(4:6,7)),height=c(1,1,1,0.3))
} else {
pdf(file=fout,width=12,height=9)
layout(cbind(c(1:3,10),c(4:6,10),7:10),height=c(1,1,1,0.3))
}

par(mar=c(2,2,4,1))
for(n in startN:3)
{
  a=nc_open(paste0(dir,"/",reanals[n],"_UM_global",type,"s_",proj,type2,fend[n]))
  systems=ncvar_get(a,varnames[n])
  years2=ncvar_get(a,"year")

  I=which(years2%in%years)
  J=which(years%in%years2)
  freq=array(NaN,c(length(lon),length(lat),length(I)))  
  
   if(season[2]>=season[1])
   {
    for(y in 1:length(I)) freq[,,y]=makesmooth(apply(systems[,,I[y],season[1]:season[2]],c(1,2),sum,na.rm=T),5)    
   } else {
    tmp=abind(systems[,,I[-length(I)],season[1]:12],systems[,,I[-1],1:season[2]],along=4)
    for(y in 1:(length(I)-1)) freq[,,y]=makesmooth(apply(tmp[,,y,],c(1,2),sum,na.rm=T),5)
   }

   meanfreq=apply(freq,c(1,2),mean,na.rm=T)
   index=read.csv(getInd(ind))
   if(min(index$Year)>year1 | max(index$Year)<year2) stop("Index has shorter length of record than selected years")

   K=which(index$Year%in%years2[I])
    if(season[2]>=season[1])
    {
      index2=apply(index[K,(season[1]:season[2])+1],1,mean)
    } else {
      tmp=cbind(index[K[-length(K)],(season[1]+1):13],index[K[-1],2:(season[2]+1)])
      index2=c(apply(tmp,1,mean),NaN)
   }

   trends<-array(NaN,c(length(lon),length(lat),5,2))
   dimnames(trends)[[3]]<-titnames<-c("Trend","Detrended.reg","Index.trend",paste("Trend congruent with",iname),"Remainder")
   dimnames(trends)[[4]]<-c("Value","Significance")

   index_detrend=detrend(index2)

# Index Trend
   a=lm(index2~years)
   b=summary(a)$coefficients
   trends[,,3,1]=a$coefficients[2]
   trends[,,3,2]=b[2,4]

   for(i in 1:length(lon))
    for(j in 1:length(lat))
     if(!is.na(meanfreq[i,j]))
     if(meanfreq[i,j]>=0.1)
     {
      count_detrend=detrend(freq[i,j,])

# Frequency trend
      a=lm(freq[i,j,]~years)
      b=summary(a)$coefficients
      trends[i,j,1,1]=a$coefficients[2]
      trends[i,j,1,2]=b[2,4]

# Relationship
      a=lm(count_detrend~index_detrend)
      b=summary(a)$coefficients
      trends[i,j,2,1]=a$coefficients[2]
      trends[i,j,2,2]=b[2,4]
    }

# Contribution
  trends[,,4,1]=trends[,,2,1]*trends[,,3,1]
  trends[,,4,2]=apply(trends[,,2:3,2],c(1,2),max) # Only sig if both sig
  trends[,,5,1]=trends[,,1,1]-trends[,,4,1] # Contribution not included
 
## Now, make everything back into %/year

for(i in c(1,4,5))
{  
   tmp=100*trends[,,i,1]/meanfreq

   image(lon,lat,tmp,breaks=breaks,col=col1,xlab="",ylab="",
          main=paste(titnames[i],"-",reanals[n]))
   map('world2',add=T)   
   if(i<5) contour(lon,lat,trends[,,i,2],levels=c(-100,0.05,100),add=T,lwd=2,col="black",drawlabels=F)
 } 

}

ColorBar(breaks,col1,subsampleg=1,vert=F)
dev.off()
}

slist=c("MAM","JJA","SON","DJF","MJJASO","NDJFMA")
smons=rbind(c(3,5),c(6,8),c(9,11),c(12,2),c(5,10),c(11,4))
inames=c("SAM_M","STRI","STRP")

for(s in 5:6)
 for(iname in inames)
{


plot_indcorr_panel(1960,2014,season=smons[s,],sname=slist,ind=iname,
       dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf",
       type="anticyclone",proj="proj100_rad10cv0.075",type2="_500km",
       fout=paste0("paperfig_anticyclintrend_",slist[s],"6014_cont",iname,"_proj100_rad10cv0.075_500km.pdf"))

plot_indcorr_panel(1960,2014,season=smons[s,],sname=slist,ind=iname,
       dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf",
       type="cyclone",proj="proj100_rad5cv0.15",type2="_500km",
       fout=paste0("paperfig_cyclintrend_",slist[s],"6014_cont",iname,"_proj100_rad5cv0.15_500km.pdf"))
}



