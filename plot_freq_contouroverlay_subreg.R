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
#col_val <- color.palette(c("white","blue","darkblue","black"),c(20,10,5))
col_val <- color.palette(c("white","cyan","blue","darkblue"), c(10,20,5))

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
makesmooth<-function(data,winwid=2,makesum=F)
{
  a=dim(data)
  
  if(lat[2]<lat[1]) ll=a[2]:1 else ll=1:a[2]

  m1=abind(data[(a[1]-winwid+1):a[1],ll],data[,ll],data[1:winwid,ll],along=1) 
  m2=raster(t(m1),xmn=min(lon)-winwid,xmx=max(lon)+winwid,ymn=min(lat),ymx=max(lat))
  w2=focalWeight(m2,winwid,type="circle") 
  if(makesum) w2[w2>0]=1
  m3=focal(m2,w2)
  tmp=t(as.matrix(m3)[ll,(winwid+1):(a[1]+winwid)])
  return(tmp)
}



plot_freq_panel<-function(year1,year2,season=c(1,12),
       dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/",reanal="ERAI",latlim=c(-90,90),lonlim=c(0,360),ozmap=F,
       type="cyclone",proj="proj100_rad5cv0.15",pnames="rad5cv0.15",fout="output",breaks=c(0,0.05,seq(0.25,1,0.25),1.5,2,1000))
{

years=seq(year1,year2,1)

if(year1<1979) startN=2 else startN=1

if(reanal=="20CR") {
varname="EnsMean"
fend="_EnsMean.nc" 
} else {
varname="systems"
fend=".nc"
}

lat=seq(-89.5,89.5)
lon=seq(0,359.5)  ### Can always combine into bigger cells later

ss=length(pnames)
#col1=col_val(length(breaks)-1)
library(viridisLite)
library(oz)
col1=viridis(length(pnames)+1)[1:length(pnames)]
#col1=c("black","darkgray")
#col1=c("blue","red")

pnum=1
pdf(file=fout,width=4,height=3)
par(mar=rep(0.2,4))

map('world2',xlim=lonlim,ylim=latlim,xlab="",ylab="",
    main=paste0("Mean ",type," frequency"))
axis(1,cex.axis=0.7)
axis(2,cex.axis=0.7)
if(ozmap) oz(states=T,add=T)
box()

for(n in 1:length(pnames))
{
  a=nc_open(paste0(dir,"/",reanal,"_UM_global",type,"s_",proj[n],fend))
  systems=ncvar_get(a,varname)
  years2=ncvar_get(a,"year")

  I=which(years2%in%years)
  J=which(years%in%years2)
  freq=array(NaN,c(length(lon),length(lat),length(I)))  

  if(season[2]>=season[1])
  {
    for(y in 1:length(I)) freq[,,y]=apply(systems[,,I[y],season[1]:season[2]],c(1,2),sum,na.rm=T)    
   } else {
    tmp=abind(systems[,,I[-length(I)],season[1]:12],systems[,,I[-1],1:season[2]],along=4)
    for(y in 1:(length(I)-1)) freq[,,y]=apply(tmp[,,y,],c(1,2),sum,na.rm=T)
   }

   if(lonlim[1]<0) {
   lon2=seq(-179.5,179.5,1)
   meanfreq=abind(meanfreq[181:360,],meanfreq[1:180,],along=1)
   meanfreq2=makesmooth(meanfreq,makesum=T,winwid=2)
   contour(lon2,lat,makesmooth(meanfreq2,winwid=2),levels=breaks,add=T,lwd=3,col=col1[n],labcex=0.8)
   } else {
   meanfreq=apply(freq,c(1,2),mean,na.rm=T)
   meanfreq2=makesmooth(meanfreq,makesum=T,winwid=2)
   contour(lon,lat,makesmooth(meanfreq2,winwid=2),levels=breaks,add=T,lwd=3,col=col1[n],labcex=0.8)
   }

}

legend("topright",pnames,lwd=3,col=col1,ncol=2,bg="white")
dev.off()
}

#plot_freq_panel(1980,2016,season=c(1,12),breaks=c(0,0.05,seq(0.2,2,0.2),1000),
#        dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf",
#        type="cyclone",reanal="ERAI",pnames=c("MSLP","850hPa","700hPa","500hPa"),
#        proj=c("proj100_rad5cv0.15_500km",paste0(c(850,700,500),"hPa_proj100_rad5cv1_500km")),latlim=c(-50,-10),lonlim=c(100,180),
#        fout="paperfig_cycfreqvheight_500km_ERAI_aust.pdf")

#plot_freq_panel(1980,2016,season=c(1,12),breaks=c(0,0.05,seq(0.1,1,0.1),1000),
#        dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf",
#        type="cyclone",reanal="ERAI",pnames=c("MSLP","500hPa"),
#        proj=c("proj100_rad2cv1_D2",paste0(500,"hPa_proj100_rad5cv4_D2")),latlim=c(-50,-10),lonlim=c(100,180),
#        fout="paperfig_cycfreqvheight_500km_ERAI_aust_rad2cv1_500hParad5cv4.pdf")

#plot_freq_panel(1980,2016,season=c(1,12),breaks=c(0,0.05,seq(0.25,3,0.25),1000),
#        dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf",
#        type="anticyclone",reanal="ERAI",pnames=c("MSLP","850hPa","700hPa","500hPa"),
#        proj=c("proj100_rad10cv0.075_500km",paste0(c(850,700,500),"hPa_proj100_rad10cv0.5_500km")),latlim=c(-50,-10),lonlim=c(100,180),
#        fout="paperfig_anticycfreqvheight_500km_ERAI_aust2.pdf")

#plot_freq_panel(1980,2016,season=c(1,12),breaks=c(0,0.05,seq(0.5,3,0.5),1000),
#        dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf",
#        type="anticyclone",reanal="ERAI",pnames=c("All MSLP anticyclones","Length>=6 hours","500km movement","All 850hPa anticyclones"),
#        proj=c("proj100_rad10cv0.075","proj100_rad10cv0.075_D2","proj100_rad10cv0.075_500km","850hPa_proj100_rad10cv0.5"),#latlim=c(-50,-10),lonlim=c(100,180),
#        fout="paperfig_anticycfreqvduration_ERAI_global_Bu.pdf")

plot_freq_panel(1980,2016,season=c(1,12),breaks=c(1,5,10,15,20),
        dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf",
        type="cyclone",reanal="ERAI",pnames=c("MSLP","500 hPa"),
        proj=c("proj240_rad2cv1.5",paste0(500,"hPa_proj240_rad2cv12")),
        latlim=c(-45,-20),lonlim=c(130,180),ozmap=T,
        fout="cyclonefreq_contours_rad2vheight_MSLP500_ERAI_SEA_viridis.pdf")

