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


plot_freq_panel<-function(year1,year2,seasons=rbind(c(5,10),c(11,4)),snames=c("MJJASO","NDJFMA"),
       dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/",reanal="ERAI",boxes=F,
       latlim=c(-90,90),lonlim=c(0,360),breaks=c(0,0.05,seq(0.1,1,0.1),1000),
       type="cyclone",proj="proj100_rad5cv0.15",fout="output")
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

ss=length(snames)
col1=col_val(length(breaks)-1)
pnum=1
pdf(file=fout,width=10,height=ss*1.5)
if(ss==2) layout(cbind(1,2,3),width=c(1,1,0.25)) else
  layout(rbind(c(seq(1,ss/2),ss+1),c(seq((ss/2)+1,ss),ss+1)),width=c(1,1,0.25))

par(mar=c(2,2,4,1))

  a=nc_open(paste0(dir,"/",reanal,"_UM_global",type,"s_",proj,fend))
  systems=ncvar_get(a,varname)
  years2=ncvar_get(a,"year")

  I=which(years2%in%years)
  J=which(years%in%years2)
  freq=array(NaN,c(length(lon),length(lat),length(I)))  

  for(s in 1:ss)
  {
   if(seasons[s,2]>=seasons[s,1])
   {
    for(y in 1:length(I)) freq[,,y]=apply(systems[,,I[y],seasons[s,1]:seasons[s,2]],c(1,2),sum,na.rm=T)    
   } else {
    tmp=abind(systems[,,I[-length(I)],seasons[s,1]:12],systems[,,I[-1],1:seasons[s,2]],along=4)

    for(y in 1:(length(I)-1)) freq[,,y]=apply(tmp[,,y,],c(1,2),sum,na.rm=T)
   }

   meanfreq=apply(freq,c(1,2),mean,na.rm=T)
   meanfreq2=makesmooth(meanfreq)

   if(lonlim[1]<0) {
   lon2=seq(-179.5,179.5,1)
   library(abind)   
   meanfreq=abind(meanfreq[181:360,],meanfreq[1:180,],along=1)
   meanfreq2=makesmooth(meanfreq)
   image(lon2,lat,meanfreq,breaks=breaks,col=col1,xlab="",ylab="",xlim=lonlim,ylim=latlim,
          main=paste0(letters[pnum],") ",snames[s]," mean ",type," frequency"))
   map('world',add=T)
   contour(lon2,lat,meanfreq2,levels=breaks[seq(3,length(breaks),2)],add=T,lwd=2,col="black",drawlabels=F)
   } else {
   image(lon,lat,meanfreq,breaks=breaks,col=col1,xlab="",ylab="",xlim=lonlim,ylim=latlim,
          main=paste0(letters[pnum],") ",snames[s]," mean ",type," frequency"))
   map('world2',add=T)   
   contour(lon,lat,meanfreq2,levels=breaks[seq(3,length(breaks),2)],add=T,lwd=2,col="black",drawlabels=F)

   if(boxes) 
   {
#    polygon(x=c(110,155,155,110,110),y=c(-25,-25,-40,-40,-25),border="black",lwd=3)
#    polygon(x=c(135,147.5,147.5,135,135),y=c(-32.5,-32.5,-42.5,-42.5,-32.5),border="red",lwd=3)
#    polygon(x=c(160,147.5,147.5,160,160),y=c(-25,-25,-40,-40,-25),border="orange",lwd=3)
#    polygon(x=c(110,122.5,122.5,110,110),y=c(-40,-40,-27.5,-27.5,-40),border="purple",lwd=3)
#    if(s==1) legend("topleft",legend=c("S.Aust","SEA","Tasman","SWWA"),col=c("black","red","orange","purple"),lwd=3,ncol=2,bg="white")

#    polygon(x=c(90,110,110,90,90),y=c(-42.5,-42.5,-27.5,-27.5,-42.5),border="red",lwd=3)
#    polygon(x=c(120,145,145,120,120),y=c(-42.5,-42.5,-27.5,-27.5,-42.5),border="orange",lwd=3)
    polygon(x=c(90,115,115,90,90),y=c(-42.5,-42.5,-27.5,-27.5,-42.5),border="red",lwd=3)
    polygon(x=c(115,150,150,115,115),y=c(-42.5,-42.5,-27.5,-27.5,-42.5),border="orange",lwd=3)
    polygon(x=c(150,170,170,150,150),y=c(-42.5,-42.5,-27.5,-27.5,-42.5),border="purple",lwd=3)
    if(s==1) legend("topleft",legend=c("SEIO","SAust","Tasman"),col=c("red","orange","purple"),lwd=3,bg="white")

   }

   }
   pnum=pnum+1
  }


ColorBar(breaks,col1,subsampleg=1,vert=T)
dev.off()
}

#plot_freq_panel(1980,2016,seasons=rbind(c(5,10),c(11,4)),snames=c("MJJASO","NDJFMA"),
#        dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf",
#        type="anticyclone",reanal="ERAI",pnames=c("2 fixes","500km movement"),
#        proj=c("proj100_rad10cv0.075_D2","proj100_rad10cv0.075_500km"),latlim=c(20,60),lonlim=c(-30,50),
#        fout="paperfig_anticycfreq_ERAI_rad10cv0.075_threshcomp_medit.pdf")

#plot_freq_panel(1980,2016,seasons=cbind(c(12,3,6,9),c(2,5,8,11)),snames=c("DJF","MAM","JJA","SON"),
#        dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf",
#        type="anticyclone",reanal="ERAI",proj="proj100_rad10cv0.075_500km",
#        latlim=c(-50,-10),lonlim=c(90,180),breaks=c(0,0.05,seq(0.1,1,0.1),1000),
#        fout="paperfig_anticycfreq_ERAI_rad10cv0.075_4season_500km_australia.pdf")

#for(lev in c(500,700,850,925))
#{
#plot_freq_panel(1980,2016,seasons=cbind(c(5,11),c(10,4)),snames=c("May-October","November-April"),
#        dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf",
#        type="cyclone",reanal="ERAI",proj=paste0(lev,"hPa_proj100_rad5cv1.5_500km"),
#        latlim=c(-50,-10),lonlim=c(90,180),breaks=c(0,0.05,seq(0.1,1,0.1),1000),
#        fout=paste0("paperfig_cycfreq",lev,"hPa_ERAI_rad5cv1.5_2season_500km_australia.pdf"))
#}

#plot_freq_panel(1980,2016,seasons=cbind(c(5,11),c(10,4)),snames=c("May-October","November-April"),
#        dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf",
#        type="cyclone",reanal="ERAI",proj="proj100_rad5cv0.5_500km",
#        latlim=c(-50,-10),lonlim=c(90,180),breaks=c(0,0.05,seq(0.1,0.5,0.05),1000),
#        fout=paste0("paperfig_cycfreq_ERAI_rad5cv0.5_2season_500km_australia2.pdf"))

plot_freq_panel(1980,2016,seasons=cbind(c(5,11),c(10,4)),snames=c("May-October","November-April"),
        dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf",boxes=T,
        type="anticyclone",reanal="ERAI",proj="proj100_rad10cv0.075_500km",
        latlim=c(-50,-10),lonlim=c(90,180),breaks=c(0,0.05,seq(0.2,2,0.2),1000),
        fout="paperfig_anticycfreq_ERAI_rad10cv0.075_2season_500km_australia_boxes3.pdf")

