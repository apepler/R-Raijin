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
  tmp=t(as.matrix(m3)[a[2]:1,(winwid+1):(a[1]+winwid)])
  return(tmp)
}


plot_freq_panel<-function(season=c(1,12),sname="Annual",
       leads=cbind(c(1,7),c(1,30)),lnames=c("Days 1-7","Days 0-30"),
       lonlim=c(0,360),latlim=c(-90,90),
       dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/",
       efile=NaN,afile=NaN,
       fout="output")
{

lat=seq(-89.5,89.5)
lon=seq(0,359.5)  ### Can always combine into bigger cells later
mdays=c(31,28.25,31,30,31,30,31,31,30,31,30,31)
mapname="world2"
if(season[2]>=season[1]) mlist=season[1]:season[2] else mlist=c(season[2]:12,1:season[1])

breaks1=c(0,0.05,seq(0.2,2,0.2),1000)
col1=col_val(length(breaks1)-1)
breaks2=c(-1000,seq(-100,100,20),1000)
col2=col_anom(length(breaks2)-1)

pdf(file=paste0(fout,".pdf"),width=13,height=6)
layout(cbind(c(1,8),c(2,3),c(4,5),c(6,7)),width=c(1,1,1,0.25))

par(mar=c(2,2,4,1))

# First, get files

a=nc_open(paste0(dir,efile)) ## All, including open
erai=ncvar_get(a,"systems")

a=nc_open(paste0(dir,afile)) ## All, including open
access=ncvar_get(a,"cyclones")

## Reshape if necessary
if(lonlim[1]<0) {
   lon=seq(-179.5,179.5,1)
   erai=abind(erai[181:360,,,],erai[1:180,,,],along=1)
   access=abind(access[181:360,,,,],access[1:180,,,,],along=1)
   mapname="world"
   }

# First, plot ERAI frequency

tmp=apply(erai[,,1:23,mlist],c(1,2),sum)
densE=makesmooth(tmp,5)*1000/(23*sum(mdays[mlist]))

image(lon,lat,densE,breaks=breaks1,col=col1,
      xlab="",ylab="",xlim=lonlim,ylim=latlim,
      main=paste0("a) ERAI system density: ",sname," (10^-3)"))
map(mapname,add=T)  
contour(lon,lat,densE,levels=breaks1[seq(3,length(breaks1),2)],add=T,lwd=2,col="black",drawlabels=F)


for(n in 1:length(lnames))
{
  llist=seq(leads[1,n],leads[2,n])
  tmp=apply(access[,,,mlist,llist],c(1,2),sum) ## Density - likelihood of cyclone
  densA=makesmooth(tmp,5)*1000/(23*11*length(mlist)*length(llist))

  image(lon,lat,densA,breaks=breaks1,col=col1,
      xlab="",ylab="",xlim=lonlim,ylim=latlim,
      main=paste0(letters[n*2],") ACCESS-S system density: ",sname," ",lnames[n]," (10^-3)"))
  map(mapname,add=T) 
  contour(lon,lat,densA,levels=breaks1[seq(3,length(breaks1),2)],add=T,lwd=2,col="black",drawlabels=F)

  percdiff=100*((densA/densE)-1)
  percdiff[densE<0.05 | densA<0.05]=NaN

  image(lon,lat,percdiff,breaks=breaks2,col=col2,
      xlab="",ylab="",xlim=lonlim,ylim=latlim,
      main=paste0(letters[n*2+1],") ACCESS-S bias: ",sname," ",lnames[n]," (%)"))
  map(mapname,add=T)
  contour(lon,lat,percdiff,levels=c(-80,-40,0,40,80),add=T,lwd=2,col="black",drawlabels=F)

}

ColorBar(breaks1,col1,subsampleg=1,vert=T)
ColorBar(breaks2,col2,subsampleg=1,vert=T)
dev.off()
}

snames=c("Annual","DJF","MAM","JJA","SON","MJJASO","NDJFMA")
seasons=cbind(c(1,12,3,6,9,5,11),c(12,2,5,8,11,10,4))

for(s in 1:length(snames))
{
plot_freq_panel(season=seasons[s,],sname=snames[s],
        dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/",
        efile="ERAIdaily_UM_globalcyclones_proj240_rad5cv0.25.nc",
        afile="ACCESS_globalcyclones_proj240_rad5cv0.25_40daylead.nc",
        fout=paste0("cycfreq_ACCESSvERAI_proj240_rad5cv0.25_global_",snames[s]))
}
