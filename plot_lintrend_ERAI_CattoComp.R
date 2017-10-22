## This is code that takes the 11 ensemble members of a POAMA run
## Started on a particular day
## And outputs some interested results for further analysis

## I care about ECLs, so main result is a daily timeseries of # of members that
## have an ECL on that day

## Also outputs a figure for each of the months in the file
## That shows the number of cyclones centred in that month
## Option to set a higher intensity threshold

library(maps)
library(ncdf4)
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


filled.contour3 <-
  function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
            col = color.palette(length(levels) - 1), plot.title, plot.axes, 
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
            axes = TRUE, frame.plot = axes,mar, ...) 
  {
    # modification by Ian Taylor of the filled.contour function
    # to remove the key and facilitate overplotting with contour()
    # further modified by Carey McGilliard and Bridget Ferris
    # to allow multiple plots on one page
    
    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq.int(0, 1, length.out = nrow(z))
        }
      }
      else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
      stop("increasing 'x' and 'y' values expected")
    # mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    # on.exit(par(par.orig))
    # w <- (3 + mar.orig[2]) * par("csi") * 2.54
    # par(las = las)
    # mar <- mar.orig
    plot.new()
    # par(mar=mar)
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
      stop("no proper 'z' matrix specified")
    if (!is.double(z)) 
      storage.mode(z) <- "double"
    .filled.contour(as.double(x), as.double(y), z, as.double(levels), 
                            col = col)
    if (missing(plot.axes)) {
      if (axes) {
        title(main = "", xlab = "", ylab = "")
        Axis(x, side = 1)
        Axis(y, side = 2)
      }
    }
    else plot.axes
    if (frame.plot) 
      box()
    if (missing(plot.title)) 
      title(...)
    else plot.title
    invisible()
  }



plot_counts_many<-function(year1=1980,year2=2015,cv=NA,dur=NA,slp=NA,month1=1,month2=12,fout=NA,mydir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_lows_rad5cv0.15/",Jdir="/g/data/eg3/ajd548/vicci/cyclone_data_JC/")
{
years=seq(year1,year2,1)

### Actually, use Jen's latitudes?
a=nc_open(paste(Jdir,"cyclones_1979_0.75.nc",sep=""))
lat=ncvar_get(a,"latitude")
lon=ncvar_get(a,"longitude")

## For interval stuff, need a different lat/lon that is inclusive
lat2=c(lat-((lat[2]-lat[1])/2),lat[length(lat)]+((lat[2]-lat[1])/2))
lon2=c(lon-((lon[2]-lon[1])/2),lon[length(lon)]+((lon[2]-lon[1])/2))


if(month2>=month1) systems<-array(0,c(length(lon),length(lat),length(years),2)) else systems<-array(0,c(length(lon),length(lat),length(years)-1,2))

### Load my version like normal

for(y in 1:length(years))
{
fname=paste(mydir,"/tracks_",years[y],".dat",sep="")
read.table(fname, sep="",skip=1)->fixes
colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Depth","Radius","Up","Vp")
#fixes=fixes[fixes$Lat<0,]
fixes$Year=floor(fixes$Date/10000)
fixes$CV=abs(fixes$CV)
fixes$Depth=abs(fixes$Depth)
if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
fixes$Month=floor(fixes$Date/100)%%100
fixes$Lon=fixes$Lon%%360
fixes$Lon[fixes$Lon>180]=fixes$Lon[fixes$Lon>180]-360

fixes$Lat2=findInterval(fixes$Lat,lat2,all.inside=TRUE)
fixes$Lon2=findInterval(fixes$Lon,lon2,all.inside=TRUE)

 if(!is.na(dur))
  {
  x<-rle(fixes$ID)
  events<-cbind(x$values,x$lengths,matrix(data=0,nrow=length(x$values),ncol=1))
  events=events[events[,2]>=dur,]
  include<-match(fixes[,1],events[,1])
  J<-which(is.na(include)==0)
  fixes=fixes[J,]
  }
if(!is.na(cv)) fixes=fixes[fixes$CV>=cv,]
if(!is.na(slp)) fixes=fixes[fixes$MSLP>=slp,]

if(month2>=month1)
 {
  I=which(fixes$Month>=month1 & fixes$Month<=month2)
  systems[,,y,1]=table(factor(fixes$Lon2[I],levels=1:length(lon)),factor(fixes$Lat2[I],levels=1:length(lat)))

 } else {
  if(y==1)
  {
   I=which(fixes$Month>=month1)
   store=fixes[I,]
  } else {
   I=which(fixes$Month<=month2)
   fixes2=rbind(store,fixes[I,])
   systems[,,y-1,1]=table(factor(fixes2$Lon2[I],levels=1:length(lon)),factor(fixes2$Lat2[I],levels=1:length(lat)))

   I=which(fixes$Month>=month1)
   store=fixes[I,]
  }
}
} # End year loop

### Make sum within +- 5 degrees, only way to even remotely compare

systems2=systems

for(i in 1:length(lon))
  for(j in 1:length(lat))
  {
     I=which(abs(lon-lon[i])<5 | abs(lon-lon[i])>355)
     J=which(abs(lat-lat[j])<5)
     systems[i,j,,1]=apply(systems2[I,J,,1],3,sum,na.rm=T)
  }

### Now, Jen's data

for(y in 1:length(years))
{
## First make list of dates
a=nc_open(paste(Jdir,"cyclones_",years[y],"_0.75.nc",sep=""))

from=as.POSIXlt(paste(years[y],"-01-01 00:00",sep=""),format="%Y-%m-%d %H:%M",tz="GMT")
to=as.POSIXlt(paste(years[y],"-12-31 18:00",sep=""),format="%Y-%m-%d %H:%M",tz="GMT")

time=seq.POSIXt(from=from,to=to,by="6 hours")
month=as.numeric(format(time,"%m"))

if(month2>=month1)
{
 I=which(month>=month1 & month<=month2)
 systems[,,y,2]=apply(ncvar_get(a,"FLAG",start=c(1,1,I[1]),count=c(length(lon),length(lat),length(I))),c(1,2),sum)
} else {
 if(y==1)
 {
 I=which(month>=month1)
 store=apply(ncvar_get(a,"FLAG",start=c(1,1,I[1]),count=c(length(lon),length(lat),length(I))),c(1,2),sum)
 } else {
 I=which(month<=month2)
 systems[,,y-1,2]=store+apply(ncvar_get(a,"FLAG",start=c(1,1,I[1]),count=c(length(lon),length(lat),length(I))),c(1,2),sum)

 I=which(month>=month1)
 store=apply(ncvar_get(a,"FLAG",start=c(1,1,I[1]),count=c(length(lon),length(lat),length(I))),c(1,2),sum)
 }
}
}

if(month2<month1) years=seq(year1,year2-1,1)

### Linear trend
meanfreq=apply(systems,c(1,2,4),mean)
cyctrend<-array(NaN,c(length(lon),length(lat),2,2))

for(i in 1:length(lon))
  for(j in 1:length(lat))
   for(em in 1:2)
     if(meanfreq[i,j,em]>=5)
     {
      a=lm(systems[i,j,,em]~years)
      b=summary(a)$coefficients
      cyctrend[i,j,em,1]=100*a$coefficients[2]/meanfreq[i,j,em] 
      cyctrend[i,j,em,2]=b[2,4]          
     }

### Plot


breaks1=c(0,1,2,5,10,20,50,100,200,500,1000,10000)
breaks2=c(-100,seq(-3,3,0.5),100)

print("Plotting")
col1=col_val(length(breaks1)-1)
col2=col_anom(length(breaks2)-1)

if(is.na(fout)) fname=paste("ERAI_trendcomp_",year1,"_",year2,"_global.pdf",sep="") else fname=fout

pdf(file=fname,width=15,height=9)
layout(cbind(c(1,2),c(3,4),c(5,6)),width=c(1,1,0.15))

par(mar=c(2,2,4,1))

names=c("UM cyclone centres w/in 5 degrees","Jen Catto cyclone area")

for(s in 1:2)
{
image(lon,lat,meanfreq[,,s],breaks=breaks1,col=col1,xlab="",ylab="",
          main=paste(names[s],"Cyclones:",year1,"-",year2))
#filled.contour3(lon,lat,meanfreq[,,s],levels=breaks1,col=col1,xlab="",ylab="",
#          main=names[s])
contour(lon,lat,meanfreq[,,s],levels=breaks1[seq(3,length(breaks1),by=2)],add=T,lwd=2,col="black",drawlabels=F)
map('world',add=T)

image(lon,lat,cyctrend[,,s,1],breaks=breaks2,col=col2,xlab="",ylab="",
          main=paste("Linear trend:",names[s]))
map('world',add=T)
#sigmask=which(cyctrend[,,s,2]<0.05 & meanfreq[,,s]>=5,arr.ind=T)
#points(lon[sigmask[,1]],lat[sigmask[,2]],col="black",pch=19,cex=0.6)

contour(lon,lat,cyctrend[,,s,2],levels=c(-100,0.05,100),add=T,lwd=2,col="black",drawlabels=F)

}

ColorBar(breaks1,col1,subsampleg=1)
ColorBar(breaks2,col2,subsampleg=1)
dev.off()

}

seas=c("","_MAM","_JJA","_SON","_DJF","_MJJ","_MJJASO","_NDJFMA")
m1=c(1,3,6,9,12,5,5,11)
m2=c(12,5,8,11,2,7,10,4)

for(s in 1:8)
{
  
plot_counts_many(1980,2015,fout=paste("ERAI_UM_cyclonetrend_comp_1980-2016_rad5cv15_D2",seas[s],".pdf",sep=""),month1=m1[s],month2=m2[s],dur=2)
plot_counts_many(1980,2015,fout=paste("ERAI_UM_cyclonetrend_comp_1980-2016_rad5cv25_D2",seas[s],".pdf",sep=""),month1=m1[s],month2=m2[s],dur=2,cv=0.25)
}
