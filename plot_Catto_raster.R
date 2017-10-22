## This is code that takes the 11 ensemble members of a POAMA run
## Started on a particular day
## And outputs some interested results for further analysis

## I care about ECLs, so main result is a daily timeseries of # of members that
## have an ECL on that day

## Also outputs a figure for each of the months in the file
## That shows the number of cyclones centred in that month
## Option to set a higher intensity threshold
library(raster)
library(maps)
library(abind)
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

plot_counts<-function(year1,year2,month1=1,month2=12,fout=NA,Jdir="/g/data/eg3/ajd548/vicci/cyclone_data_JC/")
{
years=seq(year1,year2,1)

### Actually, use Jen's latitudes?
a=nc_open(paste(Jdir,"cyclones_1979_0.75.nc",sep=""))
lat=ncvar_get(a,"latitude")
lon=ncvar_get(a,"longitude")

if(month2>=month1) systems<-array(0,c(length(lon),length(lat),length(years))) else systems<-array(0,c(length(lon),length(lat),length(years)-1))


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
 systems[,,y]=apply(ncvar_get(a,"FLAG",start=c(1,1,I[1]),count=c(length(lon),length(lat),length(I))),c(1,2),sum)
} else {
 if(y==1)
 {
 I=which(month>=month1)
 store=apply(ncvar_get(a,"FLAG",start=c(1,1,I[1]),count=c(length(lon),length(lat),length(I))),c(1,2),sum)
 } else {
 I=which(month<=month2)
 systems[,,y-1]=store+apply(ncvar_get(a,"FLAG",start=c(1,1,I[1]),count=c(length(lon),length(lat),length(I))),c(1,2),sum)

 I=which(month>=month1)
 store=apply(ncvar_get(a,"FLAG",start=c(1,1,I[1]),count=c(length(lon),length(lat),length(I))),c(1,2),sum)
 }
}
}


print("Calculating means")
if(month2<month1) years=seq(year1,year2-1,1)

### Plot the average & the linear trend

### Mean frequency

meanfreq=apply(systems,c(1,2),mean,na.rm=T)

### Plot


if(month1==1 & month2==12) breaks=c(0,0.05,seq(0.5,2,0.5),3:5,1000) else breaks=c(0,0.05,seq(0.25,1,0.25),1.5,2,1000)
print("Plotting")
col1=col_val(length(breaks)-1)

if(is.na(fout)) fname=paste("Systems_mean_",year1,"_",year2,".pdf",sep="") else fname=paste(fout,".pdf",sep="")

pdf(file=fname,width=7,height=5.5)
layout(c(1,2),height=c(1,0.15))
par(mar=c(2,2,4,1))

image(lon,lat,meanfreq,breaks=breaks,col=col1,xlab="",ylab="",
          main=type[n])
map('world',add=T)
contour(raster(t(meanfreq[,length(lat):1]),xmn=-180,xmx=180,ymn=-89.5,ymx=89.5),levels=breaks,add=T,lwd=2,col="black",drawlabels=F)
ColorBar(breaks,col1,subsampleg=1,vert=F)
dev.off()
}


seas=c("","_MAM","_JJA","_SON","_DJF","_MJJASO","_NDJFMA")
m1=c(1,3,6,9,12,5,11)
m2=c(12,5,8,11,2,10,4)

for(s in c(1,6:7))
{

plot_counts(1980,2015,fout=paste("ERAI_Catto_meanfreq_1980-2015",seas[s],sep=""),month1=m1[s],month2=m2[s])
}

