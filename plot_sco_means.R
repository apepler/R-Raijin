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

years=1990:2012

lat=seq(-89.5,89.5)
lon=seq(0,359.5)  ### Can always combine into bigger cells later
systems<-array(0,c(length(lon),length(lat),length(years),12,3)) 
dimnames(systems)[[5]]<-type<-c("ERAI","POAMA","ACCESS")
setwd("/short/eg3/asp561/cts.dir/gcyc_out/")

cv=0.075
proj="proj240_highs_rad10cv0.075"
### First, ERAI

for(y in 1:length(years))
{
print(years[y])
fname=paste("ERAI/daily_",proj,"/tracks_",years[y],".dat",sep="")
read.table(fname, sep="",skip=1)->fixes
colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Depth","Radius","Up","Vp")
fixes$Year=floor(fixes$Date/10000)
if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
fixes$Month=floor(fixes$Date/100)%%100
fixes$Lat2=floor(fixes$Lat) ## Let's skip the first step, and just make it the 10 degrees
fixes$Lon2=floor((fixes$Lon%%360))
### Make table of events to combine with DJF for exclusion
fixes$CV=abs(fixes$CV)

fixes=fixes[fixes$CV>=cv,]
systems[,,y,,1]=table(factor(fixes$Lon2,levels=seq(0,359)),factor(fixes$Lat2,levels=seq(-90,89)),factor(fixes$Month,levels=1:12))
} # End year loop

fname=paste("POAMA/",proj,"/e24a/globalanticyclones_4monthlead_day1_cv",cv,".nc",sep="")
a=nc_open(fname)
ylist=ncvar_get(a,"year")
members=ncvar_get(a,"members")
Y=which(ylist%in%years)
tmp=ncvar_get(a,"cyclones",start=c(1,1,Y[1],1,1,1),count=c(length(lon),length(lat),length(Y),12,length(members),1)) ## Only leadtime 1 
systems[,,,,2]=apply(tmp,c(1,2,3,4),mean)
rm(tmp)

fname=paste("access-s1/",proj,"/globalanticyclones_4monthlead_day1_cv",cv,".nc",sep="")
a=nc_open(fname)
ylist=ncvar_get(a,"year")
members=ncvar_get(a,"members")
Y=which(ylist%in%years)
tmp=ncvar_get(a,"cyclones",start=c(1,1,Y[1],1,1,1),count=c(length(lon),length(lat),length(Y),12,length(members),1)) ## Only leadtime 1 
systems[,,,,3]=apply(tmp,c(1,2,3,4),mean)

## Now make the +- 5 degrees version
systems2=systems

for(i in 1:length(lon))
  for(j in 1:length(lat))
  {
     I=which((lon>lon[i]-5 & lon<lon[i]+5) | lon>lon[i]+355 | lon<lon[i]-355)
     J=which(lat>lat[j]-5 & lat<lat[j]+5)
     systems[i,j,,,]=apply(systems2[I,J,,,],c(3,4,5),sum,na.rm=T)
  }


meanfreq=apply(systems,c(1,2,4,5),mean,na.rm=T)
rm(tmp)

### Plot
breaks=c(0,0.1,1,2,4,6,8,10,15,20,10000)
print("Plotting")
col1=col_val(length(breaks)-1)

for(m in 1:12)
{

fname=paste("Anticycfreq_ACCESSvPOAMA_",sprintf("%2.2d",m),"01_lead0_cv",cv,".pdf",sep="")

pdf(file=fname,width=14,height=10)
layout(cbind(c(1,2,4),c(5,3,4)),height=c(1,1,0.15)) ## ERAI on top, both models below
par(mar=c(2,2,4,1))

for(n in 1:3)
{
image(lon,lat,meanfreq[,,m,n],breaks=breaks,col=col1,xlab="",ylab="",
          main=type[n])
map('world2',add=T)
contour(lon,lat,meanfreq[,,m,n],levels=breaks,add=T,lwd=2,col="black",drawlabels=F)
}

ColorBar(breaks,col1,subsampleg=1,vert=F)
dev.off()
}

meanfreq=apply(apply(systems,c(1,2,3,5),sum),c(1,2,4),mean,na.rm=T)

### Plot
breaks=c(0,1,5,10,25,50,75,100,10000)
print("Plotting")
col1=col_val(length(breaks)-1)

fname=paste("Anticycfreq_ACCESSvPOAMA_annual_lead0_cv",cv,".pdf",sep="")

pdf(file=fname,width=14,height=10)
layout(cbind(c(1,2,4),c(5,3,4)),height=c(1,1,0.15)) ## ERAI on top, both models below
par(mar=c(2,2,4,1))

for(n in 1:3)
{
image(lon,lat,meanfreq[,,n],breaks=breaks,col=col1,xlab="",ylab="",
          main=type[n])
map('world2',add=T)
contour(lon,lat,meanfreq[,,n],levels=breaks,add=T,lwd=2,col="black",drawlabels=F)
}
ColorBar(breaks,col1,subsampleg=1,vert=F)
dev.off()


