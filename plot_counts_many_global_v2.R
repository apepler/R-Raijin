## This is code that takes the 11 ensemble members of a POAMA run
## Started on a particular day
## And outputs some interested results for further analysis

## I care about ECLs, so main result is a daily timeseries of # of members that
## have an ECL on that day

## Also outputs a figure for each of the months in the file
## That shows the number of cyclones centred in that month
## Option to set a higher intensity threshold

library(maps)
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

plot_counts_many<-function(year1,year2,dir="gcyc_out",name="proj100_highs_rad10cv0.075",cv=NA,dur=NA,slp=NA,month1=1,month2=12,fout=NA,hilo="low",members=56)
{
years=seq(year1,year2,1)

##Only include ERAI if year1>=1980
if(year1>=1980)
{
sdirs=c("20CR/","20CR/EnsMean/","NCEP1/","ERAI/")
snames=c("20CRensemble","20CRmean","NCEP1","ERAI")
} else {
sdirs=c("20CR/","20CR/EnsMean/","NCEP1/")
snames=c("20CRensemble","20CRmean","NCEP1")
}

lat=seq(-89.5,89.5)
lon=seq(0.5,359.5)  ### Can always combine into bigger cells later
if(month2>=month1) systems<-array(0,c(length(lon),length(lat),length(years),members,4)) else systems<-array(0,c(length(lon),length(lat),length(years)-1,members,4))
dimnames(systems)[[5]]=c("Freq","MSLP","CV","LonMove")

systems2=systems[,,,1:length(snames),]

for(y in 1:length(years))
  for(em in 1:members)
{
print(paste(years[y],em))
fname=paste(dir,sdirs[1],name,"/tracks_",years[y],"_",em,".dat",sep="") 
read.table(fname, sep="",skip=1)->fixes
colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Depth","Radius","Up","Vp")
#fixes=fixes[fixes$Lat<0,]
fixes$Year=floor(fixes$Date/10000)
fixes$CV=abs(fixes$CV)
fixes$Depth=abs(fixes$Depth)
if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
fixes$Month=floor(fixes$Date/100)%%100
fixes$Lat2=floor(fixes$Lat) ## Let's skip the first step, and just make it the 10 degrees
fixes$Lon2=floor((fixes$Lon%%360))

  fixes$LonMove<-NaN
  I=which(fixes$Fix!=1)
  if(I[1]==1) I=I[-1]
  fixes$LonMove[I]=(fixes$Lon[I]-fixes$Lon[I-1])%%360

### Make table of events to combine with DJF for exclusion
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
  systems[,,y,em,1]=table(factor(fixes$Lon2[I],levels=seq(0,359)),factor(fixes$Lat2[I],levels=seq(-90,89)))

  tmp2=aggregate(fixes[I,8:9],by=list(fixes$Lon2[I],fixes$Lat2[I]),FUN=mean)
  jlat=match(tmp2[,2],seq(-90,89))
  ilon=match(tmp2[,1],seq(0,359))
  for(nn in which(!is.na(jlat))) systems[ilon[nn],jlat[nn],y,em,2:3]=as.numeric(tmp2[nn,3:4])

  I=which(fixes$Month>=month1 & fixes$Month<=month2 & !is.na(fixes$LonMove))
  tmp2=aggregate(fixes$LonMove[I],by=list(fixes$Lon2[I],fixes$Lat2[I]),FUN=mean)
  jlat=match(tmp2[,2],seq(-90,89))
  ilon=match(tmp2[,1],seq(0,359))
  for(nn in which(!is.na(jlat))) systems[ilon[nn],jlat[nn],y,em,4]=as.numeric(tmp2[nn,3])


 } else {
  if(y==1)
  {
   I=which(fixes$Month>=month1)
   store=fixes[I,]
  } else {
   I=which(fixes$Month<=month2)
   fixes2=rbind(store,fixes[I,])

   systems[,,y-1,em,1]=table(factor(fixes2$Lon2,levels=seq(0,359)),factor(fixes2$Lat2,levels=seq(-90,89)))

tmp2=aggregate(fixes2[,8:9],by=list(fixes2$Lon2,fixes2$Lat2),FUN=mean)
   jlat=match(tmp2[,2],seq(-90,89))
   ilon=match(tmp2[,1],seq(0,359))
   for(nn in which(!is.na(jlat))) systems[ilon[nn],jlat[nn],y-1,em,2:3]=as.numeric(tmp2[nn,3:4])

  I=which(!is.na(fixes2$LonMove))
  tmp2=aggregate(fixes2$LonMove[I],by=list(fixes2$Lon2[I],fixes2$Lat2[I]),FUN=mean)
  jlat=match(tmp2[,2],seq(-90,89))
  ilon=match(tmp2[,1],seq(0,359))
  for(nn in which(!is.na(jlat))) systems[ilon[nn],jlat[nn],y-1,em,4]=as.numeric(tmp2[nn,3])

   I=which(fixes$Month>=month1)
   store=fixes[I,]
  }
}
} # End year loop

print("Calculating means")
#if(month2<month1) years=seq(year1,year2-1,1)

### Plot the average & the linear trend

### Now, add the other ones

systems2[,,,1,]=apply(systems[,,,1:members,],c(1,2,3,5),mean)

for(y in 1:length(years))
for(em in 2:length(sdirs))
{
if(em==3) fname=paste(dir,sdirs[em],name,"_v2/tracks_",years[y],".dat",sep="") else fname=paste(dir,sdirs[em],name,"/tracks_",years[y],".dat",sep="")
read.table(fname, sep="",skip=1)->fixes
colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Depth","Radius","Up","Vp")
#fixes=fixes[fixes$Lat<0,]
fixes$Year=floor(fixes$Date/10000)
fixes$CV=abs(fixes$CV)
fixes$Depth=abs(fixes$Depth)
if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
fixes$Month=floor(fixes$Date/100)%%100
fixes$Lat2=floor(fixes$Lat) ## Let's skip the first step, and just make it the 10 degrees
fixes$Lon2=floor((fixes$Lon%%360))

  fixes$LonMove<-NaN
  I=which(fixes$Fix!=1)
  if(I[1]==1) I=I[-1]
  fixes$LonMove[I]=(fixes$Lon[I]-fixes$Lon[I-1])%%360


### Make table of events to combine with DJF for exclusion
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
  systems2[,,y,em,1]=table(factor(fixes$Lon2[I],levels=seq(0,359)),factor(fixes$Lat2[I],levels=seq(-90,89)))

  tmp2=aggregate(fixes[I,8:9],by=list(fixes$Lon2[I],fixes$Lat2[I]),FUN=mean)
  jlat=match(tmp2[,2],seq(-90,89))
  ilon=match(tmp2[,1],seq(0,359))
  for(nn in which(!is.na(jlat))) systems2[ilon[nn],jlat[nn],y,em,2:3]=as.numeric(tmp2[nn,3:4])

  I=which(fixes$Month>=month1 & fixes$Month<=month2 & !is.na(fixes$LonMove))
  tmp2=aggregate(fixes$LonMove[I],by=list(fixes$Lon2[I],fixes$Lat2[I]),FUN=mean)
  jlat=match(tmp2[,2],seq(-90,89))
  ilon=match(tmp2[,1],seq(0,359))
  for(nn in which(!is.na(jlat))) systems2[ilon[nn],jlat[nn],y,em,4]=as.numeric(tmp2[nn,3])


 } else {
  if(y==1)
  {
   I=which(fixes$Month>=month1)
   store=fixes[I,]
  } else {
   I=which(fixes$Month<=month2)
   fixes2=rbind(store,fixes[I,])

   systems2[,,y-1,em,1]=table(factor(fixes2$Lon2,levels=seq(0,359)),factor(fixes2$Lat2,levels=seq(-90,89)))
   tmp2=aggregate(fixes2[,8:9],by=list(fixes2$Lon2,fixes2$Lat2),FUN=mean)
   jlat=match(tmp2[,2],seq(-90,89))
   ilon=match(tmp2[,1],seq(0,359))
   for(nn in which(!is.na(jlat))) systems2[ilon[nn],jlat[nn],y-1,em,2:3]=as.numeric(tmp2[nn,3:4])

  I=which(!is.na(fixes2$LonMove))
  tmp2=aggregate(fixes2$LonMove[I],by=list(fixes2$Lon2[I],fixes2$Lat2[I]),FUN=mean)
  jlat=match(tmp2[,2],seq(-90,89))
  ilon=match(tmp2[,1],seq(0,359))
  for(nn in which(!is.na(jlat))) systems2[ilon[nn],jlat[nn],y-1,em,4]=as.numeric(tmp2[nn,3])

   I=which(fixes$Month>=month1)
   store=fixes[I,]
  }
}
} # End year loop

if(month2<month1) years=seq(year1,year2-1,1)

### Now to do some averaging across the region

systems=systems2

for(i in 1:length(lon))
  for(j in 1:length(lat))
  {
     I=which((lon>lon[i]-5 & lon<lon[i]+5) | lon>lon[i]+355 | lon<lon[i]-355)
     J=which(lat>lat[j]-5 & lat<lat[j]+5)
     systems[i,j,,,1]=apply(systems2[I,J,,,1],c(3,4),sum,na.rm=T)
     for(n in 2:4) systems[i,j,,,n]=apply(systems2[I,J,,,1]*systems2[I,J,,,n],c(3,4),sum,na.rm=T)/apply(systems2[I,J,,,1],c(3,4),sum,na.rm=T)
  }

systems2=systems

### Linear trend
meanfreq=apply(systems2,c(1,2,4,5),mean)
cyctrend<-array(NaN,c(length(lon),length(lat),length(snames),4,2))

for(i in 1:length(lon))
  for(j in 1:length(lat))
   for(em in 1:length(snames))
    for(k in 1:4)
     if(!is.na(meanfreq[i,j,em,1]))
     if(meanfreq[i,j,em,1]>=5)
     {
      a=lm(systems2[i,j,,em,k]~years)
      b=summary(a)$coefficients
      if(k==1) cyctrend[i,j,em,k,1]=100*a$coefficients[2]/meanfreq[i,j,em,k] else cyctrend[i,j,em,k,1]=10*a$coefficients[2]
      cyctrend[i,j,em,k,2]=b[2,4]          
     }

### Plot

type=c("freq","MSLP","CV","LonMove")
units1=c("(Systems/deg^2)","(hPa)","(hPa/(deg.lat)^2)","(degrees/6h)")
units2=c("(%/year)","(hPa/decade)","(Laplacian/decade)","(degrees/decade)")

if(hilo=="low")
{
blist1=list(c(0,0.05,seq(0.25,1,0.25),2:5,1000),
            c(900,seq(970,1020,5),1100),
            c(0,seq(0.15,0.75,0.05),100),
            c(0:10,100),
            c(0,seq(3,7,0.5),100))
} else {
blist1=list(c(0,0.05,seq(0.25,1,0.25),2:5,1000),
            c(900,seq(1010,1050,5),10000),
            c(0,0.075,seq(0.1,0.4,0.05),100),
            c(0:10,100),
            c(0,seq(3,7,0.5),100))
}
blist2=list(c(-100,seq(-3,3,0.5),100),
            c(-100,seq(-3,3,0.5),100),
            c(-100,seq(-0.1,0.1,0.02),100),
            c(-100,seq(-1,1,0.2),100),
            c(-100,seq(-0.25,0.25,0.05),100))

for(k in 1:4)
{
breaks1=blist1[[k]]
breaks2=blist2[[k]]

print("Plotting")
col1=col_val(length(breaks1)-1)
col2=col_anom(length(breaks2)-1)

if(k==2) col1=col1[length(col1):1]

if(is.na(fout)) fname=paste("Systems_",type[k],"_trend_",year1,"_",year2,"_global.pdf",sep="") else fname=paste(fout,"_all_",type[k],"_global2.pdf",sep="")


pdf(file=fname,width=14,height=10)
layout(cbind(c(1,2,5),c(3,4,5)),height=c(1,1,0.15))

par(mar=c(2,2,4,1))

for(s in 1:length(sdirs))
{
image(lon,lat,cyctrend[,,s,k,1],breaks=breaks2,col=col2,xlab="",ylab="",
          main=paste("Trend in",snames[s],":",year1,"-",year2,units2[k]))
map('world2',add=T)
contour(lon,lat,cyctrend[,,s,k,2],levels=c(-100,0.05,100),add=T,lwd=2,col="black",drawlabels=F)
}
if(year1<1980) plot.new()
ColorBar(breaks2,col2,subsampleg=1,vert=F)
dev.off()

if(is.na(fout)) fname=paste("Systems_",type[k],"_trend_",year1,"_",year2,"_global.pdf",sep="") else fname=paste(fout,"_mean_",type[k],"_global2.pdf",sep="")
pdf(file=fname,width=8,height=5)
layout(cbind(1,2),width=c(1,0.13))
par(mar=c(3,3,4,1))
image(lon,lat,apply(cyctrend[,,,k,1],c(1,2),mean),breaks=breaks2,col=col2,xlab="",ylab="",
          main=paste("Mean trend:",year1,"-",year2,units2[k]))
map('world2',add=T)
contour(lon,lat,apply(cyctrend[,,,k,2]<0.05,c(1,2),sum),levels=c(-100,length(sdirs)-1,100),add=T,lwd=2,col="black",drawlabels=F)
ColorBar(breaks2,col2,subsampleg=1)
dev.off()
}
}

seas=c("","_MAM","_JJA","_SON","_DJF")
m1=c(1,3,6,9,12)
m2=c(12,5,8,11,2)

for(s in 1)
{
  
plot_counts_many(1950,2014,dir="/short/eg3/asp561/cts.dir/gcyc_out/",name="proj100_highs_rad10cv0.075",fout=paste("UM_anticyclonetrend_1950-2014_rad10cv075_D2",seas[s],sep=""),hilo="high",month1=m1[s],month2=m2[s],dur=2)
  plot_counts_many(1950,2014,dir="/short/eg3/asp561/cts.dir/gcyc_out/",name="proj100_lows_rad5cv0.15",fout=paste("UM_cyclonetrend_1950-2014_rad5cv15_D2",seas[s],sep=""),hilo="low",month1=m1[s],month2=m2[s],dur=2)  



}
