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
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(1.5, length(brks) - 1.5, subsampleg), tick = TRUE, 
       labels = brks[seq(2, length(brks)-1, subsampleg)])
}

plot_counts<-function(year1,year2,dir="gcyc_out",cv=NA,dur=NA,slp=NA,month1=1,month2=12,fout=NA,hilo="low")
{
years=seq(year1,year2,1)

lat=seq(-85,-5,10)
lon=seq(5,355,10)  ### Can always combine into bigger cells later
if(month2>=month1) systems<-array(0,c(length(lon),length(lat),length(years),5)) else systems<-array(0,c(length(lon),length(lat),length(years)-1,5))
dimnames(systems)[[4]]=c("Freq","MSLP","CV","Depth","Radius")

for(y in 1:length(years))
{
print(years[y])
fname=paste(dir,"/tracks_",years[y],".dat",sep="")
read.table(fname, sep="",skip=1)->fixes
colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Depth","Radius","Up","Vp")
fixes=fixes[fixes$Lat<0,]
fixes$Year=floor(fixes$Date/10000)
fixes$CV=abs(fixes$CV)
fixes$Depth=abs(fixes$Depth)
if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
fixes$Month=floor(fixes$Date/100)%%100
fixes$Lat2=10*floor(fixes$Lat/10) ## Let's skip the first step, and just make it the 10 degrees
fixes$Lon2=10*floor((fixes$Lon%%360)/10)
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
  systems[,,y,1]=table(factor(fixes$Lon2[I],levels=seq(0,350,10)),factor(fixes$Lat2[I],levels=seq(-90,-10,10)))

  tmp2=aggregate(fixes[I,8:11],by=list(fixes$Lon2[I],fixes$Lat2[I]),FUN=mean)
  jlat=match(tmp2[,2],seq(-90,-10,10))
  ilon=match(tmp2[,1],seq(0,350,10))
  for(nn in 1:length(tmp2[,1])) systems[ilon[nn],jlat[nn],y,2:5]=as.numeric(tmp2[nn,3:6])
 } else {
  if(y==1)
  {
   I=which(fixes$Month>=month1)
   store=fixes[I,]
  } else {
   I=which(fixes$Month<=month2)
   fixes2=rbind(store,fixes[I,])

   systems[,,y-1,1]=table(factor(fixes2$Lon2,levels=seq(0,350,10)),factor(fixes2$Lat2,levels=seq(-90,-10,10)))
   tmp2=aggregate(fixes2[,8:11],by=list(fixes2$Lon2,fixes2$Lat2),FUN=mean)
   jlat=match(tmp2[,2],seq(-90,-10,10))
   ilon=match(tmp2[,1],seq(0,350,10))
   for(nn in 1:length(tmp2[,1])) systems[ilon[nn],jlat[nn],y-1,2:5]=as.numeric(tmp2[nn,3:6])

   I=which(fixes$Month>=month1)
   store=fixes[I,]
  }
}
} # End year loop

print("Calculating means")
if(month2<month1) years=seq(year1,year2-1,1)

### Plot the average & the linear trend

### Mean frequency

meanfreq=apply(systems,c(1,2,4),mean)

### Linear trend

cyctrend<-array(NaN,c(length(lon),length(lat),5,2))

for(i in 1:length(lon))
  for(j in 1:length(lat))
    for(k in 1:5)
     if(meanfreq[i,j,1]>=5)
     {
      a=lm(systems[i,j,,k]~years)
      b=summary(a)$coefficients
      if(k==1) cyctrend[i,j,k,1]=100*a$coefficients[2]/meanfreq[i,j,k] else cyctrend[i,j,k,1]=10*a$coefficients[2]
      cyctrend[i,j,k,2]=b[2,4]          
     }

### Plot

type=c("freq","MSLP","CV","Depth","Radius")
units1=c("Cyclones/deg^2)","(hPa)","(hPa/(deg.lat)^2)","hPa","(degrees)")
units2=c("(%/year)","(hPa/decade)","(Laplacian/decade)","(hPa/decade)","(degrees/decade)")

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

for(k in 1:5)
{
breaks1=blist1[[k]]
breaks2=blist2[[k]]

print("Plotting")
col1=col_val(length(breaks1)-1)
col2=col_anom(length(breaks2)-1)

if(k==2) col1=col1[length(col1):1]

if(is.na(fout)) fname=paste("Systems_",type[k],"_trend_",year1,"_",year2,".pdf",sep="") else fname=paste(fout,"_",type[k],".pdf",sep="")

pdf(file=fname,width=10,height=7)
layout(cbind(c(1,2),c(3,4)),width=c(1,0.15))
par(mar=c(3,3,4,1))

if(type[k]=="freq") tmp=meanfreq[,,k]/100 else tmp=meanfreq[,,k] # per degree
image(lon,lat,tmp,breaks=breaks1,col=col1,xlab="",ylab="",
          main=paste("Mean:",year1,"-",year2,units1[k]))
map('world2',add=T)

image(lon,lat,cyctrend[,,k,1],breaks=breaks2,col=col2,xlab="",ylab="",
          main=paste("Linear trend:",year1,"-",year2,units2[k]))
map('world2',add=T)
sigmask=which(cyctrend[,,k,2]<0.05 & meanfreq[,,1]>=5,arr.ind=T)
points(lon[sigmask[,1]],lat[sigmask[,2]],col="black",pch=19,cex=0.6)

ColorBar(breaks1,col1,subsampleg=1)
ColorBar(breaks2,col2,subsampleg=1)
dev.off()
}
}

#plot_counts(1950,2016,dir="/short/eg3/asp561/cts.dir/gcyc_out/NCEP1/proj100_rad5cv0.15",fout="NCEP1_UM_cyclonetrend_1950-2016_rad5cv25_D2_v2_MAM",hilo="low",month1=3,month2=5,cv=0.25,dur=2)

#plot_counts(1950,2016,dir="/short/eg3/asp561/cts.dir/gcyc_out/NCEP1/proj100_highs_rad10cv0.075",fout="NCEP1_UM_anticyclonetrend_1950-2016_rad10cv075_slp1020_D2_v2_MAM",hilo="high",month1=3,month2=5,slp=1020,dur=2)


#plot_counts(1950,2016,dir="/short/eg3/asp561/cts.dir/gcyc_out/NCEP1/proj100_rad5cv0.15",fout="NCEP1_UM_cyclonetrend_1950-2016_rad5cv25_D2_v2_JJA",hilo="low",month1=6,month2=8,cv=0.25,dur=2)

#plot_counts(1950,2016,dir="/short/eg3/asp561/cts.dir/gcyc_out/NCEP1/proj100_highs_rad10cv0.075",fout="NCEP1_UM_anticyclonetrend_1950-2016_rad10cv075_slp1020_D2_v2_JJA",hilo="high",month1=6,month2=8,slp=1020,dur=2)


#plot_counts(1950,2016,dir="/short/eg3/asp561/cts.dir/gcyc_out/NCEP1/proj100_rad5cv0.15",fout="NCEP1_UM_cyclonetrend_1950-2016_rad5cv25_D2_v2_SON",hilo="low",month1=9,month2=11,cv=0.25,dur=2)

#plot_counts(1950,2016,dir="/short/eg3/asp561/cts.dir/gcyc_out/NCEP1/proj100_highs_rad10cv0.075",fout="NCEP1_UM_anticyclonetrend_1950-2016_rad10cv075_slp1020_D2_v2_SON",hilo="high",month1=9,month2=11,slp=1020,dur=2)


#plot_counts(1950,2016,dir="/short/eg3/asp561/cts.dir/gcyc_out/NCEP1/proj100_rad5cv0.15",fout="NCEP1_UM_cyclonetrend_1950-2016_rad5cv25_D2_v2_DJF",hilo="low",month1=12,month2=2,cv=0.25,dur=2)

#plot_counts(1950,2016,dir="/short/eg3/asp561/cts.dir/gcyc_out/NCEP1/proj100_highs_rad10cv0.075",fout="NCEP1_UM_anticyclonetrend_1950-2016_rad10cv075_slp1020_D2_v2_DJF",hilo="high",month1=12,month2=2,slp=1020,dur=2)


seas=c("","_MAM","_JJA","_SON","_DJF")
m1=c(1,3,6,9,12)
m2=c(12,5,8,11,2)

for(s in 1:5)
{
#plot_counts(1950,2014,dir="/short/eg3/asp561/cts.dir/gcyc_out/20CR/EnsMean/proj100_highs_rad10cv0.075",fout=paste("20CR_UM_anticyclonetrend_1950-2014_rad10cv075_v2",seas[s],sep=""),hilo="high",month1=m1[s],month2=m2[s])

#plot_counts(1950,2014,dir="/short/eg3/asp561/cts.dir/gcyc_out/20CR/EnsMean/proj100_highs_rad10cv0.075",fout=paste("20CR_UM_anticyclonetrend_1950-2014_rad10cv075_slp1020_D2_v2",seas[s],sep=""),hilo="high",month1=m1[s],month2=m2[s],slp=1020,dur=2)

plot_counts(1980,2014,dir="/short/eg3/asp561/cts.dir/gcyc_out/20CR/EnsMean/proj100_highs_rad10cv0.075",fout=paste("20CR_UM_anticyclonetrend_1980-2014_rad10cv075_slp1020_D2_1910_v2",seas[s],sep=""),hilo="high",month1=m1[s],month2=m2[s],slp=1020,dur=2)

plot_counts(1980,2016,dir="/short/eg3/asp561/cts.dir/gcyc_out/NCEP1/proj100_highs_rad10cv0.075",fout=paste("NCEP1_anticyclonetrend_1980-2016_rad10cv075_slp1020_D2_1910_v2",seas[s],sep=""),hilo="high",month1=m1[s],month2=m2[s],slp=1020,dur=2)

plot_counts(1980,2016,dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_highs_rad10cv0.075",fout=paste("ERAI_UM_anticyclonetrend_1980-2016_rad10cv075_slp1020_D2_1910_v2",seas[s],sep=""),hilo="high",month1=m1[s],month2=m2[s],slp=1020,dur=2)
}

