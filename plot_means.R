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

plot_counts<-function(year1,year2,dir="gcyc_out",anti="proj100_highs_rad10cv0.075",cyc="proj100_lows_rad5cv0.15",hcv=NA,lcv=NA,dur=NA,hslp=NA,lslp=NA,month1=1,month2=12,fout=NA)
{
years=seq(year1,year2,1)

lat=seq(-89.5,89.5)
lon=seq(0,359.5)  ### Can always combine into bigger cells later
if(month2>=month1) systems<-array(0,c(length(lon),length(lat),length(years),2)) else systems<-array(0,c(length(lon),length(lat),length(years)-1,2))
dimnames(systems)[[4]]=c("Cyclones","Anticyclones")
name=c(cyc,anti)
cv=c(lcv,hcv)

for(n in 1:2)
for(y in 1:length(years))
{
print(years[y])
fname=paste(dir,name[n],"/tracks_",years[y],".dat",sep="")
read.table(fname, sep="",skip=1)->fixes
colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Depth","Radius","Up","Vp")
fixes$Year=floor(fixes$Date/10000)
fixes$CV=abs(fixes$CV)
fixes$Depth=abs(fixes$Depth)
if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
fixes$Month=floor(fixes$Date/100)%%100
fixes$Lat2=floor(fixes$Lat) ## Let's skip the first step, and just make it the 10 degrees
fixes$Lon2=floor((fixes$Lon%%360))
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
if(!is.na(cv[n])) fixes=fixes[fixes$CV>=cv[n],]

if(n==1 & !is.na(lslp)) fixes=fixes[fixes$MSLP<=lslp,]
if(n==2 & !is.na(hslp)) fixes=fixes[fixes$MSLP>=hslp,]

if(month2>=month1)
 {
  I=which(fixes$Month>=month1 & fixes$Month<=month2)
  systems[,,y,n]=table(factor(fixes$Lon2[I],levels=seq(0,359)),factor(fixes$Lat2[I],levels=seq(-90,89)))
 } else {
  if(y==1)
  {
   I=which(fixes$Month>=month1)
   store=fixes[I,]
  } else {
   I=which(fixes$Month<=month2)
   fixes2=rbind(store,fixes[I,])
   systems[,,y-1,n]=table(factor(fixes2$Lon2,levels=seq(0,359)),factor(fixes2$Lat2,levels=seq(-90,89)))

   I=which(fixes$Month>=month1)
   store=fixes[I,]
  }
}
} # End year loop

print("Calculating means")
if(month2<month1) years=seq(year1,year2-1,1)

### Plot the average & the linear trend

### Mean frequency

meanfreq=apply(systems,c(1,2,4),mean,na.rm=T)

### Plot

type=c("Cyclone Frequency (Cyclones/deg^2)","Anticyclone Frequency (Anticyclones/deg^2)")

if(month1==1 & month2==12) breaks=c(0,0.05,seq(0.5,2,0.5),3:5,1000) else breaks=c(0,0.05,seq(0.25,1,0.25),seq(1.5,3,0.5),1000)
print("Plotting")
col1=col_val(length(breaks)-1)

if(is.na(fout)) fname=paste("Systems_mean_",year1,"_",year2,".pdf",sep="") else fname=paste(fout,".pdf",sep="")

pdf(file=fname,width=7,height=10)
layout(c(1,2,3),height=c(1,1,0.15))
par(mar=c(2,2,4,1))

for(n in 1:2)
{
image(lon,lat,meanfreq[,,n],breaks=breaks,col=col1,xlab="",ylab="",
          main=type[n])
map('world2',add=T)
#contour(lon,lat,meanfreq[,,n],levels=breaks,add=T,lwd=2,col="black",drawlabels=F)
}

ColorBar(breaks,col1,subsampleg=1,vert=F)
dev.off()
}


seas=c("","_MAM","_JJA","_SON","_DJF","_MJJASO","_NDJFMA")
m1=c(1,3,6,9,12,5,11)
m2=c(12,5,8,11,2,10,4)

for(s in 1:7)
{
#plot_counts(1980,2016,dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/",anti="proj100_highs_rad10cv0.075",cyc="proj100_lows_rad5cv0.15",fout=paste("ERAI_UM_meanfreq_1980-2016_D5",seas[s],sep=""),month1=m1[s],month2=m2[s],dur=5)

#plot_counts(1980,2016,dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/",anti="proj100_highs_rad10cv0.075",cyc="proj100_lows_rad5cv0.15",fout=paste("ERAI_UM_meanfreq_1980-2016",seas[s],sep=""),month1=m1[s],month2=m2[s],dur=2)

#plot_counts(1980,2016,dir="/short/eg3/asp561/cts.dir/gcyc_out/NCEP1/",anti="proj100_highs_rad10cv0.075_v2",cyc="proj100_lows_rad5cv0.15_v2",fout=paste("NCEP1_UM_meanfreq_1980-2016_D5",seas[s],sep=""),month1=m1[s],month2=m2[s],dur=5)

plot_counts(1980,2016,dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/",anti="proj100_highs_rad10cv0.075_notopo",cyc="proj100_lows_rad5cv0.15_notopo",fout=paste("ERAI_UMnotopo_meanfreq_1980-2016a",seas[s],sep=""),month1=m1[s],month2=m2[s],dur=2)
}

