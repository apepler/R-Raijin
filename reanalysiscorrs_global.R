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
sdirs=c("20CR/","NCEP1/","ERAI/")
snames=c("20CRensemble","NCEP1","ERAI")
} else {
sdirs=c("20CR/","NCEP1/")
snames=c("20CRensemble","NCEP1")
}

lat=seq(-85,85,10)
lon=seq(5,355,10)  ### Can always combine into bigger cells later
if(month2>=month1) systems<-array(0,c(length(lon),length(lat),length(years),members)) else systems<-array(0,c(length(lon),length(lat),length(years)-1,members))

systems2=systems[,,,1:length(snames)]

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
fixes$Lat2=10*floor(fixes$Lat/10) ## Let's skip the first step, and just make it the 10 degrees
fixes$Lon2=10*floor((fixes$Lon%%360)/10)

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
  systems[,,y,em]=table(factor(fixes$Lon2[I],levels=seq(0,350,10)),factor(fixes$Lat2[I],levels=seq(-90,80,10)))

 } else {
  if(y==1)
  {
   I=which(fixes$Month>=month1)
   store=fixes[I,]
  } else {
   I=which(fixes$Month<=month2)
   fixes2=rbind(store,fixes[I,])

   systems[,,y-1,em]=table(factor(fixes2$Lon2,levels=seq(0,350,10)),factor(fixes2$Lat2,levels=seq(-90,80,10)))

   I=which(fixes$Month>=month1)
   store=fixes[I,]
  }
}
} # End year loop

print("Calculating means")
#if(month2<month1) years=seq(year1,year2-1,1)

### Plot the average & the linear trend

### Now, add the other ones

systems2[,,,1]=apply(systems[,,,1:members],c(1,2,3),mean)

for(y in 1:length(years))
for(em in 2:length(sdirs))
{
if(em==2) fname=paste(dir,sdirs[em],name,"_v2/tracks_",years[y],".dat",sep="") else fname=paste(dir,sdirs[em],name,"/tracks_",years[y],".dat",sep="")
read.table(fname, sep="",skip=1)->fixes
colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Depth","Radius","Up","Vp")
#fixes=fixes[fixes$Lat<0,]
fixes$Year=floor(fixes$Date/10000)
fixes$CV=abs(fixes$CV)
fixes$Depth=abs(fixes$Depth)
if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
fixes$Month=floor(fixes$Date/100)%%100
fixes$Lat2=10*floor(fixes$Lat/10) ## Let's skip the first step, and just make it the 10 degrees
fixes$Lon2=10*floor((fixes$Lon%%360)/10)

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
  systems2[,,y,em]=table(factor(fixes$Lon2[I],levels=seq(0,350,10)),factor(fixes$Lat2[I],levels=seq(-90,80,10)))

 } else {
  if(y==1)
  {
   I=which(fixes$Month>=month1)
   store=fixes[I,]
  } else {
   I=which(fixes$Month<=month2)
   fixes2=rbind(store,fixes[I,])

   systems2[,,y-1,em]=table(factor(fixes2$Lon2,levels=seq(0,350,10)),factor(fixes2$Lat2,levels=seq(-90,80,10)))

   I=which(fixes$Month>=month1)
   store=fixes[I,]
  }
}
} # End year loop

if(month2<month1) years=seq(year1,year2-1,1)

### Linear trend
meanfreq=apply(systems2,c(1,2,4),mean)

if(length(snames)>2)
{
cyccorr<-array(NaN,c(length(lon),length(lat),3,2))
dimnames(cyccorr)[[3]]<-cornames<-c(paste(snames[1],"vs",snames[2]),paste(snames[1],"vs",snames[3]),paste(snames[2],"vs",snames[3]))
} else cyccorr<-array(NaN,c(length(lon),length(lat),1,2))

for(i in 1:length(lon))
  for(j in 1:length(lat))
     {

     if(meanfreq[i,j,1]>=5 & meanfreq[i,j,2]>=5)
     {
      a=cor.test(systems2[i,j,,1],systems2[i,j,,2])
      cyccorr[i,j,1,1]=a$estimate
      cyccorr[i,j,1,2]=a$p.value      
     }
     if(length(snames)>2)
     {
     if(meanfreq[i,j,1]>=5 & meanfreq[i,j,3]>=5)
     {
      a=cor.test(systems2[i,j,,1],systems2[i,j,,3])
      cyccorr[i,j,2,1]=a$estimate
      cyccorr[i,j,2,2]=a$p.value
     }
     if(meanfreq[i,j,2]>=5 & meanfreq[i,j,3]>=5)
     {
      a=cor.test(systems2[i,j,,2],systems2[i,j,,3])
      cyccorr[i,j,3,1]=a$estimate
      cyccorr[i,j,3,2]=a$p.value
     }
     }
     }


### Plot

breaks=c(-1,0,seq(0.4,0.9,0.1),0.95,1)
cols=col_val(length(breaks)-1)

if(is.na(fout)) fname=paste("Systems_frqtrend_",year1,"_",year2,"_global.pdf",sep="") else fname=paste(fout,"_all_freq_global.pdf",sep="")

if(length(snames)==2)
{
pdf(file=fname,width=8,height=5)
layout(cbind(1,2),width=c(1,0.13))
par(mar=c(3,3,4,1))
image(lon,lat,cyccorr[,,1,1],breaks=breaks,col=cols,xlab="",ylab="",
          main=paste(snames[1],"vs",snames[2]))
map('world2',add=T)
#sigmask=which(cyccorr[,,1,2]<0.05,arr.ind=T)
#points(lon[sigmask[,1]],lat[sigmask[,2]],col="black",pch=19,cex=0.6)
ColorBar(breaks,cols,subsampleg=1,vert=F)
dev.off()

} else {

pdf(file=fname,width=14,height=10)
layout(cbind(c(1,2,5),c(3,4,5)),height=c(1,1,0.15))
par(mar=c(2,2,4,1))

for(s in 1:3)
{
image(lon,lat,cyccorr[,,s,1],breaks=breaks,col=cols,xlab="",ylab="",
          main=cornames[s])
map('world2',add=T)
#sigmask=which(cyccorr[,,s,2]<0.05,arr.ind=T)
#points(lon[sigmask[,1]],lat[sigmask[,2]],col="black",pch=19,cex=0.6)
}
plot.new()
ColorBar(breaks,cols,subsampleg=1,vert=F)
dev.off()
}

}

seas=c("","_MAM","_JJA","_SON","_DJF","_NDJFMA","_MJJASO")
m1=c(1,3,6,9,12,11,5)
m2=c(12,5,8,11,2,4,10)

for(s in 1:7)
{
  
plot_counts_many(1980,2014,dir="/short/eg3/asp561/cts.dir/gcyc_out/",name="proj100_highs_rad10cv0.075",fout=paste("UM_anticyclone_reanalcomp_1980-2014_rad10cv075_D2",seas[s],sep=""),hilo="high",month1=m1[s],month2=m2[s],dur=2)

plot_counts_many(1980,2014,dir="/short/eg3/asp561/cts.dir/gcyc_out/",name="proj100_lows_rad5cv0.15",fout=paste("UM_cyclone_reanalcomp_1980-2014_rad5cv15_D2",seas[s],sep=""),hilo="low",month1=m1[s],month2=m2[s],dur=2)

}
