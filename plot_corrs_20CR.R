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

plot_corrs_many<-function(year1,year2,dir="gcyc_out",name="proj100_highs_rad10cv0.075",cv=NA,dur=NA,slp=NA,month1=1,month2=12,fout=NA,hilo="low",members=56,
                          inames="SOI",ifiles=NaN)
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

### Step one: get the index file

lat=seq(-85,85,10)
lon=seq(5,355,10)  ### Can always combine into bigger cells later
if(month2>=month1) systems<-array(0,c(length(lon),length(lat),length(years),members+1,4)) else systems<-array(0,c(length(lon),length(lat),length(years)-1,members+1,4))
dimnames(systems)[[5]]=c("Freq","MSLP","CV","LonMove")

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
  systems[,,y,em,1]=table(factor(fixes$Lon2[I],levels=seq(0,350,10)),factor(fixes$Lat2[I],levels=seq(-90,80,10)))
  tmp2=aggregate(fixes[I,8:9],by=list(fixes$Lon2[I],fixes$Lat2[I]),FUN=mean)
  jlat=match(tmp2[,2],seq(-90,80,10))
  ilon=match(tmp2[,1],seq(0,350,10))
  for(nn in which(!is.na(jlat))) systems[ilon[nn],jlat[nn],y,em,2:3]=as.numeric(tmp2[nn,3:4])

#  I=which(fixes$Month>=month1 & fixes$Month<=month2 & !is.na(fixes$LonMove))
#  tmp2=aggregate(fixes$LonMove[I],by=list(fixes$Lon2[I],fixes$Lat2[I]),FUN=mean)
#  jlat=match(tmp2[,2],seq(-90,80,10))
#  ilon=match(tmp2[,1],seq(0,350,10))
#  for(nn in which(!is.na(jlat))) systems[ilon[nn],jlat[nn],y,em,4]=as.numeric(tmp2[nn,3])


 } else {
  if(y==1)
  {
   I=which(fixes$Month>=month1)
   store=fixes[I,]
  } else {
   I=which(fixes$Month<=month2)
   fixes2=rbind(store,fixes[I,])

   systems[,,y-1,em,1]=table(factor(fixes2$Lon2,levels=seq(0,350,10)),factor(fixes2$Lat2,levels=seq(-90,80,10)))
   tmp2=aggregate(fixes2[,8:9],by=list(fixes2$Lon2,fixes2$Lat2),FUN=mean)
   jlat=match(tmp2[,2],seq(-90,80,10))
   ilon=match(tmp2[,1],seq(0,350,10))
   for(nn in which(!is.na(jlat))) systems[ilon[nn],jlat[nn],y-1,em,2:3]=as.numeric(tmp2[nn,3:4])

#  I=which(!is.na(fixes2$LonMove))
#  tmp2=aggregate(fixes2$LonMove[I],by=list(fixes2$Lon2[I],fixes2$Lat2[I]),FUN=mean)
#  jlat=match(tmp2[,2],seq(-90,80,10))
#  ilon=match(tmp2[,1],seq(0,350,10))
#  for(nn in which(!is.na(jlat))) systems[ilon[nn],jlat[nn],y-1,em,4]=as.numeric(tmp2[nn,3])

   I=which(fixes$Month>=month1)
   store=fixes[I,]
  }
}
} # End year loop

systems[,,,members+1,]=apply(systems[,,,1:members,],c(1,2,3,5),mean)

print("Calculating means")

if(month2<month1) years=seq(year1,year2-1,1)

### Correlations
meanfreq=apply(systems[,,,members+1,1],c(1,2),mean)

for(ii in 1:length(inames))
{
iname=inames[ii]

if(!is.na(ifiles)) ifile=ifiles[ii] else {
ifile=switch(iname,
                              SOI = "/short/eg3/asp561/Timeseries/SOI_18762016.csv",
                              SAM = "/short/eg3/asp561/Timeseries/SAM_NOAA_19792016.csv",
                              STRI = "/short/eg3/asp561/Timeseries/STRI_18902016.csv",
                              STRP = "/short/eg3/asp561/Timeseries/STRP_18902016.csv",
                              Hadley.SH = "/short/eg3/asp561/Timeseries/HadleyCell.SH.csv",
                              Hadley.AsiaPac = "/short/eg3/asp561/Timeseries/HadleyCell.AsiaPac.csv",
                              DMI = "/short/eg3/asp561/Timeseries/NOAA.DMI.monthly.csv",
                              AOI = "/short/eg3/asp561/Timeseries/NOAA.AOI.csv",
                              NAO = "/short/eg3/asp561/Timeseries/NOAA.NAO.csv",
                              PNA = "/short/eg3/asp561/Timeseries/NOAA.PNA.csv",
                              IPO = "/short/eg3/asp561/Timeseries/Henley.IPO.ERSST4.csv",
                              IPOf = "/short/eg3/asp561/Timeseries/Henley.IPO.ERSST4.filt.csv",
                              NaN)
if(is.na(ifile)) stop("No file provided for index name") 
}

index=read.csv(ifile)
if(min(index$Year)>year1 | max(index$Year)<year2) stop("Index has shorter length of record than selected years")

if(month2>=month1) index2=apply(index[which(index$Year>=year1 & index$Year<=year2),(month1:month2)+1],1,mean) else {
  tmp=cbind(index[which(index$Year>=year1 & index$Year<year2),(month1+1):13],index[which(index$Year>year1 & index$Year<=year2),2:(month2+1)])
  index2=apply(tmp,1,mean)
}
cyccorr<-array(NaN,c(length(lon),length(lat),members+1,3,2))

for(i in 1:length(lon))
  for(j in 1:length(lat))
   for(em in 1:(members+1))
    for(k in 1:3)
     if(meanfreq[i,j]>=5)
     {
      a=cor.test(systems[i,j,,em,k],index2,use="pairwise.complete.obs")
      cyccorr[i,j,em,k,1]=a$estimate
      cyccorr[i,j,em,k,2]=a$p.value      
     }

### Plot

type=c("freq","MSLP","CV","LonMove")
breaks=seq(-1,1,0.1)
cols=col_anom(length(breaks)-1)

tmp=abs(cyccorr[,,members+1,1,1])
sigthresh=min(tmp[cyccorr[,,members+1,1,2]<0.05],na.rm=T)

for(k in 1:3)
{
print("Plotting")

if(is.na(fout)) fname=paste(hilo,"s_",iname,"corr_all_",type[k],"_",year1,"_",year2,"_global.pdf",sep="") else 
  fname=paste(fout,"_",iname,"corr_",type[k],"_global.pdf",sep="")

pdf(file=fname,width=14,height=10)
layout(cbind(c(1,3,5),c(2,4,5)),height=c(1,1,0.15))
par(mar=c(2,2,4,1))

image(lon,lat,cyccorr[,,members+1,k,1],breaks=breaks,col=cols,xlab="",ylab="",
          main=paste("Correlation of Ensemble Mean"))
map('world2',add=T)
sigmask=which(cyccorr[,,members+1,k,2]<0.05 & meanfreq>=5,arr.ind=T)
points(lon[sigmask[,1]],lat[sigmask[,2]],col="black",pch=19,cex=0.6)

tmp=apply(cyccorr[,,1:members,k,1],c(1,2),mean)
image(lon,lat,tmp,breaks=breaks,col=cols,xlab="",ylab="",
          main=paste("Average correlation"))
map('world2',add=T)
sigmask=which(apply(cyccorr[,,1:members,k,2]<0.05,c(1,2),sum)>=(0.95*members) & meanfreq>=5,arr.ind=T)
points(lon[sigmask[,1]],lat[sigmask[,2]],col="black",pch=19,cex=0.6)

tmp=apply(cyccorr[,,1:members,k,1],c(1,2),max)
image(lon,lat,tmp,breaks=breaks,col=cols,xlab="",ylab="",
          main=paste("Most positive correlation"))
map('world2',add=T)
sigmask=which(abs(tmp)>sigthresh & meanfreq>=5,arr.ind=T)
points(lon[sigmask[,1]],lat[sigmask[,2]],col="black",pch=19,cex=0.6)

tmp=apply(cyccorr[,,1:members,k,1],c(1,2),min)
image(lon,lat,tmp,breaks=breaks,col=cols,xlab="",ylab="",
          main=paste("Most negative correlation"))
map('world2',add=T)
sigmask=which(abs(tmp)>sigthresh & meanfreq>=5,arr.ind=T)
points(lon[sigmask[,1]],lat[sigmask[,2]],col="black",pch=19,cex=0.6)

ColorBar(breaks,cols,subsampleg=1,vert=F)
dev.off()

}
}
}
seas=c("","_MAM","_JJA","_SON","_DJF","_NDJFMA","_MJJASO")
m1=c(1,3,6,9,12,11,5)
m2=c(12,5,8,11,2,4,10)

for(s in c(1,6,7))
    {

plot_corrs_many(1950,2010,dir="/short/eg3/asp561/cts.dir/gcyc_out/",name="proj100_lows_rad5cv0.15",fout=paste("UM_cyclonecorr_20CRens_1950-2010_rad5cv15_D2",seas[s],sep=""),hilo="low",month1=m1[s],month2=m2[s],dur=2,inames=c("IPO","IPOf","SOI"))

plot_corrs_many(1950,2010,dir="/short/eg3/asp561/cts.dir/gcyc_out/",name="proj100_highs_rad10cv0.075",fout=paste("UM_anticyclonecorr_20CRens_1950-2010_rad10cv075_D2",seas[s],sep=""),hilo="high",month1=m1[s],month2=m2[s],dur=2,inames=c("IPO","IPOf","SOI"))

plot_corrs_many(1910,2010,dir="/short/eg3/asp561/cts.dir/gcyc_out/",name="proj100_lows_rad5cv0.15",fout=paste("UM_cyclonecorr_20CRens_1910-2010_rad5cv15_D2",seas[s],sep=""),hilo="low",month1=m1[s],month2=m2[s],dur=2,inames=c("IPO","IPOf","SOI"))

plot_corrs_many(1910,2010,dir="/short/eg3/asp561/cts.dir/gcyc_out/",name="proj100_highs_rad10cv0.075",fout=paste("UM_anticyclonecorr_20CRens_1910-2010_rad10cv075_D2",seas[s],sep=""),hilo="high",month1=m1[s],month2=m2[s],dur=2,inames=c("IPO","IPOf","SOI"))
    }
