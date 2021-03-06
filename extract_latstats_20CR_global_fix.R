years=seq(1851,2014)
library(ncdf4)
library(sp)

stypes=c("Fixes","Fix.MSLP","Fix.CV","Fix.Radius","Fix.Depth","Fixes.1000")
lats=seq(-90,89,1)

stats=array(NaN,c(length(lats),length(years),12,length(stypes),59))
dimnames(stats)[[1]]=seq(-89.5,89.5,1)
dimnames(stats)[[2]]=years
dimnames(stats)[[3]]=1:12
dimnames(stats)[[4]]=stypes
dimnames(stats)[[5]]<-names<-c(1:56,"EnsMean","NCEP1","ERAI")

dir="/short/eg3/asp561/cts.dir/gcyc_out/"
type="proj100_lows_rad5cv0.15/"

infile="UM_20CR_globalcyclones_statsbylat_proj100_rad5cv0.15_500km_BAD.RData"
outfile="UM_20CR_globalcyclones_statsbylat_proj100_rad5cv0.15_500km.RData"
load(paste0(dir,infile))
dimnames(stats)[[4]]=stypes

for(y in 1:length(years))
{

for(n in 58:59)
{
  fname=paste(dir,"/",names[n],"/",type,"/tracks_",years[y],".dat",sep="")
  if(!file.exists(fname)) next
  read.table(fname, sep="",skip=1)->fixes
  colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Depth","Radius","Up","Vp")
  fixes=data.frame(fixes)
  fixes$Year=floor(fixes$Date/10000)
  fixes$CV=abs(fixes$CV)
  fixes$Depth=abs(fixes$Depth)
  if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
  fixes$Date=(fixes$Date%%10000) + years[y]*10000
  fixes$Lat2=floor(fixes$Lat)
  
  fixes$Month=floor(fixes$Date/100)%%100
  fixes$Move<-NaN
  I=which(fixes$Fix>1)
  if(I[1]==1) I=I[-1]
  for(i in 1:length(I)) fixes$Move[I[i]]=spDistsN1(cbind(fixes$Lon[I[i]],fixes$Lat[I[i]]),cbind(fixes$Lon[I[i]-1],fixes$Lat[I[i]-1]),longlat=T)

  fixes$Closed=0
  fixes$Closed[fixes$Open==0 | fixes$Open==10]=1

  x<-rle(fixes$ID)
    events<-data.frame(ID=x$values,Length=x$lengths,Date1=rep(0,length(x$values)),Move=rep(0,length(x$values)),Closed=rep(0,length(x$values)))
    for(i in 1:length(events$ID))
    {
    events$Date1[i]=min(fixes$Date[fixes$ID==events$ID[i] & fixes$Location==1])
    I=which(fixes$ID==events[i,1])
    events$Move[i]=spDistsN1(cbind(fixes$Lon[min(I)],fixes$Lat[min(I)]),
                           cbind(fixes$Lon[max(I)],fixes$Lat[max(I)]),longlat=T)

    events$Closed[i]=max(fixes$Closed[I])
    events$Length[i]=sum(fixes$Closed[I]) ## Only count closed fixes in length
    }

    events=events[events$Move>=500 & events$Closed==1,]
    include<-match(fixes[,1],events[,1])
    J<-which(is.na(include)==0 & fixes$Closed==1)
    fixes=fixes[J,]

    for(m in 1:12)
     for(ll in 1:length(lats))
    {
      J=which(fixes$Month==m & fixes$Lat2==lats[ll])
      stats[ll,y,m,1,n]=length(J)
      if(length(I)>0)
      {
        stats[ll,y,m,2,n]=mean(fixes$MSLP[J])
        stats[ll,y,m,3,n]=mean(fixes$CV[J])
        stats[ll,y,m,4,n]=mean(fixes$Radius[J])
        stats[ll,y,m,5,n]=mean(fixes$Depth[J])
        stats[ll,y,m,6,n]=length(which(fixes$MSLP[J]<=1000))
      }

    }  
    
  }
}

##In case netcdf doesn't work
save(stats,file=paste(dir,outfile,sep=""))

