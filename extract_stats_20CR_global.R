years=seq(1851,2014)
library(ncdf4)
library(sp)

stypes1=c("Events","Event.Length","Event.Move","Event.MSLP","Event.CV","Fixes","Fix.MSLP","Fix.CV","Fix.Radius","Fix.Depth","Events.1000","Events.D8")

stypes=c(stypes1,paste0("NH.",stypes1),paste0("SH.",stypes1))

stats=array(NaN,c(length(years),12,length(stypes),59))
dimnames(stats)[[1]]=years
dimnames(stats)[[2]]=1:12
dimnames(stats)[[3]]=stypes
dimnames(stats)[[4]]<-names<-c(1:56,"EnsMean","NCEP1","ERAI")

dir="/short/eg3/asp561/cts.dir/gcyc_out/"
type="proj100_lows_rad5cv0.15/"

outfile="UM_20CR_globalcyclonestats_proj100_rad5cv0.15_500km.RData"

for(y in 1:length(years))
{
  for(n in 1:56)
  {
    print(paste(years[y],n))
    fname=paste(dir,"/20CR/",type,"/tracks_",years[y],"_",n,".dat",sep="")
    read.table(fname, sep="",skip=1)->fixes
    colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Depth","Radius","Up","Vp")
    fixes=data.frame(fixes)
    fixes$Year=floor(fixes$Date/10000)
    fixes$CV=abs(fixes$CV)
    fixes$Depth=abs(fixes$Depth)
    if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
    fixes$Date=(fixes$Date%%10000) + years[y]*10000
    fixes$Month=floor(fixes$Date/100)%%100

    fixes$Move<-NaN
    I=which(fixes$Fix>1)
    if(I[1]==1) I=I[-1]
    for(i in 1:length(I)) fixes$Move[I[i]]=spDistsN1(cbind(fixes$Lon[I[i]],fixes$Lat[I[i]]),cbind(fixes$Lon[I[i]-1],fixes$Lat[I[i]-1]),longlat=T)
    
    x<-rle(fixes$ID)
    events<-data.frame(ID=x$values,Length=x$lengths,Date1=rep(0,length(x$values)),Move=rep(0,length(x$values)))
    for(i in 1:length(events$ID)) 
    {
    events$Date1[i]=min(fixes$Date[fixes$ID==events$ID[i]])
    I=which(fixes$ID==events[i,1])
    events$Move[i]=spDistsN1(cbind(fixes$Lon[min(I)],fixes$Lat[min(I)]),
                           cbind(fixes$Lon[max(I)],fixes$Lat[max(I)]),longlat=T)

    }
    events$Month=floor(events$Date1/100)%%100
    events$MSLP=aggregate(fixes$MSLP,by=list(fixes$ID),FUN=max)[,2]
    events$CV=aggregate(fixes$CV,by=list(fixes$ID),FUN=max)[,2]   
    events$Lat=aggregate(fixes$Lat,by=list(fixes$ID),FUN=mean)[,2]
 
    events=events[events$Move>=500,]
    include<-match(fixes[,1],events[,1])
    J<-which(is.na(include)==0)
    fixes=fixes[J,]
    
    for(m in 1:12)
     for(tt in 0:2)
    {
      I=switch(as.character(tt),
         "0"=which(events$Month==m),
         "1"=which(events$Month==m & events$Lat>0),
         "2"=which(events$Month==m & events$Lat<0))
      lntt=length(stypes1)*tt

      stats[y,m,1+lntt,n]=length(I)
      if(length(I)>0)
      {
        stats[y,m,2+lntt,n]=mean(events$Length[I])
        stats[y,m,3+lntt,n]=mean(events$Move[I])
        stats[y,m,4+lntt,n]=mean(events$MSLP[I])
        stats[y,m,5+lntt,n]=mean(events$CV[I])

        J=switch(as.character(tt),
           "0"=which(fixes$Month==m),
           "1"=which(fixes$Month==m & fixes$Lat>0),
           "2"=which(fixes$Month==m & fixes$Lat<0))

        stats[y,m,6+lntt,n]=length(J)
        stats[y,m,7+lntt,n]=mean(fixes$MSLP[J])
        stats[y,m,8+lntt,n]=mean(fixes$CV[J])
        stats[y,m,9+lntt,n]=mean(fixes$Radius[J])
        stats[y,m,10+lntt,n]=mean(fixes$Depth[J])

        stats[y,m,11+lntt,n]=length(which(events$MSLP[I]<=1000))
        stats[y,m,12+lntt,n]=length(which(events$Length[I]>=8))

      }
      
    }
  }

  stats[y,,,57]=apply(stats[y,,,1:56],c(1,2),mean)

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
  
  fixes$Month=floor(fixes$Date/100)%%100
  fixes$Move<-NaN
  I=which(fixes$Fix>1)
  if(I[1]==1) I=I[-1]
  for(i in 1:length(I)) fixes$Move[I[i]]=spDistsN1(cbind(fixes$Lon[I[i]],fixes$Lat[I[i]]),cbind(fixes$Lon[I[i]-1],fixes$Lat[I[i]-1]),longlat=T)
  
    x<-rle(fixes$ID)
    events<-data.frame(ID=x$values,Length=x$lengths,Date1=rep(0,length(x$values)),Move=rep(0,length(x$values)))
    for(i in 1:length(events$ID))
    {
    events$Date1[i]=min(fixes$Date[fixes$ID==events$ID[i]])
    I=which(fixes$ID==events[i,1])
    events$Move[i]=spDistsN1(cbind(fixes$Lon[min(I)],fixes$Lat[min(I)]),
                           cbind(fixes$Lon[max(I)],fixes$Lat[max(I)]),longlat=T)

    }
    events$Month=floor(events$Date1/100)%%100
    events$MSLP=aggregate(fixes$MSLP,by=list(fixes$ID),FUN=max)[,2]
    events$CV=aggregate(fixes$CV,by=list(fixes$ID),FUN=max)[,2]
    events$Lat=aggregate(fixes$Lat,by=list(fixes$ID),FUN=mean)[,2]

    events=events[events$Move>=500,]
    include<-match(fixes[,1],events[,1])
    J<-which(is.na(include)==0)
    fixes=fixes[J,]
  
  for(m in 1:12)
     for(tt in 0:2)
    {
      I=switch(as.character(tt),
         "0"=which(events$Month==m),
         "1"=which(events$Month==m & events$Lat>0),
         "2"=which(events$Month==m & events$Lat<0))
      lntt=length(stypes1)*tt

      stats[y,m,1+lntt,n]=length(I)
      if(length(I)>0)
      {
        stats[y,m,2+lntt,n]=mean(events$Length[I])
        stats[y,m,3+lntt,n]=mean(events$Move[I])
        stats[y,m,4+lntt,n]=mean(events$MSLP[I])
        stats[y,m,5+lntt,n]=mean(events$CV[I])

        J=switch(as.character(tt),
           "0"=which(fixes$Month==m),
           "1"=which(fixes$Month==m & fixes$Lat>0),
           "2"=which(fixes$Month==m & fixes$Lat<0))

        stats[y,m,6+lntt,n]=length(J)
        stats[y,m,7+lntt,n]=mean(fixes$MSLP[J])
        stats[y,m,8+lntt,n]=mean(fixes$CV[J])
        stats[y,m,9+lntt,n]=mean(fixes$Radius[J])
        stats[y,m,10+lntt,n]=mean(fixes$Depth[J])

        stats[y,m,11+lntt,n]=length(which(events$MSLP[I]<=1000))
        stats[y,m,12+lntt,n]=length(which(events$Length[I]>=8))

      }

    }
    
  }
}

##In case netcdf doesn't work
save(stats,file=paste(dir,outfile,sep=""))

