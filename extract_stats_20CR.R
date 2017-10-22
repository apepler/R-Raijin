years=seq(1851,2014)
library(ncdf4)
library(sp)
stypes=c("Aust.Events.D2","Aust.Length","Aust.Hours","Aust.Hours.1020","Aust.Lat","Aust.MSLP","Aust.CV",
         "Aust.Depth","Aust.Radius","Aust.Up","Aust.Vp","Aust.LonMove","Aust.Move",
         "SAust.Events.D2","SAust.Length","SAust.Hours","SAust.Hours.1020","SAust.Lat","SAust.MSLP","SAust.CV","SAust.LonMove","Aust.Events.D8","SAust.Events.D8")

stats=array(NaN,c(length(years),12,length(stypes),57))
dimnames(stats)[[1]]=years
dimnames(stats)[[2]]=1:12
dimnames(stats)[[3]]=stypes
dimnames(stats)[[4]]=c(1:56,"EnsMean")

dir="/short/eg3/asp561/cts.dir/gcyc_out/20CR/"
type="proj100_highs_rad10cv0.075/"

outfile="UM_20CR_anticyclonestats_proj100_rad10cv0.075_500km_long.RData"

for(y in 1:length(years))
{
  for(n in 1:56)
  {
    print(paste(years[y],n))
    fname=paste(dir,type,"/tracks_",years[y],"_",n,".dat",sep="")
    read.table(fname, sep="",skip=1)->fixes
    colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Depth","Radius","Up","Vp")
    fixes=data.frame(fixes)
    fixes$Year=floor(fixes$Date/10000)
    fixes$CV=abs(fixes$CV)
    fixes$Depth=abs(fixes$Depth)
    if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
    fixes$Date=(fixes$Date%%10000) + years[y]*10000

    fixes$Month=floor(fixes$Date/100)%%100
    fixes$Location2<-fixes$Location<-0
    I<-which(fixes$Lon>=110 & fixes$Lon<=155 & fixes$Lat<(-10) & fixes$Lat>=-45)
    fixes$Location[I]<-1
    I=which(fixes$Lat>=-40 & fixes$Lat<=-30 & fixes$Lon>=110 & fixes$Lon<=155)
    fixes$Location2[I]=1

    fixes$Move<-fixes$LonMove<-NaN
    I=which(fixes$Fix!=1)
    if(I[1]==1) I=I[-1]
    fixes$LonMove[I]=fixes$Lon[I]-fixes$Lon[I-1]
    for(i in 1:length(I)) fixes$Move[I[i]]=spDistsN1(cbind(fixes$Lon[I[i]],fixes$Lat[I[i]]),cbind(fixes$Lon[I[i]-1],fixes$Lat[I[i]-1]),longlat=T)
    
    x<-rle(fixes$ID)
    events<-data.frame(ID=x$values,Length=x$lengths,Location=rep(0,length(x$values)),Location2=rep(0,length(x$values)),Date1=rep(0,length(x$values)),Move=rep(0,length(x$values)))
    for(i in 1:length(events$ID)) 
    {
    events$Location[i]=sum(fixes$Location[fixes$ID==events$ID[i]])
    events$Location2[i]=sum(fixes$Location2[fixes$ID==events$ID[i]])
    events$Date1[i]=min(fixes$Date[fixes$ID==events$ID[i] & fixes$Location==1])
    I=which(fixes$ID==events[i,1])
    events$Move[i]=spDistsN1(cbind(fixes$Lon[min(I)],fixes$Lat[min(I)]),
                           cbind(fixes$Lon[max(I)],fixes$Lat[max(I)]),longlat=T)

    }
    events$Month=floor(events$Date1/100)%%100
    events$MSLP=aggregate(fixes$MSLP,by=list(fixes$ID),FUN=max)
    
    events=events[events$Move>=500 & events$Location>=1,]
    include<-match(fixes[,1],events[,1])
    J<-which(is.na(include)==0)
    fixes=fixes[J,]
    
    for(m in 1:12)
    {
      I=which(events$Month==m)
      
      stats[y,m,1,n]=length(I)
      if(length(I)>0)
      {
        stats[y,m,2,n]=mean(events$Location[I])
        J=which(fixes$ID%in%events$ID[I] & fixes$Location==1)
        stats[y,m,3,n]=length(J)
        J2=which(fixes$ID%in%events$ID[I] & fixes$Location==1 & fixes$MSLP>=1020)
        stats[y,m,4,n]=length(J2)
        stats[y,m,5,n]=mean(fixes$Lat[J])
        stats[y,m,6,n]=mean(fixes$MSLP[J])
        stats[y,m,7,n]=mean(fixes$CV[J])
        stats[y,m,8,n]=mean(fixes$Depth[J])
        stats[y,m,9,n]=mean(fixes$Radius[J])
        stats[y,m,10,n]=mean(fixes$Up[J])
        stats[y,m,11,n]=mean(fixes$Vp[J])
        stats[y,m,12,n]=mean(fixes$LonMove[J],na.rm=T)
        stats[y,m,13,n]=mean(fixes$Move[J],na.rm=T)
        
        I2=which(events$Location2[I]>0)
        stats[y,m,14,n]=length(I2)
        if(length(I)>0)
        {
          stats[y,m,15,n]=mean(events$Location2[I[I2]])
          J=which(fixes$ID%in%events$ID[I[I2]] & fixes$Location2==1)
          stats[y,m,16,n]=length(J)
          J2=which(fixes$ID%in%events$ID[I[I2]] & fixes$Location2==1 & fixes$MSLP>=1020)
          stats[y,m,17,n]=length(J2)
          
          stats[y,m,18,n]=mean(fixes$Lat[J])
          stats[y,m,19,n]=mean(fixes$MSLP[J])
          stats[y,m,20,n]=mean(fixes$CV[J])
          stats[y,m,21,n]=mean(fixes$LonMove[J],na.rm=T)
        }

       stats[y,m,22,n]=length(which(events$Length[I]>=8))
       stats[y,m,23,n]=length(which(events$Length[I]>=8 & events$Location2[I]==1))
      }
      
    }
  }
  n=57
  fname=paste(dir,"EnsMean/",type,"/tracks_",years[y],".dat",sep="")
  read.table(fname, sep="",skip=1)->fixes
  colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Depth","Radius","Up","Vp")
  fixes=data.frame(fixes)
  fixes$Year=floor(fixes$Date/10000)
  fixes$CV=abs(fixes$CV)
  fixes$Depth=abs(fixes$Depth)
  if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
  fixes$Date=(fixes$Date%%10000) + years[y]*10000
  
  fixes$Month=floor(fixes$Date/100)%%100
  fixes$Location2<-fixes$Location<-0
  I<-which(fixes$Lon>=110 & fixes$Lon<=155 & fixes$Lat<(-10) & fixes$Lat>=-45)
  fixes$Location[I]<-1
  I=which(fixes$Lat>=-40 & fixes$Lat<=-30 & fixes$Lon>=110 & fixes$Lon<=155)
  fixes$Location2[I]=1

    fixes$Move<-fixes$LonMove<-NaN
    I=which(fixes$Fix!=1)
    if(I[1]==1) I=I[-1]
    fixes$LonMove[I]=fixes$Lon[I]-fixes$Lon[I-1]
    for(i in 1:length(I)) fixes$Move[I[i]]=spDistsN1(cbind(fixes$Lon[I[i]],fixes$Lat[I[i]]),cbind(fixes$Lon[I[i]-1],fixes$Lat[I[i]-1]),longlat=T)
  
  x<-rle(fixes$ID)
  events<-data.frame(ID=x$values,Length=x$lengths,Location=rep(0,length(x$values),Location2=rep(0,length(x$values))))
  for(i in 1:length(events$ID)) events$Location[i]=sum(fixes$Location[fixes$ID==events$ID[i]])
  for(i in 1:length(events$ID)) events$Location2[i]=sum(fixes$Location2[fixes$ID==events$ID[i]])
  events$Date1<-0
  for(i in 1:length(events$ID)) events$Date1[i]=min(fixes$Date[fixes$ID==events$ID[i] & fixes$Location==1])
  events$Month=floor(events$Date1/100)%%100
  
  events=events[events$Length>=2 & events$Location>=1,]
  include<-match(fixes[,1],events[,1])
  J<-which(is.na(include)==0)
  fixes=fixes[J,]
  
  for(m in 1:12)
  {
    I=which(events$Month==m)
    
    stats[y,m,1,n]=length(I)
    if(length(I)>0)
    {
      stats[y,m,2,n]=mean(events$Location[I])
      J=which(fixes$ID%in%events$ID[I] & fixes$Location==1)
      stats[y,m,3,n]=length(J)
      J2=which(fixes$ID%in%events$ID[I] & fixes$Location==1 & fixes$MSLP>=1020)
      stats[y,m,4,n]=length(J2)
      stats[y,m,5,n]=mean(fixes$Lat[J])
      stats[y,m,6,n]=mean(fixes$MSLP[J])
      stats[y,m,7,n]=mean(fixes$CV[J])
      stats[y,m,8,n]=mean(fixes$Depth[J])
      stats[y,m,9,n]=mean(fixes$Radius[J])
      stats[y,m,10,n]=mean(fixes$Up[J])
      stats[y,m,11,n]=mean(fixes$Vp[J])
      stats[y,m,12,n]=mean(fixes$LonMove[J],na.rm=T)
      stats[y,m,13,n]=mean(fixes$Move[J],na.rm=T)
      
      I2=which(events$Location2[I]>0)
      stats[y,m,14,n]=length(I2)
      if(length(I)>0)
      {
        stats[y,m,15,n]=mean(events$Location2[I[I2]])
        J=which(fixes$ID%in%events$ID[I[I2]] & fixes$Location2==1)
        stats[y,m,16,n]=length(J)
        J2=which(fixes$ID%in%events$ID[I[I2]] & fixes$Location2==1 & fixes$MSLP>=1020)
        stats[y,m,17,n]=length(J2)
        
        stats[y,m,18,n]=mean(fixes$Lat[J])
        stats[y,m,19,n]=mean(fixes$MSLP[J])
        stats[y,m,20,n]=mean(fixes$CV[J])
        stats[y,m,21,n]=mean(fixes$LonMove[J],na.rm=T)
      }
       stats[y,m,22,n]=length(which(events$Length[I]>=8))
       stats[y,m,23,n]=length(which(events$Length[I]>=8 & events$Location2[I]==1))

    }
    
  }
}

##In case netcdf doesn't work
save(stats,file=paste(dir,outfile,sep=""))

