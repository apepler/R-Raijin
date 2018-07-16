years=seq(1851,2016)
library(ncdf4)
library(sp)

stypes=c("Events","Events.1020","Event.Length","Event.Move","Fixes","Fix.MSLP","Fix.CV","Fix.Radius","Fix.Depth","Fix.Move","Fix.LonMove","Fix.Latitude")
regnames=c("Aust","BigTas","GAB","SEIO")

dir="/short/eg3/asp561/cts.dir/gcyc_out/"
type="proj100_highs_rad10cv0.075/"

infile="UM_20CR_anticyclonestats_austregions3_proj100_rad10cv0.075.RData"
outfile="UM_20CR_anticyclonestats_austregions3_proj100_rad10cv0.075_fix.RData"

load(file=paste(dir,infile,sep=""))
stats[,,,,58:59]=NaN

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
  fixes$Month=floor(fixes$Date/100)%%100
  fixes$Closed=0
  fixes$Closed[fixes$Open==0 | fixes$Open==10]=1
  fixes$Radius[fixes$Radius==0 | fixes$Radius > 100] = NaN
  fixes$Depth[fixes$Depth==0 | fixes$Depth > 100] = NaN

  fixes$Move<-NaN
  I=which(fixes$Fix>1)
  if(I[1]==1) I=I[-1]
  for(i in 1:length(I)) fixes$Move[I[i]]=spDistsN1(cbind(fixes$Lon[I[i]],fixes$Lat[I[i]]),cbind(fixes$Lon[I[i]-1],fixes$Lat[I[i]-1]),longlat=T)
  fixes$LonMove=NaN
  fixes$LonMove[I]=fixes$Lon[I]-fixes$Lon[I-1] 
  fixes$LonMove[abs(fixes$LonMove>50)]=NaN

    numcols=dim(fixes)[2]

    fixes$BigTas<-fixes$GAB<-fixes$SEIO<-fixes$Aust<-0
    I<-which(fixes$Lon>=110 & fixes$Lon<=155 & fixes$Lat<(-10) & fixes$Lat>=-45)
    fixes$Aust[I]<-1
    I<-which(fixes$Lon>=120 & fixes$Lon<=145 & fixes$Lat<=(-27.5) & fixes$Lat>=-42.5)
    fixes$GAB[I]<-1
    I<-which(fixes$Lon>=90 & fixes$Lon<=110 & fixes$Lat<=(-27.5) & fixes$Lat>=-42.5)
    fixes$SEIO[I]<-1
    I<-which(fixes$Lon>=150 & fixes$Lon<=170 & fixes$Lat<(-27.5) & fixes$Lat>=-42.5)
    fixes$BigTas[I]<-1

    x<-rle(fixes$ID)
    events<-data.frame(ID=x$values,Length=x$lengths,Date1=rep(0,length(x$values)),Move=rep(0,length(x$values)),Closed=rep(0,length(x$values)))
    enumcols=dim(events)[2]
    events$BigTas<-events$GAB<-events$SEIO<-events$Aust<-0

    for(i in 1:length(events$ID))
    {
    events$Date1[i]=min(fixes$Date[fixes$ID==events$ID[i]])
    I=which(fixes$ID==events[i,1])
#    events$Move[i]=spDistsN1(cbind(fixes$Lon[min(I)],fixes$Lat[min(I)]),
#                           cbind(fixes$Lon[max(I)],fixes$Lat[max(I)]),longlat=T)

    events$Closed[i]=max(fixes$Closed[I])
    events[i,(enumcols+1):(enumcols+length(regnames))]=apply(fixes[I,(numcols+1):(numcols+length(regnames))],2,max)
    }
    events$Month=floor(events$Date1/100)%%100
    events$MSLP=aggregate(fixes$MSLP,by=list(fixes$ID),FUN=max)[,2]
    events$CV=aggregate(fixes$CV,by=list(fixes$ID),FUN=max)[,2]
    events$Lat=aggregate(fixes$Lat,by=list(fixes$ID),FUN=mean)[,2]

#    events=events[events$Move>=500 & events$Closed==1,]
#    include<-match(fixes[,1],events[,1])
#    J<-which(is.na(include)==0)
#    fixes=fixes[J,]
  
  for(m in 1:12)
     for(tt in 1:length(regnames))
    {
      J=which(fixes$Month==m & fixes[,regnames[tt]]==1)
      I=which(events$Month==m & events[,regnames[tt]]==1)

      stats[y,m,1,tt,n]=length(I)
      stats[y,m,5,tt,n]=length(J)
      if(length(I)>0)
      {
        stats[y,m,2,tt,n]=length(which(events$MSLP[I]>=1020))
        stats[y,m,3,tt,n]=mean(events$Length[I])
        stats[y,m,4,tt,n]=mean(events$Move[I])
      }
      if(length(J)>0)
      {
        stats[y,m,6,tt,n]=mean(fixes$MSLP[J])
        stats[y,m,7,tt,n]=mean(fixes$CV[J])
        stats[y,m,8,tt,n]=mean(fixes$Radius[J],na.rm=T)
        stats[y,m,9,tt,n]=mean(fixes$Depth[J],na.rm=T)
        stats[y,m,10,tt,n]=mean(fixes$Move[J],na.rm=T)
        stats[y,m,11,tt,n]=mean(fixes$LonMove[J],na.rm=T)
        stats[y,m,12,tt,n]=mean(fixes$Lat[J],na.rm=T)
      }

    }
    
  }
}

##In case netcdf doesn't work
save(stats,file=paste(dir,outfile,sep=""))

