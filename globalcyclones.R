##This is code to take the output from the low tracking software 
##and turn it into a .csv file of ECL track data.

##Note that this script performs only the initial filtering 
##to restrict ECLs based on location, track length and at least one closed fix
##To obtain a better database of ECLs with impacts you can filter for 
##a fix curvature of >= 0.25 & the fix being in the ECL region

##This version of the code requires both the annual version of the data
##as well as that run for the previous summer, to collate the tracks for events
##over January 1.

##Note also that R has issues with large datasets - 
##I tend to separate into sets of 20 years or fewer. Pretty quick. 

##yearS and yearE are start & end years. Must have DJF data for yearS-1
##output is a string to help label the fix & event files (so use "name")

##You will need to set the location to where the track files are stored.
library(sp)
get_fixes<-function(year,output,dir,move=NA,dur=NA,cv=NA)
{
fname=paste(dir,"/tracks_",year,".dat",sep="")
read.table(fname, sep="",skip=1)->fixes
colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Depth","Radius","Up","Vp")
tmp=floor(fixes$Date/10000)
if(length(unique(tmp))>1) fixes=fixes[tmp==unique(tmp)[2],]
fixes$Date=year*10000+fixes$Date%%10000
fixes$Lon=fixes$Lon%%360
fixes$Year=year
fixes$Month=floor(fixes$Date/100)%%100
fixes$CV=abs(fixes$CV)

 if(!is.na(move))
 {
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

  events=events[events$Move>=move,]
  include<-match(fixes[,1],events[,1])
  J<-which(is.na(include)==0)
  fixes=fixes[J,]
 }

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
 
 write.csv(fixes,file=paste0(dir,output))

}

for(year in 1980:2016)
  get_fixes(year,move=500,
  out=paste0("ERAI_UM_highs_ERAI_proj100_rad10cv0.075_500km_",year,".csv"),
  dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_highs_rad10cv0.075/")

