## This is code that takes the 11 ensemble members of a POAMA run
## Started on a particular day
## And outputs some interested results for further analysis

## I care about ECLs, so main result is a daily timeseries of # of members that
## have an ECL on that day

## Also outputs a figure for each of the months in the file
## That shows the number of cyclones centred in that month
## Option to set a higher intensity threshold
library(sp)
library(ncdf4)
extract_counts<-function(year1,year2,type="low",dir="/short/eg3/asp561/cts.dir/gcyc_out",cvthresh=seq(0.15,3,0.05),dur=NA,outf=NA,outdir=dir,move=NA,closed=T,lonlim=c(0,360),latlim=c(-90,90))
{
years=seq(year1,year2,1)
months=1:12
freqdist<-array(0,c(length(cvthresh),12))
rownames(freqdist)=cvthresh
colnames(freqdist)=month.name
cvfact=1/(cvthresh[2]-cvthresh[1])

for(y in 1:length(years))
{
print(years[y])
fname=paste(dir,"/tracks_",years[y],".dat",sep="")
read.table(fname, sep="",skip=1)->fixes
colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Meh")
fixes$Year=floor(fixes$Date/10000)
if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
fixes$Month=floor(fixes$Date/100)%%100
fixes$Lon=fixes$Lon%%360
fixes$CV=abs(fixes$CV)
fixes$CV2=floor(cvfact*fixes$CV)/cvfact

### Make table of events to combine with DJF for exclusion
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
if(closed) fixes=fixes[(fixes$Open==0 | fixes$Open==10),] # Remove all open systems

I=which(fixes$Lat>=latlim[1] & fixes$Lat<=latlim[2] & 
        fixes$Lon>=lonlim[1] & fixes$Lon<=lonlim[2])
tmp=table(factor(fixes$CV2[I],levels=cvthresh),fixes$Month[I])
freqdist=freqdist+tmp
} # End year loop

save(freqdist,file=paste0(outdir,outf,".RData"))

} # End function

extract_counts(1990,2012,type="low",cvthresh=seq(0.15,3,0.05),closed=T,
     dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/daily_proj240_lows_rad5cv0.15/",
     outdir="/short/eg3/asp561/cts.dir/gcyc_out/",lonlim=c(0,360),latlim=c(-90,90),
     outf="ERAIdaily_UM_globalcyclones_CVdist_proj240_rad5cv0.15")

extract_counts(1990,2012,type="low",cvthresh=seq(0.15,3,0.05),closed=F,
     dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/daily_proj240_lows_rad5cv0.15/",
     outdir="/short/eg3/asp561/cts.dir/gcyc_out/",lonlim=c(0,360),latlim=c(-90,90),
     outf="ERAIdaily_UM_globalcyclones_CVdist_proj240_rad5cv0.15_open")

