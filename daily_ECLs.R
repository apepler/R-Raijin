## This is code that takes the 11 ensemble members of a POAMA run
## Started on a particular day
## And outputs some interested results for further analysis

## I care about ECLs, so main result is a daily timeseries of # of members that
## have an ECL on that day

## Also outputs a figure for each of the months in the file
## That shows the number of cyclones centred in that month
## Option to set a higher intensity threshold

library(ncdf4)
daily_ECLs<-function(year1,year2,type="low",dir="/short/eg3/asp561/cts.dir/gcyc_out",cvthresh=c(0.15,0.2,0.25,0.3),dur=NA,outfile=NA)
{
years=seq(year1,year2,1)
datelist=as.numeric(format(seq.Date(as.Date(paste(year1,"-01-01",sep="")),
    as.Date(paste(year2,"-12-31",sep="")),by="1 day"),"%Y%m%d"))
yearlist=floor(datelist/10000)

ECLlist=as.data.frame(cbind(datelist,matrix(0,length(datelist),length(cvthresh))))
colnames(ECLlist)=c("Date",paste("CV",cvthresh,sep=""))

for(y in 1:length(years))
{
print(years[y])
fname=paste(dir,"/tracks_",years[y],".dat",sep="")
read.table(fname, sep="",skip=1)->fixes
colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Meh")
fixes$Year=floor(fixes$Date/10000)
fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
fixes$Date=years[y]*10000+fixes$Date%%10000
if(type=="high") fixes$CV=-fixes$CV

fixes$Location=0
  I<-which(fixes$Lon>=149 & fixes$Lon<=161 & fixes$Lat<(-37) & fixes$Lat>=-41)
  fixes$Location[I]<-1
  I<-which(fixes$Lon>=(149+(37+fixes$Lat)/2) & fixes$Lon<=161 & fixes$Lat<(-31) & fixes$Lat>=-37)
  fixes$Location[I]<-1
  I<-which(fixes$Lon>=152 & fixes$Lon<=161 & fixes$Lat<=(-24) & fixes$Lat>=-31)
  fixes$Location[I]<-1

 if(!is.na(dur))
  {
  x<-rle(fixes$ID)
  events<-cbind(x$values,x$lengths,matrix(data=0,nrow=length(x$values),ncol=1))
  events=events[events[,2]>=dur,]
  include<-match(fixes[,1],events[,1])
  J<-which(is.na(include)==0)
  fixes=fixes[J,]
  }

I=which(yearlist==years[y])

for(t in 1:length(cvthresh))
{
J=which(fixes$Location==1 & fixes$CV>=cvthresh[t])
tmp=table(factor(fixes$Date[J],levels=datelist[I]))
ECLlist[I,t+1]=tmp
}
} # End year loop

if(is.na(outfile)) outfile="Daily_ECLs.csv"
outfile=paste(dir,outfile,sep="/")
write.csv(ECLlist,file=outfile)

} # End function

daily_ECLs(1990,2013,type="low",dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/daily_proj240_rad5cv0.15/",outfile="ERAI_ECLs_proj240_rad5cv0.15.csv",cvthresh=seq(0.15,0.5,0.05))

