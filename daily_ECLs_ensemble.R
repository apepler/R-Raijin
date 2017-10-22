## This is code that takes the 11 ensemble members of a POAMA run
## Started on a particular day
## And outputs some interested results for further analysis

## I care about ECLs, so main result is a daily timeseries of # of members that
## have an ECL on that day

## Also outputs a figure for each of the months in the file
## That shows the number of cyclones centred in that month
## Option to set a higher intensity threshold

library(ncdf4)
daily_ECLs<-function(year1,year2,inday,type="low",dir="/short/eg3/asp561/cts.dir/gcyc_out",cvthresh=c(0.15,0.2,0.25,0.3),dur=NA,outfile=NA)
{
years=seq(year1,year2,1)
months=1:12
numdays=125 ## Number of days the scheme was run for
members=paste("e",sprintf(1:11,fmt="%2.2d"),sep="")

ECLlist=array(0,c(length(years),length(months),length(members),numdays,length(cvthresh)))

outdates=array(0,c(length(years),length(months),numdays))

for(y in 1:length(years))
for(m in 1:length(months))
if(m!=4) # Bcause still missing
for(e in 1:length(members))
{
print(years[y])
indate=paste(years[y],sprintf("%02d",months[m]),sprintf("%02d",inday),sep="")
fname=paste(dir,members[e],"/tracks_",indate,"_4month.dat",sep="")
read.table(fname, sep="",skip=1)->fixes
colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Meh")
fixes$Year=floor(fixes$Date/10000)
  
yy2=unique(fixes$Year)
if(length(yy2)>1)
  {
    I=which(fixes$Year==yy2[1])
    fixes$Date[I]=(fixes$Date[I]%%10000) + years[y]*10000
    fixes$Date[-I]=(fixes$Date[-I]%%10000) + (years[y]+1)*10000
  } else fixes$Date=(fixes$Date%%10000) + years[y]*10000

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

datelist=as.numeric(format(seq.Date(as.Date(indate,format="%Y%m%d"),by="1 day",length.out=numdays),"%Y%m%d"))

if(e==1) outdates[y,m,]=datelist

for(t in 1:length(cvthresh))
{
J=which(fixes$Location==1 & fixes$CV>=cvthresh[t])
tmp=table(factor(fixes$Date[J],levels=datelist))
ECLlist[y,m,e,,t]=tmp
}
} # End year loop

dimT1<-ncdim_def("year","years",years)
dimT2<-ncdim_def("month","months",months)
dimN<-ncdim_def("members","count",1:length(members))
dimM<-ncdim_def("leadtime","days",0:(numdays-1))
dimC<-ncdim_def("threshold","lap",cvthresh)

fillvalue <- 1e32

## Write a netcdf file with cyclones - so can read later

cyc_def <- ncvar_def("ECLs","count",list(dimT1,dimT2,dimN,dimM,dimC),fillvalue,"Number of ECLs counted at different leads",prec="single")
day_def <- ncvar_def("Day","daynum",list(dimT1,dimT2,dimM),99999,"Date in YYMMDDDD",prec="integer")

# create netCDF file and put arrays
if(is.na(outfile)) outfile=paste(dir,"/ECLs_4monthlead_day",inday,".nc",sep="") else outfile=paste(dir,outfile,sep="/")
ncout <- nc_create(outfile,list(cyc_def,day_def),force_v4=T)

# put variables
ncvar_put(ncout,cyc_def,ECLlist)
ncvar_put(ncout,day_def,outdates)


} # End function

daily_ECLs(1990,2012,inday=1,dir="/short/eg3/asp561/cts.dir/gcyc_out/access-s1/proj240_rad5cv0.15/",outfile="ACCESS_ECLs_proj240_rad5cv0.15.nc",cvthresh=seq(0.15,0.5,0.05))

