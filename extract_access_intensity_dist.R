## This is code that takes the 11 ensemble members of a POAMA run
## Started on a particular day
## And outputs some interested results for further analysis

## I care about ECLs, so main result is a daily timeseries of # of members that
## have an ECL on that day

## Also outputs a figure for each of the months in the file
## That shows the number of cyclones centred in that month
## Option to set a higher intensity threshold

library(ncdf4)

extract_counts<-function(year1,year2,inday=1,numdays=40,dir="../gcyc_out/access-s1/proj240_rad2cv1/",cvthresh=seq(0.15,3,0.05),outfile="ACCESS_austcyclones_40daylead",closed=T,lonlim=c(0,360),latlim=c(-90,90),outdir=".")
{
years=seq(year1,year2,1)
months=1:12
members=paste("e",sprintf(1:11,fmt="%2.2d"),sep="")
freqdist<-array(0,c(length(cvthresh),12,numdays))
dimnames(freqdist)[[1]]=cvthresh
dimnames(freqdist)[[2]]=month.name
dimnames(freqdist)[[3]]=paste0("lead.",seq(1,numdays))
cvfact=1/(cvthresh[2]-cvthresh[1])

for(y in 1:length(years))
for(m in 1:length(months))
for(e in 1:length(members))
{
indate=paste(years[y],sprintf("%02d",months[m]),sprintf("%02d",inday),sep="")
fname=paste(dir,members[e],"/tracks_",indate,"_4month.dat",sep="")
print(fname)
read.table(fname, sep="",skip=1)->fixes
colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Meh")
fixes$Date=fixes$Date+10000*(years[y]-floor(fixes$Date[1]/10000))
fixes$Lon=fixes$Lon%%360
fixes$CV=abs(fixes$CV)
fixes$CV2=floor(cvfact*fixes$CV)/cvfact

dates2=format.Date(seq.Date(as.Date(indate,format="%Y%m%d"),by="1 day",length.out=numdays),format="%Y%m%d")

if(closed) fixes=fixes[(fixes$Open==0 | fixes$Open==10),] # Remove all open systems
I=which(fixes$Date%in%dates2 & fixes$Lat>=latlim[1] & fixes$Lat<=latlim[2] & 
        fixes$Lon>=lonlim[1] & fixes$Lon<=lonlim[2])

tmp=table(factor(fixes$CV2[I],levels=cvthresh),factor(fixes$Date[I],levels=dates2))
freqdist[,m,]=freqdist[,m,]+tmp

}

## Next step - write ECL file
save(freqdist,file=paste0(outdir,outfile,".RData"))

} # End function


extract_counts(1990,2012,inday=1,numdays=40,cvthresh=seq(0.15,3,0.05),closed=T,
    dir="/short/eg3/asp561/cts.dir/gcyc_out/access-s1/proj240_lows_rad5cv0.15/",
    outdir="/short/eg3/asp561/cts.dir/gcyc_out/",lonlim=c(0,360),latlim=c(-90,90),
    outfile="ACCESS_globalcyclones_CVdist_proj240_rad5cv0.15_40daylead")

extract_counts(1990,2012,inday=1,numdays=40,cvthresh=seq(0.15,3,0.05),closed=F,
    dir="/short/eg3/asp561/cts.dir/gcyc_out/access-s1/proj240_lows_rad5cv0.15/",
    outdir="/short/eg3/asp561/cts.dir/gcyc_out/",lonlim=c(0,360),latlim=c(-90,90),
    outfile="ACCESS_globalcyclones_CVdist_proj240_rad5cv0.15_40daylead_open")
