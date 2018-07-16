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
extract_counts<-function(year1,year2,type="low",dir="/short/eg3/asp561/cts.dir/gcyc_out",thresh=0,dur=NA,outf=NA,outdir=dir,move=NA,closed=F)
{
years=seq(year1,year2,1)
months=1:12

lat=seq(-89.5,89.5,1)
lon=seq(0.5,359.5,1)  ### Can always combine into bigger cells later
systems<-array(0,c(length(lon),length(lat),length(years),12))

for(y in 1:length(years))
{
print(years[y])
fname=paste(dir,"/tracks_",years[y],".dat",sep="")
read.table(fname, sep="",skip=1)->fixes
colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Meh")
fixes$Year=floor(fixes$Date/10000)
if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
fixes$Month=floor(fixes$Date/100)%%100
fixes$Lat2=floor(fixes$Lat)
fixes$Lon2=floor(fixes$Lon)%%360
fixes$CV=abs(fixes$CV)

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
fixes=fixes[fixes$CV>=thresh,]
if(closed) fixes=fixes[(fixes$Open==0 | fixes$Open==10),] # Remove all open systems
tmp=table(factor(fixes$Lon2,levels=0:359),factor(fixes$Lat2,levels=-90:89),fixes$Month)
systems[,,y,]=tmp
print(mean(systems[,,y,],na.rm=T))
} # End year loop

## Write a netcdf file with cyclones - so can read later

#dimX<-dim.def.ncdf("lon","degrees_E",lon)
#dimY<-dim.def.ncdf("lat","degrees_N",lat)
#dimT1<-dim.def.ncdf("year","years",years)
#dimT2<-dim.def.ncdf("month","months",1:12)

dimX<-ncdim_def("lon","degrees_E",lon)
dimY<-ncdim_def("lat","degrees_N",lat)
dimT1<-ncdim_def("year","years",years)
dimT2<-ncdim_def("month","months",months)

fillvalue <- 1e32
cyc_def <- ncvar_def("systems","count",list(dimX,dimY,dimT1,dimT2),fillvalue,paste("Number of",type,"pressure systems during each month for each location in the Australian region"),prec="single")

# create netCDF file and put arrays
if(!is.na(outf)) ncfname <- paste(outdir,outf,sep="") else
  if(is.na(dur)) ncfname <- paste(outdir,"/count_",type,"s","_cv",thresh,".nc",sep="") else ncfname<-paste(outdir,"/count_",type,"s","_cv",thresh,"_D",dur,".nc",sep="")
ncout <- nc_create(ncfname,cyc_def) #force_v4=T)

# put variables
ncvar_put(ncout,cyc_def,systems)

# put additional attributes into dimension and data variables
ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout,"lat","axis","Y")
ncatt_put(ncout,"year","axis","T1")
ncatt_put(ncout,"month","axis","T2")

nc_close(ncout)

} # End function

#extract_counts(1980,2016,type="high",move=500,
#     dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_highs_rad10cv0.075_notopo/",
#     outdir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/",
#     outf="ERAI_UM_globalanticyclones_proj100_rad10cv0.75_notopo_500km.nc")

for(lev in 500)
{
extract_counts(1980,2016,type="low",thresh=4,dur=2,
     dir=paste0("/short/eg3/asp561/cts.dir/gcyc_out/ERAI/",lev,"hPa_z/proj100_lows_rad5cv1/"),
     outdir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/",
     outf=paste0("ERAI_UM_globalcyclones_",lev,"hPa_proj100_rad5cv4_D2.nc"))
}
extract_counts(1980,2016,type="low",thresh=1,dur=2,
     dir=paste0("/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_lows_rad2cv1/"),
     outdir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/",
     outf=paste0("ERAI_UM_globalcyclones_proj100_rad2cv1_D2.nc"))


#extract_counts(1980,2016,type="high",thresh=0.075,
#     dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_highs_rad10cv0.075/",
#     outdir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/",
#     outf="ERAI_UM_globalanticyclones_proj100_rad10cv0.075.nc")

#for(cv in c(0.25,0.5))
#{
#extract_counts(1980,2016,type="low",thresh=cv,closed=T,move=500,
#     dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_lows_rad5cv0.15/",
#     outdir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/",
#     outf=paste0("ERAI_UM_globalcyclones_proj100_rad5cv",cv,"_500km.nc"))
#}
#extract_counts(1950,2016,type="low",thresh=0.15,move=500,closed=T,
#     dir="/short/eg3/asp561/cts.dir/gcyc_out/NCEP1/proj100_lows_rad5cv0.15/",
#     outdir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/",
#     outf="NCEP1_UM_globalcyclones_proj100_rad5cv0.15_500km.nc")

#name="ACCESS1-3"
#basedir="/short/eg3/asp561/cts.dir/gcyc_out/CMIP5/"

#thresh="rad5cv0.15"
#for(thresh in c("rad2cv1","rad5cv0.15"))
#{
#extract_counts(1950,2005,type="low",move=500,
#    dir=paste0(basedir,name,"/historical/r1i1p1/proj100_lows_",thresh),
#    outdir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/",
#    outf=paste0(name,"_historical_r1i1p1_globalcyclones_proj100_",thresh,"_500km.nc"))

#extract_counts(2006,2100,type="low",move=500,
#    dir=paste0(basedir,name,"/rcp85/r1i1p1/proj100_lows_",thresh),
#    outdir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/",
#    outf=paste0(name,"_rcp85_r1i1p1_globalcyclones_proj100_",thresh,"_500km.nc"))
#}

