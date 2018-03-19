## This is code that takes the 11 ensemble members of a POAMA run
## Started on a particular day
## And outputs some interested results for further analysis

## I care about ECLs, so main result is a daily timeseries of # of members that
## have an ECL on that day

## Also outputs a figure for each of the months in the file
## That shows the number of cyclones centred in that month
## Option to set a higher intensity threshold

library(ncdf4)

ACCESS_ECLs_4month<-function(year1,year2,indays=1,numdays=40,dir="../gcyc_out/access-s1/proj240_rad2cv1/",thresh=1,outfile="ACCESS_austcyclones_40daylead.nc",closed=T,outdir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/")
{
years=seq(year1,year2,1)
months=1:12
members=paste("e",sprintf(1:11,fmt="%2.2d"),sep="")

lat=seq(-87.5,87.5,5)
lon=seq(2.5,357.5,5)

cyclones<-array(0,c(length(lon),length(lat),length(years),length(months),length(indays),numdays,length(members)))

for(y in 1:length(years))
for(m in 1:length(months))
for(ii in 1:length(indays))
for(e in 1:length(members))
{
indate=paste(years[y],sprintf("%02d",months[m]),sprintf("%02d",indays[ii]),sep="")
fname=paste(dir,members[e],"/tracks_",indate,"_4month.dat",sep="")
print(fname)
read.table(fname, sep="",skip=1)->fixes
colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Meh")
fixes$Date=fixes$Date+10000*(years[y]-floor(fixes$Date[1]/10000))
fixes$Lat2=5*floor(fixes$Lat/5)
fixes$Lon2=5*floor((fixes$Lon%%360)/5)
dates2=format.Date(seq.Date(as.Date(indate,format="%Y%m%d"),by="1 day",length.out=numdays),format="%Y%m%d")

if(closed) fixes=fixes[(fixes$Open==0 | fixes$Open==10),] # Remove all open systems
I=which(fixes$Date%in%dates2)

tmp2=aggregate(fixes$CV[I],by=list(fixes$Lon2[I],fixes$Lat2[I],fixes$Date[I]),FUN=max)
jlat=match(tmp2[,2],seq(-90,85,5))
ilon=match(tmp2[,1],seq(0,355,5))
idate=match(tmp2[,3],dates2)
for(nn in which(!is.na(jlat))) cyclones[ilon[nn],jlat[nn],y,m,ii,idate[nn],e]=as.numeric(tmp2[nn,4])
}

## Next step - write ECL file

dimT1<-ncdim_def("year","years",years)
dimT2<-ncdim_def("month","months",months)
dimT3<-ncdim_def("startday","days",indays)
dimM<-ncdim_def("leadtime","days",0:(numdays-1))
dimE<-ncdim_def("member","counts",1:11)

fillvalue <- 1e32

## Write a netcdf file with cyclones - so can read later

dimX<-ncdim_def("lon","degrees_E",lon)
dimY<-ncdim_def("lat","degrees_N",lat)
cyc_def <- ncvar_def("cyclones","count",list(dimX,dimY,dimT1,dimT2,dimT3,dimM,dimE),fillvalue,"Strongest cyclone identified per day in the Australian region",prec="single")

# create netCDF file and put arrays
ncfname <- paste(outdir,outfile,sep="")
ncout <- nc_create(ncfname,cyc_def,force_v4=T)

# put variables
ncvar_put(ncout,cyc_def,cyclones)

# put additional attributes into dimension and data variables
ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout,"lat","axis","Y")
ncatt_put(ncout,"year","axis","T1")
ncatt_put(ncout,"month","axis","T2")
ncatt_put(ncout,"startday","axis","T3")
ncatt_put(ncout,"leadtime","axis","M")
ncatt_put(ncout,"member","axis","E")

nc_close(ncout)


} # End function


ACCESS_ECLs_4month(1990,2012,indays=c(1,9,17,25),numdays=40,dir="/short/eg3/asp561/cts.dir/gcyc_out/access-s1/proj240_lows_rad5cv0.15/",outfile="ACCESS_globalcyclones_proj240_rad5cv0.15_40daylead_strongestcyc_allleads.nc",closed=T)
