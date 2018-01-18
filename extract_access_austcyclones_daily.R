## This is code that takes the 11 ensemble members of a POAMA run
## Started on a particular day
## And outputs some interested results for further analysis

## I care about ECLs, so main result is a daily timeseries of # of members that
## have an ECL on that day

## Also outputs a figure for each of the months in the file
## That shows the number of cyclones centred in that month
## Option to set a higher intensity threshold

library(ncdf4)

ACCESS_ECLs_4month<-function(year1,year2,inday=1,numdays=40,dir="../gcyc_out/access-s1/proj240_rad2cv1/",thresh=1,outfile="ACCESS_austcyclones_40daylead.nc",closed=T)
{
years=seq(year1,year2,1)
months=1:12
members=paste("e",sprintf(0:10,fmt="%2.2d"),sep="")

lat=seq(-89.5,89.5,1)
lon=seq(0.5,359.5,1)
members=paste("e",sprintf(1:11,fmt="%2.2d"),sep="")

cyclones<-array(0,c(length(lon),length(lat),length(years),length(months),numdays))

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
fixes$Lat2=floor(fixes$Lat)
fixes$Lon2=floor(fixes$Lon)%%360
dates2=format.Date(seq.Date(as.Date(indate,format="%Y%m%d"),by="1 day",length.out=numdays),format="%Y%m%d")

if(closed) fixes=fixes[(fixes$Open==0 | fixes$Open==10),] # Remove all open systems
I=which(fixes$Date%in%dates2 & fixes$CV>=thresh &
        fixes$Lat2%in%floor(lat) & fixes$Lon2%in%floor(lon))

tmp=table(factor(fixes$Lon2[I],levels=floor(lon)),factor(fixes$Lat2[I],levels=floor(lat)),factor(fixes$Date[I],levels=dates2))
cyclones[,,y,m,]=cyclones[,,y,m,]+tmp # So ends up number of cyclones on that day across all members
}

## Next step - write ECL file

dimT1<-ncdim_def("year","years",years)
dimT2<-ncdim_def("month","months",months)
dimM<-ncdim_def("leadtime","days",0:(numdays-1))

fillvalue <- 1e32

## Write a netcdf file with cyclones - so can read later

dimX<-ncdim_def("lon","degrees_E",lon)
dimY<-ncdim_def("lat","degrees_N",lat)
cyc_def <- ncvar_def("cyclones","count",list(dimX,dimY,dimT1,dimT2,dimM),fillvalue,"Number of cyclones identified per day in the Australian region",prec="single")

# create netCDF file and put arrays
ncfname <- paste(dir,outfile,sep="")
ncout <- nc_create(ncfname,cyc_def,force_v4=T)

# put variables
ncvar_put(ncout,cyc_def,cyclones)

# put additional attributes into dimension and data variables
ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout,"lat","axis","Y")
ncatt_put(ncout,"year","axis","T1")
ncatt_put(ncout,"month","axis","T2")
ncatt_put(ncout,"leadtime","axis","M")

nc_close(ncout)


} # End function


ACCESS_ECLs_4month(1990,2012,inday=1,numdays=40,dir="/short/eg3/asp561/cts.dir/gcyc_out/access-s1/proj240_lows_rad5cv0.15/",thresh=0.25,outfile="ACCESS_globalcyclones_proj240_rad5cv0.25_40daylead.nc",closed=T)
