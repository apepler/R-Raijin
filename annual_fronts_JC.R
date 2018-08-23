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
library(abind)

spreadeffect<-function(data,winwid=3,circ=T,lon=NaN,lat=NaN)
{
  a=dim(data)

  if(is.na(lon[1]))
  {
    lon=seq(1,a[1])
    lat=seq(1,a[2])
  }

  if(lat[2]<lat[1]) ll=a[2]:1 else ll=1:a[2]

  ## Only make cyclic if it's actually global, -180 to 180 or 0 to 360
  if(lon[2]-lon[1]>=350) m1=abind(data[(a[1]-winwid+1):a[1],ll],data[,ll],data[1:winwid,ll],along=1) else
    m1=abind(matrix(0,winwid,a[2]),data[,ll],matrix(0,winwid,a[2]),along=1)

  m2=raster(t(m1),xmn=min(lon)-winwid,xmx=max(lon)+winwid,ymn=min(lat),ymx=max(lat))
  if(circ) w2=focalWeight(m2,winwid,type="circle") else w2=matrix(1,1+2*winwid,1+2*winwid)
  m3=focal(m2,w2,pad=T,padValue=0)
  tmp=t(as.matrix(m3)[ll,(winwid+1):(a[1]+winwid)])
  tmp[tmp>0]=1
  return(tmp)
}


jen_compare<-function(year1,year2,dir="/short/eg3/asp561/cts.dir/gcyc_out/",outf=NA,inf=NA)
{

years=seq(year1,year2)
a=nc_open(paste0(dir,inf,years[1],".nc"))
lon=ncvar_get(a,"lon")
lat=ncvar_get(a,"lat")

def_f<-cat_f<-rud_f<-array(0,c(length(lon),length(lat),length(years),12))

# Load front data

for(y in 1:length(years))
{
datelist=seq.POSIXt(as.POSIXct(paste0(years[y],"0101 00:00"),format="%Y%m%d %H:%M",tz="GMT"),as.POSIXct(paste0(years[y]+1,"0101 00:00"),format="%Y%m%d %H:%M",tz="GMT"),by="6 hours")
datelist=datelist[-length(datelist)]
monthlist=as.numeric(format(datelist,"%m"))

a=nc_open(paste0(dir,inf,years[y],".nc"))
fronts=ncvar_get(a,"FRONTFLAG")

for(m in 1:12)
{
I=which(monthlist==m)
tmp=fronts[,,I]

def_f[,,y,m]=apply(tmp==1,c(1,2),sum)
cat_f[,,y,m]=apply(tmp==2,c(1,2),sum)
rud_f[,,y,m]=apply(tmp==3,c(1,2),sum)
}
}

dimX<-ncdim_def("lon","degrees_E",lon)
dimY<-ncdim_def("lat","degrees_N",lat)
dimT1<-ncdim_def("year","years",years)
dimT2<-ncdim_def("month","months",1:12)

fillvalue <- -9999
def <- ncvar_def("comb_fronts","count",list(dimX,dimY,dimT1,dimT2),fillvalue,"Hours with a cold front using both Catto & Rudeva methods",prec="integer")
def2 <- ncvar_def("catto_fronts","count",list(dimX,dimY,dimT1,dimT2),fillvalue,"Hours with a cold front in Catto data only",prec="integer")
def3 <- ncvar_def("rudeva_fronts","count",list(dimX,dimY,dimT1,dimT2),fillvalue,"Hours with a cold front in Rudeva data only",prec="integer")


# create netCDF file and put arrays
ncout <- nc_create(paste0(dir,outf),list(def,def2,def3)) #force_v4=T)

# put variables
ncvar_put(ncout,def,def_f)
ncvar_put(ncout,def2,cat_f)
ncvar_put(ncout,def3,rud_f)

# put additional attributes into dimension and data variables
ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout,"lat","axis","Y")
ncatt_put(ncout,"year","axis","T1")
ncatt_put(ncout,"month","axis","T2")

nc_close(ncout)

} # End function

jen_compare(2005,2015,
     dir="/g/data/eg3/asp561/Fronts/",
     inf="frontcomp_cattocold_rudevalen3_plus6h_3deg_",
     outf="frontcomp_cattocold_rudevalen3_plus6h_3deg_20052015.nc")


