library(ncdf4)
library(raster)
library(abind)
library(sp)

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

expandgrid<-function(dir,file1,file2,sname,widtype="deg",winwid=3,wincirc=T,type="system")
{

a=nc_open(paste0(dir,file1))
data1=ncvar_get(a,sname)
lat=ncvar_get(a,"lat")
lon=ncvar_get(a,"lon")
datelist=as.POSIXlt(ncvar_get(a,"time")*6*60*60,origin=paste0(year,"-01-01 00:00"),tz="UTC")

data2=array(0,dim(data1))
for(i in 1:length(datelist)) data2[,,i]=spreadeffect(data1[,,i],lon=lon,lat=lat,circ=wincirc,winwid=winwid)

## Save as a netcdf

dimX<-ncdim_def("lon","degrees_E",lon)
dimY<-ncdim_def("lat","degrees_N",lat)
dimT<-ncdim_def("time","hours since 1970-1-1 00:00:00",as.numeric(datelist)/(60*60))

fillvalue <- 1e32
cyc_def <- ncvar_def("systems","count",list(dimX,dimY,dimT),fillvalue,paste("Cells affected by a ",type),prec="single")

# create netCDF file and put arrays
ncout <- nc_create(paste0(dir,file2),cyc_def) #force_v4=T)

# put variables
ncvar_put(ncout,cyc_def,data2)

# put additional attributes into dimension and data variables
ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout,"lat","axis","Y")
ncatt_put(ncout,"time","axis","T")

nc_close(ncout)
}

for(year in 2005:2015) expandgrid("/g/data/eg3/asp561/CattoData_20052015/",paste0("cold_fronts_Aus",year,".nc"),paste0("cold_fronts_Aus",year,"_v2.nc"),sname="Cold fronts for use in collaboration with A. Dowdy & J. Catto:",type="cold front")


