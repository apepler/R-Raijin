library(ncdf4)
library(abind)
library(raster)

spreadeffect<-function(data,winwid=1,w2=matrix(1,3,3))
{
  a=dim(data)  
  if(lat[2]<lat[1]) ll=a[2]:1 else ll=1:a[2]  
  m1=abind(data[(a[1]-winwid+1):a[1],ll],data[,ll],data[1:winwid,ll],along=1)  
  m2=raster(t(m1),xmn=min(lon)-winwid,xmx=max(lon)+winwid,ymn=min(lat),ymx=max(lat))
  m3=focal(m2,w2,pad=T,padValue=0,fun=max)
  tmp=t(as.matrix(m3)[ll,(winwid+1):(a[1]+winwid)])
  return(tmp)
}

dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/"
a=nc_open(paste0(dir,"ACCESS_globalcyclones_proj240_rad5cv0.15_40daylead_strongestcyc_allleads.nc"))
lat=ncvar_get(a,"lat")
lon=ncvar_get(a,"lon")
years=1990:2012
months=1:12
indays=c(1,9,17,25)

wmatrix=matrix(1,3,3) # Max of surrounding five cells
#wmatrix=matrix(0,3,3)
#wmatrix[2:3,2:3]=1 # This version is only the centre cell and those to east and north
                     # Which means a 10x10 array, where the 0-5 and 5-10 cells combine & lon[1]=lon[1]+2.5
#lon=lon+2.5
#lat=lat+2.5


thresh=0.25
tmp=ncvar_get(a,"cyclones")
tmp2=abind(apply(tmp[,,,,,1:7,],c(1,2,3,4,5,7),max),apply(tmp[,,,,,8:14,],c(1,2,3,4,5,7),max),apply(tmp[,,,,,15:21,],c(1,2,3,4,5,7),max),
           apply(tmp[,,,,,22:28,],c(1,2,3,4,5,7),max),apply(tmp[,,,,,29:35,],c(1,2,3,4,5,7),max),along=7)

## Convert to a 15x15 array
access_weekly=array(0,dim(tmp2)[c(1,2,3,4,5,7)])
for(y in 1:length(years))
 for(m in 1:12)
  for(w in 1:5)
   for(i in 1:length(indays))
   for(e in 1:11)
     access_weekly[,,y,m,i,w]=access_weekly[,,y,m,i,w]+1*(spreadeffect(tmp2[,,y,m,i,e,w],w2=wmatrix)>=thresh)

access_weekly=access_weekly/11

a=nc_open(paste0(dir,"ERAIdaily_globalcyclones_proj240_rad5cv0.15_40daylead_strongestcyc_allleads.nc"))
tmp=ncvar_get(a,"cyclones")
tmp2=abind(apply(tmp[,,,,,1:7],c(1,2,3,4,5),max),apply(tmp[,,,,,8:14],c(1,2,3,4,5),max),apply(tmp[,,,,,15:21],c(1,2,3,4,5),max),
                     apply(tmp[,,,,,22:28],c(1,2,3,4,5),max),apply(tmp[,,,,,29:35],c(1,2,3,4,5),max),along=6)

control_weekly=array(0,dim(tmp2))

for(y in 1:length(years))
 for(m in 1:12)
  for(i in 1:length(indays))
   for(w in 1:5)
    control_weekly[,,y,m,i,w]=1*(spreadeffect(tmp2[,,y,m,i,w],w2=wmatrix)>=thresh)
   
dimnames(access_weekly)<-dimnames(control_weekly)<-list(lon,lat,years,1:12,indays,paste0("week.",1:5))

## Now we have two matrices by lon/lat/year/month
## Next step: save as .RData so I can manipulate more easily on desktop

save(lat,lon,years,months,thresh,access_weekly,indays,control_weekly,file=paste0(dir,"ACCESS_weeklyeval_cv",thresh,"_15deg_allleads.RData")) 



