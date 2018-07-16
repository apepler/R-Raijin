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


jen_compare<-function(year1,year2,dir="/short/eg3/asp561/cts.dir/gcyc_out/",outf=NA,inf=NA,winwid=0,latlim=c(-45,-10),lonlim=c(110,155))
{

years=seq(year1,year2)
a=nc_open(paste0(dir,inf,years[1],".nc"))
lon=ncvar_get(a,"lon")
lat=ncvar_get(a,"lat")

J=which(lat>=min(latlim) & lat<=max(latlim))
if(max(lonlim)<=180) I=which(lon>=min(lonlim) & lon<=max(lonlim)) else
{
if(min(lonlim)<180) I=which(lon>=min(lonlim) | lon<=(max(lonlim)-360)) else
I=which(lon>=(min(lonlim)-360) & lon<=(max(lonlim)-360))
}


# First, load Jen's cyclones

for(y in 1:length(years))
{
datelist=seq.POSIXt(as.POSIXct(paste0(years[y],"0101 00:00"),format="%Y%m%d %H:%M",tz="GMT"),as.POSIXct(paste0(years[y]+1,"0101 00:00"),format="%Y%m%d %H:%M",tz="GMT"),by="6 hours")
datelist=datelist[-length(datelist)]

a=nc_open(paste0(dir,inf,years[y],".nc"))
jc_cyc=ncvar_get(a,"CYCFLAG")

if(winwid!=0) # If we want to spread Jen cyclones out again
{
for(i in 1:length(datelist)) jc_cyc[,,i]=spreadeffect(jc_cyc[,,i],lon=lon,lat=lat,circ=T,winwid=winwid)
}

tmp=data.frame(Date=format(datelist,"%Y%m%d"),Time=format(datelist,"%H"),
                  Region.DefCyc=apply(jc_cyc[I,J,]==1,3,sum),
                  Region.MaybeCyc=apply(jc_cyc[I,J,]==2,3,sum))

if(y==1) cyc_out=tmp else cyc_out=rbind(cyc_out,tmp)
}

cyc_out[,3:4]=cyc_out[,3:4]/(length(I)*length(J))
write.csv(cyc_out,file=paste0(dir,outf))

} # End function

jen_compare(1980,2015,winwid=0,
     dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_lows_rad5cv0.15/",
     inf="Catto_cyclones_proj100_rad5cv0.15_",
     outf="Catto_viccyclones_proj100_rad5cv0.15_hourlyprop.csv",
     latlim=c(-40,-30),lonlim=c(135,155))


