## This is code that takes the 11 ensemble members of a POAMA run
## Started on a particular day
## And outputs some interested results for further analysis

## I care about ECLs, so main result is a daily timeseries of # of members that
## have an ECL on that day

## Also outputs a figure for each of the months in the file
## That shows the number of cyclones centred in that month
## Option to set a higher intensity threshold
library(raster)
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


jen_compare<-function(year,dir="/short/eg3/asp561/cts.dir/gcyc_out/",jdir="/g/data/eg3/asp561/CattoData_20052015/",outf=NA,lonlim=c(-180,180),latlim=c(-90,90),winwid=0,winwidH=7,thresh=0,dur=NA,move=NA,closed=F,cold=T,daily=F)
{

datelist=seq.POSIXt(as.POSIXct(paste0(year,"0101 00:00"),format="%Y%m%d %H:%M",tz="GMT"),as.POSIXct(paste0(year+1,"0101 00:00"),format="%Y%m%d %H:%M",tz="GMT"),by="6 hours")
datelist=datelist[-length(datelist)]

# First, load Jen's data

if(cold)
{
a=nc_open(paste0(jdir,"cold_fronts_Aus",year,".nc"))
jc_front=ncvar_get(a,"Cold fronts for use in collaboration with A. Dowdy & J. Catto:")
lat=ncvar_get(a,"lat")
lon=ncvar_get(a,"lon")

if(winwid>0) 
{
tmp=jc_front
for(i in 1:length(datelist)) jc_front[,,i]=spreadeffect(tmp[,,i],lon=lon,lat=lat,circ=T,winwid=winwid)
}

} else {
a=nc_open(paste0(jdir,"storm_types_Aus333_",year,".nc"))
tmp=ncvar_get(a,"Storm types for use in collaboration with A. Dowdy & J. Catto:")
jc_front<-array(0,dim(tmp))
jc_front[tmp%in%c(2,4,6,7)]=1 ## Fronts
lat=ncvar_get(a,"lat")
lon=ncvar_get(a,"lon")
}

dlat=abs(lat[2]-lat[1])
dlon=abs(lon[2]-lon[1])

if(min(lat)>-80 | max(lat)<=80) # Need a different latitude for fronts to allow region of effect
{

rlon=c(rev(seq(min(lon)-dlon,min(lon)-10,-dlon)),
      lon,
      seq(max(lon)+dlon,max(lon)+10,dlon))
lon2=(rlon%%360)/dlon

if(lat[2]>lat[1])
{
rlat=c(rev(seq(min(lat)-dlat,min(lat)-10,-dlat)),
      lat,
      seq(max(lat)+dlat,max(lat)+10,dlat))
lat2=rlat/dlat
} else {
rlat=c(rev(seq(max(lat)+dlat,max(lat)+10,dlat)),
      lat,
      seq(min(lat)-dlat,min(lat)-10,-dlat))
lat2=rlat/dlat
}

} else 
{
rlat=lat
rlon=lon
lat2=lat/dlat
lon2=lon/dlat
}

final_front<-array(0,dim(jc_front))
fname=paste(dir,"/tracks_",year,".dat",sep="")
read.table(fname, sep="",skip=1)->fixes
colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Meh")
fixes$Year=floor(fixes$Date/10000)
if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
fixes$Month=floor(fixes$Date/100)%%100
fixes$CV=abs(fixes$CV)

fixes$Lat2=round(fixes$Lat/dlat,0)
fixes$Lon2=fixes$Lon%%360
if(min(lon)<0) fixes$Lon2[fixes$Lon2>180]=fixes$Lon2[fixes$Lon2>180]-360
fixes$Lon2=round(fixes$Lon2/dlon,0)
fixes$Date2=as.POSIXct(paste0(year*10000+fixes$Date%%10000," ",fixes$Time),format="%Y%m%d %H:%M",tz="GMT")

 if(!is.na(move))
 {
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
  fixes=fixes[abs(fixes$CV)>=thresh,]
  if(closed) fixes=fixes[(fixes$Open==0 | fixes$Open==10),] # Remove all open sys

grid1=table(factor(fixes$Lon2,levels=lon2),factor(fixes$Lat2,levels=lat2),factor(as.character(fixes$Date2),levels=as.character(datelist)))

highs=array(0,dim(grid1))
for(i in 1:length(datelist)) highs[,,i]=spreadeffect(grid1[,,i],lon=rlon,lat=rlat,circ=T,winwid=winwidH)

## Now, need to shrink this if necessary
if(min(lat)>-80 | max(lat)<=80)
{
I=which(rlon%in%lon)
J=which(rlat%in%lat)
highs=highs[I,J,]
}

I=which(jc_front==1 & highs==0)
final_front[I]=1
I=which(jc_front==1 & highs==1)
final_front[I]=2

if(daily)
{
datelist2=as.numeric(format(datelist,"%Y%m%d"))
daylist=unique(datelist2)
tmp=array(0,c(length(lon),length(lat),length(daylist)))

for(t in 1:length(daylist))
{
tmp1=apply(final_front[,,datelist2==daylist[t]]==1,c(1,2),sum)
tmp2=apply(final_front[,,datelist2==daylist[t]]==2,c(1,2),sum)

tmp3=array(0,c(length(lon),length(lat)))
tmp3[tmp2>0]=2
tmp3[tmp1>0]=1 # If any are highs-fre it's highs-free, otherwise not
tmp[,,t]=tmp3
}
final_front=tmp
}

I=which(lon>=min(lonlim) & lon<=max(lonlim))
J=which(lat>=min(latlim) & lat<=max(latlim))

dimX<-ncdim_def("lon","degrees_E",lon[I])
dimY<-ncdim_def("lat","degrees_N",lat[J])

if(daily) {
tmp=seq.POSIXt(as.POSIXct(paste0(year,"0101 09:00"),format="%Y%m%d %H:%M",tz="GMT"),as.POSIXct(paste0(year,"1231 09:00"),format="%Y%m%d %H:%M",tz="GMT"),by="24 hours")
dimT<-ncdim_def("time","hours since 1970-1-1 00:00:00",as.numeric(tmp)/(60*60))  
} else  dimT<-ncdim_def("time","hours since 1970-1-1 00:00:00",as.numeric(datelist)/(60*60))

fillvalue <- 9999
print(unique(as.numeric(final_front),na.rm=T))
cyc_def <- ncvar_def("FRONTFLAG","count",list(dimX,dimY,dimT),fillvalue,"Front flag - 1 if a front in Catto is not overlapped by a high, 2 for if it is (so not real)",prec="integer")

# create netCDF file and put arrays
ncout <- nc_create(paste0(dir,outf),cyc_def) #force_v4=T)

# put variables
ncvar_put(ncout,cyc_def,final_front[I,J,])

# put additional attributes into dimension and data variables
ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout,"lat","axis","Y")
ncatt_put(ncout,"time","axis","T")

nc_close(ncout)

} # End function

for(years in 2015)
{
jen_compare(years,winwid=3,winwidH=8,cold=F,closed=T,daily=T,
     jdir="/g/data/eg3/asp561/CattoData_20052015/",
#     jdir="/g/data/eg3/ajd548/vicci/front_data_JC/",
     dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_highs_rad10cv0.075/",
     outf=paste0("frontcomp_cattoall_highs8deg_proj100_rad10cv0.075_",years,"_daily.nc"))
}

