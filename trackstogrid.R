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

makegrid<-function(year,dir,outf,type="low",thresh=0,dur=NA,move=NA,closed=F,lat=seq(-90,90,1),lon=seq(0,360,1),winwid=10,widtype="index",wincirc=T,len=NA,lonlim=c(0,360),latlim=c(-90,90),expandf=F)
{

datelist=seq.POSIXt(as.POSIXct(paste0(year,"0101 00:00"),format="%Y%m%d %H:%M",tz="GMT"),as.POSIXct(paste0(year+1,"0101 00:00"),format="%Y%m%d %H:%M",tz="GMT"),by="6 hours")
datelist=datelist[-length(datelist)]

## Need an expanded lat & lon for doing initial analysis
griddata<-array(0,c(length(lon),length(lat),length(datelist)))
dlat=abs(lat[2]-lat[1])
dlon=abs(lon[2]-lon[1])
lat2=lat/dlat
lon2=lon/dlat

if(type=="front") {
fixes=read.table(paste0(dir,"ERAI_fronts_",year,".dat"),header=F,stringsAsFactors=F)
colnames(fixes)=c("ID","Point","Date","Time","Lat","Lon","Vdiff")
fixes$Year=floor(fixes$Date/10000)
fixes=fixes[fixes$Year==year,]
fixes$Lat2=round(fixes$Lat/dlat,0) 
fixes$Lon2=round((fixes$Lon%%360)/dlon,0)
fixes$Date2=as.POSIXct(paste0(year*10000+fixes$Date%%10000," ",fixes$Time/100,":00"),format="%Y%m%d %H:%M",tz="GMT")

x<-rle(fixes$ID)
events<-data.frame(ID=x$values,Length=x$lengths)

if(!is.na(len)) ## Get rid of short fronts
{
  events=events[events[,2]>=len,]
  include<-match(fixes[,1],events[,1])
  J<-which(is.na(include)==0)
  fixes=fixes[J,]
}

if(expandf)
{
init=F
for(i in 1:length(events$ID))
{
I=which(fixes$ID==events$ID[i])
tmplat=lat2[lat2>=min(fixes$Lat2[I]) & lat2<=max(fixes$Lat2[I])]
if(length(tmplat)>0)
{
 tmplon=approx(fixes$Lat2[I],fixes$Lon[I],tmplat)
 tmplon2=round((tmplon$y%%360)/dlon,0)
 tmp=data.frame(Lon2=tmplon2,Lat2=tmplat,Date2=rep(fixes$Date2[I[1]],length(tmplat)))

 if(init) fixes2=rbind(fixes2,tmp) else {
   fixes2=tmp
   init=T
 }
 }
}
fixes=fixes2
}

grid1=table(factor(fixes$Lon2,levels=lon2),factor(fixes$Lat2,levels=lat2),factor(as.character(fixes$Date2),levels=as.character(datelist))) ## This is the lazy approach, no interpolation of the front to missing latitudes

} else {

fname=paste(dir,"/tracks_",year,".dat",sep="")
read.table(fname, sep="",skip=1)->fixes
colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Meh")
fixes$Year=floor(fixes$Date/10000)
if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
fixes$Month=floor(fixes$Date/100)%%100
fixes$CV=abs(fixes$CV)
fixes$Lat2=round(fixes$Lat/dlat,0) 
fixes$Lon2=round((fixes$Lon%%360)/dlon,0)
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
  fixes=fixes[fixes$CV>=thresh,]
  if(closed) fixes=fixes[(fixes$Open==0 | fixes$Open==10),] # Remove all open systems
  grid1=table(factor(fixes$Lon2,levels=lon2),factor(fixes$Lat2,levels=lat2),factor(as.character(fixes$Date2),levels=as.character(datelist))) ## This is the lazy approach, no interpolation of the front to missing latitudes

}

## Spread to a region based on winwid

if((widtype=="deg" | widtype=="degree") & wincirc==F) winwid=round(winwid/dlon,0)

for(i in 1:length(datelist)) griddata[,,i]=spreadeffect(grid1[,,i],lon=lon,lat=lat,circ=wincirc,winwid=winwid)

## Save as a netcdf

I=which(lon>=min(lonlim) & lon<=max(lonlim))
J=which(lat>=min(latlim) & lat<=max(latlim))

dimX<-ncdim_def("lon","degrees_E",lon[I])
dimY<-ncdim_def("lat","degrees_N",lat[J])
dimT<-ncdim_def("time","hours since 1970-1-1 00:00:00",as.numeric(datelist)/(60*60))

fillvalue <- 1e32
cyc_def <- ncvar_def("systems","count",list(dimX,dimY,dimT),fillvalue,paste("Cells affected by a ",type),prec="single")

# create netCDF file and put arrays
ncout <- nc_create(paste0(dir,outf),cyc_def) #force_v4=T)

# put variables
ncvar_put(ncout,cyc_def,griddata[I,J,])

# put additional attributes into dimension and data variables
ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout,"lat","axis","Y")
ncatt_put(ncout,"time","axis","T")

nc_close(ncout)
}

#for(year in 1980:2016) makegrid(year,dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_highs_rad10cv0.075/",outf=paste0("austhighs_7deg_proj100_rad10cv0.075_D2_",year,"_vic.nc"),type="high",dur=2,closed=T,lat=seq(-50,-20,0.25),lon=seq(120,160,0.25),latlim=c(-40,-32.5),lonlim=c(140,151),winwid=7,widtype="deg",wincirc=T)

#for(year in 1979:2015) makegrid(year,dir="/short/eg3/asp561/Fronts/RudevaFiles/",outf=paste0("austfronts_len3_3deg_",year,"_vic.nc"),type="front",lat=seq(-50,-20,0.25),lon=seq(120,160,0.25),latlim=c(-40,-32.5),lonlim=c(140,151),winwid=3,widtype="deg",wincirc=T,len=3,expandf=T)

#for(year in 1979:2015) makegrid(year,dir="/short/eg3/asp561/Fronts/RudevaFiles/",outf=paste0("austfronts_len3_3deg_",year,".nc"),type="front",lat=seq(-70,0,0.5),lon=seq(90,180,0.5),latlim=c(-50,-4),lonlim=c(110,161),winwid=3,widtype="deg",wincirc=T,len=3,expandf=T)


#for(year in 1979:2016) makegrid(year,dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_lows_rad5cv0.15/",outf=paste0("austlows_6deg_proj100_rad5cv0.15_open_D2_",year,"_vic.nc"),type="low",dur=2,closed=F,lat=seq(-50,-20,0.25),lon=seq(120,160,0.25),latlim=c(-40,-32.5),lonlim=c(140,151),winwid=6,widtype="deg",wincirc=T)

#for(year in 1979:2016) makegrid(year,dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/850hPa_z/proj100_lows_rad5cv1/",outf=paste0("austlows_6deg_850hPa_proj100_rad5cv1_D2_",year,"_vic.nc"),type="low",dur=2,closed=T,lat=seq(-50,-20,0.25),lon=seq(120,160,0.25),latlim=c(-40,-32.5),lonlim=c(140,151),winwid=6,widtype="deg",wincirc=T)

for(year in 1979:2016) makegrid(year,dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_lows_rad5cv0.15_notopo/",outf=paste0("austlows_6deg_proj100_rad5cv0.15_notopo_open_",year,"_vic.nc"),type="low",closed=F,lat=seq(-50,-20,0.25),lon=seq(120,160,0.25),latlim=c(-40,-32.5),lonlim=c(140,151),winwid=6,widtype="deg",wincirc=T)

