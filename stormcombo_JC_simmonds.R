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


storm_combo<-function(year,outf=NA,outdir="",fdir,fdirS,cdir,cdirS,tdir,lonlim=c(-180,180),latlim=c(-90,90),winwid=0,cold=T,daily=F,TSenv=T,closed=T,cvthresh=0,len=0,shiftR=T,expandF=F)
{

datelist=seq.POSIXt(as.POSIXct(paste0(year,"0101 00:00"),format="%Y%m%d %H:%M",tz="GMT"),as.POSIXct(paste0(year+1,"0101 00:00"),format="%Y%m%d %H:%M",tz="GMT"),by="6 hours")
datelist=datelist[-length(datelist)]

## First, load fronts

if(cold)
{
 a=nc_open(paste0(fdir,"cold_fronts_Aus",year,".nc"))
 front=ncvar_get(a,"Cold fronts for use in collaboration with A. Dowdy & J. Catto:")
 flat=ncvar_get(a,"lat")
 flon=ncvar_get(a,"lon")

 if(winwid>0) 
 {
 tmp=front
 for(i in 1:length(datelist)) front[,,i]=spreadeffect(tmp[,,i],lon=flon,lat=flat,circ=T,winwid=winwid)
 }
} else {
 a=nc_open(paste0(fdir,"storm_types_Aus333_",year,".nc"))
 tmp=ncvar_get(a,"Storm types for use in collaboration with A. Dowdy & J. Catto:")
 front<-array(0,dim(tmp))
 front[tmp%in%c(2,4,6,7)]=1 ## Fronts
 flat=ncvar_get(a,"lat")
 flon=ncvar_get(a,"lon")
}

print("C fronts done")

## And make the corresponding Rudeva fronts
dlat=abs(flat[2]-flat[1])
dlon=abs(flon[2]-flon[1])
slat=seq(-90,90,dlat)
if(min(flon)<0) slon=seq(-180,180,dlon) else slon=seq(0,360,dlon)
lat2=slat/dlat
lon2=slon/dlon

fixes=read.table(paste0(fdirS,"ERAI_fronts_",year,".dat"),header=F,stringsAsFactors=F)
colnames(fixes)=c("ID","Point","Date","Time","Lat","Lon","Vdiff")
fixes$Year=floor(fixes$Date/10000)
fixes$Lat2=round(fixes$Lat/dlat,0)
fixes$Lon2=fixes$Lon%%360
if(min(flon)<0) fixes$Lon2[fixes$Lon2>180]=fixes$Lon2[fixes$Lon2>180]-360
fixes$Lon2=round(fixes$Lon2/dlon,0)
fixes$Date2=as.POSIXct(paste0(fixes$Date," ",fixes$Time/100,":00"),format="%Y%m%d %H:%M",tz="GMT")
if(shiftR) fixes$Date2=fixes$Date2+60*60*6 # Attribute front to next 6 hour

x<-rle(fixes$ID)
events<-data.frame(ID=x$values,Length=x$lengths)

if(!is.na(len)) ## Get rid of short fronts
{
  events=events[events[,2]>=len,]
  include<-match(fixes[,1],events[,1])
  J<-which(is.na(include)==0)
  fixes=fixes[J,]
}

if(expandF) ## Interpolates between points, if we care
{
init=F
for(i in 1:length(events$ID))
{
 I=which(fixes$ID==events$ID[i])
 tmplat=lat2[lat2>=min(fixes$Lat2[I]) & lat2<=max(fixes$Lat2[I])]
 if(length(tmplat)>1)
 {
  tmplon=approx(fixes$Lat2[I],fixes$Lon[I],tmplat)
  tmplon2=round((tmplon$y%%360)/dlon,0)
  tmp=data.frame(Lon2=tmplon2,Lat2=tmplat,Date2=rep(fixes$Date2[I[1]],length(tmplat)))
 }  else tmp=data.frame(Lon2=fixes$Lon2[I],Lat2=fixes$Lat2[I],Date2=fixes$Date2[I])

 if(init) fixes2=rbind(fixes2,tmp) else {
   fixes2=tmp
   init=T
 }
}
fixes=fixes2
}

print("R fronts done")

## Only keep areas affected by both Cattop & Rudeva fronts

grid1=table(factor(fixes$Lon2,levels=lon2),factor(fixes$Lat2,levels=lat2),factor(as.character(fixes$Date2),levels=as.character(datelist)))
ir_front=array(0,dim(grid1))
if(winwid>0) for(i in 1:length(datelist)) ir_front[,,i]=spreadeffect(grid1[,,i],lon=slon,lat=slat,circ=T,winwid=winwid)

I=which(slon%in%flon)
J=which(slat%in%flat)
ir_front=ir_front[I,J,]

front[ir_front==0]=0 # If not a Rudeva front, remove

print("Combo fronts done")

## Next, load cyclones

a=nc_open(paste0(cdir,"cyclones_",year,"_0.75.nc"))
jc_cyc=ncvar_get(a,"IDCLUST")
clat=ncvar_get(a,"latitude")
clon=ncvar_get(a,"longitude")

cyclone=array(0,dim(jc_cyc))
dlat=abs(clat[2]-clat[1])
dlon=abs(clon[2]-clon[1])
lat2=clat/dlat
lon2=clon/dlon

print("C lows done")

## And Simmonds cyclones

fname=paste(cdirS,"/tracks_",year,".dat",sep="")
read.table(fname, sep="",skip=1)->fixes
colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Meh")
fixes$Year=floor(fixes$Date/10000)
if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
fixes$Month=floor(fixes$Date/100)%%100
fixes$Lat2=round(fixes$Lat/dlat,0)
fixes$Lon2=fixes$Lon%%360
if(min(clon)<0) fixes$Lon2[fixes$Lon2>180]=fixes$Lon2[fixes$Lon2>180]-360
fixes$Lon2=round(fixes$Lon2/dlon,0)
fixes$Date2=as.POSIXct(paste0(year*10000+fixes$Date%%10000," ",fixes$Time),format="%Y%m%d %H:%M",tz="GMT")
if(closed) fixes=fixes[(fixes$Open==0 | fixes$Open==10),] # Remove all open systems
fixes=fixes[fixes$CV>=cvthresh,]

  grid1=table(factor(fixes$Lon2,levels=lon2),factor(fixes$Lat2,levels=lat2),factor(as.character(fixes$Date2),levels=as.character(datelist))) ## Nearest location to cyclone centre

print("S lows done")

# Only keep Catto cyclones where there was a simmonds cyclone co-located

for(t in 1:length(datelist))
{
tmp1=jc_cyc[,,t]
tmp2=grid1[,,t] 
tmp3=array(0,dim(tmp1))
idlist=unique(tmp1[tmp1>0]) # List of IDs
for(i in idlist)
{
I=which(tmp1==i)
if(max(tmp2[I])>0) tmp3[I]=1
}
cyclone[,,t]=tmp3
}

tmp=cyclone
if(winwid>0)
{
 for(i in 1:length(datelist)) cyclone[,,i]=spreadeffect(tmp[,,i],lon=clon,lat=clat,circ=T,winwid=winwid)
} else cyclone[tmp>0]=1

print("Combo lows done")

## Finally, load thunderstorms

if(TSenv)
{
 a=nc_open(paste0(tdir,"TS_environments_",year,".nc"))
 ts=ncvar_get(a,"Thunderstorm_environments_produced_by_ADowdy_May_2018:")
 tlat=ncvar_get(a,"lat")
 tlon=ncvar_get(a,"lon")

 if(winwid>0)
 {
 tmp=ts
 for(i in 1:length(datelist)) ts[,,i]=spreadeffect(tmp[,,i],lon=tlon,lat=tlat,circ=T,winwid=winwid)
 }
} else {
 a=nc_open(paste0(tdir,"storm_types_Aus333_",year,".nc"))
 tmp=ncvar_get(a,"Storm types for use in collaboration with A. Dowdy & J. Catto:")
 ts<-array(0,dim(tmp))
 ts[tmp%in%c(3,5,6,7)]=1 
 tlat=ncvar_get(a,"lat")
 tlon=ncvar_get(a,"lon")
} 

print("TS done")

## Convert to daily data if necessary

if(daily)
{
datelist2=as.numeric(format(datelist,"%Y%m%d"))
daylist=unique(datelist2)

tmp<-array(0,c(length(flon),length(flat),length(daylist)))
for(t in 1:length(daylist)) tmp[,,t]=apply(front[,,datelist2==daylist[t]],c(1,2),max)
front=tmp

tmp<-array(0,c(length(clon),length(clat),length(daylist)))
for(t in 1:length(daylist)) tmp[,,t]=apply(cyclone[,,datelist2==daylist[t]],c(1,2),max)
cyclone=tmp

tmp<-array(0,c(length(tlon),length(tlat),length(daylist)))
for(t in 1:length(daylist)) tmp[,,t]=apply(ts[,,datelist2==daylist[t]],c(1,2),max)
ts=tmp
}

## Combine together - need to use the smallest common region
## At least we know they're all 0.75 degree resolution

lonlim2=c(max(c(min(clon),min(flon),min(tlon))),
          min(c(max(clon),max(flon),max(tlon))))
latlim2=c(max(c(min(clat),min(flat),min(tlat))),
          min(c(max(clat),max(flat),max(tlat))))

I=which(clon>=min(lonlim2) & clon<=max(lonlim2))
J=which(clat>=min(latlim2) & clat<=max(latlim2))
lon=clon[I]
lat=clat[J]
cyclone=cyclone[I,J,]

I=which(flon>=min(lonlim2) & flon<=max(lonlim2))
J=which(flat>=min(latlim2) & flat<=max(latlim2))
front=front[I,J,]

I=which(tlon>=min(lonlim2) & tlon<=max(lonlim2))
J=which(tlat>=min(latlim2) & tlat<=max(latlim2))
ts=ts[I,J,]

## Now use these to give the regional types

storm_combo<-array(0,dim(ts))
storm_combo[cyclone==1 & front==0 & ts==0]=1 # Cyclone only
storm_combo[cyclone==0 & front==1 & ts==0]=2 # Front only
storm_combo[cyclone==0 & front==0 & ts==1]=3 # TS only
storm_combo[cyclone==1 & front==1 & ts==0]=4 # CF
storm_combo[cyclone==1 & front==0 & ts==1]=5 # CT
storm_combo[cyclone==0 & front==1 & ts==1]=6 # FT
storm_combo[cyclone==1 & front==1 & ts==1]=7 # CFT


I=which(lon>=min(lonlim) & lon<=max(lonlim))
J=which(lat>=min(latlim) & lat<=max(latlim))

dimX<-ncdim_def("lon","degrees_E",lon[I])
dimY<-ncdim_def("lat","degrees_N",lat[J])

if(daily) {
tmp=seq.POSIXt(as.POSIXct(paste0(year,"0101 09:00"),format="%Y%m%d %H:%M",tz="GMT"),as.POSIXct(paste0(year,"1231 09:00"),format="%Y%m%d %H:%M",tz="GMT"),by="24 hours")
dimT<-ncdim_def("time","hours since 1970-1-1 00:00:00",as.numeric(tmp)/(60*60))  
} else  dimT<-ncdim_def("time","hours since 1970-1-1 00:00:00",as.numeric(datelist)/(60*60))

fillvalue <- 9999
print(unique(as.numeric(storm_combo),na.rm=T))
cyc_def <- ncvar_def("storm_type","count",list(dimX,dimY,dimT),fillvalue,"Storm combinations from Catto and Dowdy 2017 - Cyclone Only (1), Front Only (2), Thunderstorm Only (3), Cyclone-Front (4), Cyclone-Thunderstorm (5), Front-Thunderstorm (6), Triple Storm (7)",prec="integer")

# create netCDF file and put arrays
ncout <- nc_create(paste0(outdir,outf),cyc_def) #force_v4=T)

# put variables
ncvar_put(ncout,cyc_def,storm_combo[I,J,])

# put additional attributes into dimension and data variables
ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout,"lat","axis","Y")
ncatt_put(ncout,"time","axis","T")

nc_close(ncout)

} # End function

for(years in 1979:2014)
{
storm_combo(years,winwid=3,daily=T,
     fdir="/g/data/eg3/ajd548/vicci/front_data_JC/",cold=T,
     cdir="/g/data/eg3/ajd548/vicci/cyclone_data_JC/",
#     tdir="/g/data/eg3/asp561/CattoData_20052015/",TSenv=F,
     tdir="/g/data/eg3/ajd548/vicci/TS_environments/",TSenv=T,
     fdirS="/g/data/eg3/asp561/Fronts/",len=2,shiftR=T,expandF=T,
     cdirS="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_lows_rad5cv0.15/",closed=T,
     outdir="/g/data/eg3/asp561/CattoData_20052015/acaciacombos/",
     outf=paste0("stormcombo_cattodowdysimmonds_TSenv_vic_",years,"_daily.nc"))#,
#     lonlim=c(110,160),latlim=c(-45,-10))
}

