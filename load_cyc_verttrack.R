source("/short/eg3/asp561/R/init_cyclones.R")
library(geosphere)
library(ncdf4)

## Set up my projections and thresholds
## This should be the only thing I need to edit

sproj="proj240_rad2cv1"
sproj2="proj240_lows_rad2cv1"
uproj="proj240_rad2cv10"
uproj2="proj240_lows_rad2cv10"
lev=c(1000,925,850,700,600,500) # The 1000 is a dummy, loads surface lows

## The thresholds used. These are important, can set to vary with height or be constant
## Chosen by experimenting over the next section (when lows loaded), some examples given here

#thresh=c(1.5,11.5,11,10,11,13) # Closed + 2 fixes
thresh=c(1.5,13,11.5,10.5,12.5,17) # Closed, any length - may need altering for diff domains
#thresh=c(1.5,11,10,10,11,14.5) # Open + 2 fixes
#thresh=c(1.5,13,11.5,11,14,19.5) # Open, any length

lenlim=F # Do cyclones need to last for 2 fixes?
closed=T # Do cyclones need to be closed?
dist=500 # Distance for identifying upper cyclones within

years=1979:2016
latlim=c(-40,-25)
lonlim=c(148,160)
toplev=500 # Will usually be top level, but maybe want to track surface lows up to 300 sometimes, or track down from 700 or 850

## Done with setup, let's run this thing
## Load the cyclones for all levels

lows=list()

for(u in 1:(length(lev)))
{
  if(u==1)
  {
    tdir=paste("/short/eg3/asp561/cts.dir/gcyc_out/ERAI/",sproj2,"/",sep="")
    tmp=read.csv(paste(tdir,"UM_lows_ERAI_",sproj,"_bigaust_fixes.csv",sep=""))
    tmp$Date2=as.POSIXct(paste(as.character(tmp$Date),substr(tmp$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    tmp$Year=floor(tmp$Date/10000)
    tmp$Month=floor(tmp$Date/100)%%100
    tmp$Lon2=tmp$Lon
    tmp$Lon2[tmp$Lon>180]=tmp$Lon[tmp$Lon>180]-360
    
    if(lenlim)
    {
      tmp2=read.csv(paste(tdir,"UM_lows_ERAI_",sproj,"_bigaust_events.csv",sep=""))
      I=sort(tmp2$ID[which(tmp2$Length>1)])
      J=which(tmp$ID%in%I)
      tmp=tmp[J,]
    }
    
    if(closed) tmp=tmp[tmp$Open%%10==0,]
    lows[[u]]=tmp[tmp$CV>=thresh[u] & tmp$Year>=min(years) & tmp$Year<=max(years),]
    
  } else {
    tdir=paste("/short/eg3/asp561/cts.dir/gcyc_out/ERAI/",lev[u],"hPa_z/",uproj2,"/",sep="")
    tmp=read.csv(paste(tdir,"UM_lows_ERAI_",lev[u],"hPa_",uproj,"_bigaust_fixes.csv",sep=""))
    tmp$Date2=as.POSIXct(paste(as.character(tmp$Date),substr(tmp$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    tmp$Year=floor(tmp$Date/10000)
    tmp$Month=floor(tmp$Date/100)%%100
    tmp$Lon2=tmp$Lon
    tmp$Lon2[tmp$Lon>180]=tmp$Lon[tmp$Lon>180]-360
    
    if(lenlim)
    {
    tmp2=read.csv(paste(tdir,"UM_lows_ERAI_",lev[u],"hPa_",uproj,"_bigaust_events.csv",sep=""))
    I=sort(tmp2$ID[which(tmp2$Length>1)])
    J=which(tmp$ID%in%I)
    tmp=tmp[J,]
    }
    
    if(closed) tmp=tmp[tmp$Open%%10==0,]
    lows[[u]]=tmp[tmp$CV>=thresh[u] & tmp$Year>=min(years) & tmp$Year<=max(years),]
  }
}

for(u in 1:length(lev))
  print(length(unique(lows[[u]]$Date[lows[[u]]$Lon>=min(lonlim) & lows[[u]]$Lon<=max(lonlim) & 
                                           lows[[u]]$Lat>=min(latlim) & lows[[u]]$Lat<=max(latlim)]))/length(years))


## Now need to do two sets of vertical tracking, from surface and from top

I=which(lows[[1]]$Lon>=min(lonlim) & lows[[1]]$Lon<=max(lonlim) & 
          lows[[1]]$Lat>=min(latlim) & lows[[1]]$Lat<=max(latlim))
surfcyc=lows[[1]][I,]

surf_verttrack=array(NaN,c(length(surfcyc[,1]),length(lev),9))
dimnames(surf_verttrack)[[2]]=lev
dimnames(surf_verttrack)[[3]]=c("X","ID","Fix","Dist","Bearing","GPH","CV","Depth","Radius")

surf_verttrack[,1,1]=surfcyc[,1]
surf_verttrack[,1,2]=surfcyc$ID
surf_verttrack[,1,3]=surfcyc$Fix
surf_verttrack[,1,6]=surfcyc$MSLP*100/(9.80665*1.225) # Convert MSLP to GPH
surf_verttrack[,1,7]=surfcyc$CV
surf_verttrack[,1,8]=surfcyc$Depth
surf_verttrack[,1,9]=surfcyc$Radius

for(i in 1:length(surfcyc[,1]))
{
  surf=surfcyc[i,]
  
  for(u in 2:length(lev))
  {
    I=which(lows[[u]]$Date2==surf$Date2)
    
    if(length(I)>0)
    {
      # Find great circle distance between surf and next level, check if below threshold
      tmp=distGeo(cbind(lows[[u]]$Lon2[I],lows[[u]]$Lat[I]),cbind(surf$Lon2,surf$Lat))/1000
      J=which(tmp<=dist)
      
      # If multiple in range find closest
      if(length(J)>1)
      {
        L=which(tmp[J]==min(tmp[J]))
        J=J[L]
      }
      
      # If there is a match 
      if(length(J)==0) break else 
      {
        upper=lows[[u]][I[J],]
        surf_verttrack[i,u,1]=I[J]
        surf_verttrack[i,u,2]=upper$ID
        surf_verttrack[i,u,3]=upper$Fix
        surf_verttrack[i,u,4]=tmp[J]
        surf_verttrack[i,u,5]=bearing(cbind(surf$Lon2,surf$Lat),cbind(upper$Lon2,upper$Lat))%%360
        
        surf_verttrack[i,u,6]=upper$MSLP
        surf_verttrack[i,u,7]=upper$CV
        surf_verttrack[i,u,8]=upper$Depth
        surf_verttrack[i,u,9]=upper$Radius
        
        surf=upper # Next loop treat this one as the "base" to compare against
      }
    }
  }
}

apply(!is.na(surf_verttrack[,,1]),2,mean)

## Some extra variables for analysis

surfcyc$TopLevel=1000
for(i in 1:length(surfcyc[,1])) surfcyc$TopLevel[i]=lev[max(which(!is.na(surf_verttrack[i,,1])))]

surfcyc$TopCV<-surfcyc$TopBearing<-surfcyc$TopDist<-NaN
U=which(ulev==toplev)
I=which(!is.na(surf_verttrack[,U,1]))

for(i in I)
{
  surf=surfcyc[i,]
  upper=lows[[U]][surf_verttrack[i,U,1],]
  
  surfcyc$TopDist[i]=distGeo(cbind(surf$Lon2,surf$Lat),cbind(upper$Lon2,upper$Lat))/1000
  surfcyc$TopBearing[i]=bearing(cbind(surf$Lon2,surf$Lat),cbind(upper$Lon2,upper$Lat))%%360
  surfcyc$TopCV[i]=upper$CV
}

## And same but going down

U=which(lev==toplev)
I=which(lows[[U]]$Lon>=min(lonlim) & lows[[U]]$Lon<=max(lonlim) & 
          lows[[U]]$Lat>=min(latlim) & lows[[U]]$Lat<=max(latlim))
topcyc=lows[[U]][I,]

top_verttrack=array(NaN,c(length(topcyc[,1]),length(lev),9))
dimnames(top_verttrack)[[2]]=rev(lev)
dimnames(top_verttrack)[[3]]=c("X","ID","Fix","Dist","Bearing","GPH","CV","Depth","Radius")

top_verttrack[,1,1]=topcyc[,1]
top_verttrack[,1,2]=topcyc$ID
top_verttrack[,1,3]=topcyc$Fix
top_verttrack[,1,6]=topcyc$MSLP
top_verttrack[,1,7]=topcyc$CV
top_verttrack[,1,8]=topcyc$Depth
top_verttrack[,1,9]=topcyc$Radius

for(i in 1:length(topcyc[,1]))
{
  surf=topcyc[i,]
  
  n=1
  
  for(u in seq(U-1,1,-1))
  {
    n=n+1
    I=which(lows[[u]]$Date2==surf$Date2)
    
    if(length(I)>0)
    {
      # Find great circle distance between surf and next level, check if below threshold
      tmp=distGeo(cbind(lows[[u]]$Lon2[I],lows[[u]]$Lat[I]),cbind(surf$Lon2,surf$Lat))/1000
      J=which(tmp<=dist)
      
      # If multiple in range find closest
      if(length(J)>1)
      {
        L=which(tmp[J]==min(tmp[J]))
        J=J[L]
      }
      
      # If there is a match 
      if(length(J)==0) break else 
      {
        upper=lows[[u]][I[J],]
        top_verttrack[i,n,1]=I[J]
        top_verttrack[i,n,2]=upper$ID
        top_verttrack[i,n,3]=upper$Fix
        top_verttrack[i,n,4]=tmp[J]
        top_verttrack[i,n,5]=bearing(cbind(surf$Lon2,surf$Lat),cbind(upper$Lon2,upper$Lat))%%360
        
        if(u==1) top_verttrack[i,u,6]=upper$MSLP*100/(9.80665*1.225) else top_verttrack[i,u,6]=upper$MSLP
        top_verttrack[i,n,7]=upper$CV
        top_verttrack[i,n,8]=upper$Depth
        top_verttrack[i,n,9]=upper$Radius
        
        surf=upper # Next loop treat this one as the "base" to compare against
      }
    }
  }
}

apply(!is.na(top_verttrack[,,1]),2,mean)

## Some extra variables for analysis

topcyc$LowLevel=toplev
for(i in 1:length(topcyc[,1])) topcyc$LowLevel[i]=rev(lev)[max(which(!is.na(top_verttrack[i,,1])))]

topcyc$LowCV<-topcyc$LowBearing<-topcyc$LowDist<-NaN
I=which(!is.na(top_verttrack[,length(lev),1]))

for(i in I)
{
  surf=topcyc[i,]
  upper=lows[[1]][top_verttrack[i,length(lev),1],]
  
  topcyc$LowDist[i]=distGeo(cbind(surf$Lon2,surf$Lat),cbind(upper$Lon2,upper$Lat))/1000
  topcyc$LowBearing[i]=bearing(cbind(surf$Lon2,surf$Lat),cbind(upper$Lon2,upper$Lat))%%360
  topcyc$LowCV[i]=upper$CV
}

### Now, load the corresponding composite data for both directions
elat<-elon<-seq(-10.5,10.5,0.75)
lat<-lon<-seq(-10,10,0.25)

## Surface first

tdir=paste("/short/eg3/asp561/cts.dir/gcyc_out/ERAI/",sproj2,"/",sep="")
tmp=read.csv(paste(tdir,"UM_lows_ERAI_",sproj,"_bigaust_fixes.csv",sep=""))
tmp=tmp[tmp$Lon>=110 & tmp$Lon<=160 & tmp$Lat>=-45 & tmp$Lat<=-10,] ## This is what I did composites for
tmp$Year=floor(tmp$Date/10000)

## Now I need to get the same subset of cyclones as used for the tracking before
I=which(tmp[,1]%in%surfcyc[,1]) # Since i have that pesky first column, 

a=nc_open(paste(tdir,"composite_ERAI_",sproj,"_aust.nc",sep=""))
surf_slp=ncvar_get(a,"ECL_SLP")[,,I]
surf_ws=ncvar_get(a,"ECL_WS10")[,,I]

a=nc_open(paste(tdir,"raintemp_ERAI_",sproj,"_aust.nc",sep=""))
surf_gph500=ncvar_get(a,"ECL_GPH500")[,,I]
surf_prcp=ncvar_get(a,"ECL_PRCP")[,,I]
surf_tas=ncvar_get(a,"ECL_T2")[,,I]

a=nc_open(paste(tdir,"rain_TRMM_ERAI_",sproj,"_aust_19982015.nc",sep=""))
p=ncvar_get(a,"ECL_PRCP")
surf_prcpTRMM<-array(NaN,c(length(lat),length(lon),length(I)))
J=which(tmp$Year[I]>=1998 & tmp$Year[I]<=2015)
K=which(tmp$Year>=1998 & tmp$Year<=2015)
K2=which(tmp[K,1]%in%surfcyc[,1])
surf_prcpTRMM[,,J]=p[,,K2]

a=nc_open(paste(tdir,"lightning_ERAI_",sproj,"_aust_20052015.nc",sep=""))
p=ncvar_get(a,"ECL_lightning")
surf_lng<-array(NaN,c(length(lat),length(lon),length(I)))
J=which(tmp$Year[I]>=2005 & tmp$Year[I]<=2015)
K=which(tmp$Year>=2005 & tmp$Year<=2015)
K2=which(tmp[K,1]%in%surfcyc[,1])
surf_lng[,,J]=p[,,K2]

## Then upper

U=which(lev==toplev)
tdir=paste("/short/eg3/asp561/cts.dir/gcyc_out/ERAI/",lev[U],"hPa_z/",uproj2,"/",sep="")
tmp=read.csv(paste(tdir,"UM_lows_ERAI_",lev[U],"hPa_",uproj,"_bigaust_fixes.csv",sep=""))
tmp=tmp[tmp$Lon>=110 & tmp$Lon<=160 & tmp$Lat>=-45 & tmp$Lat<=-10,] ## This is what I did composites for
tmp$Year=floor(tmp$Date/10000)

## Now I need to get the same subset of cyclones as used for the tracking before
I=which(tmp[,1]%in%topcyc[,1]) # Since i have that pesky first column, 

a=nc_open(paste(tdir,"composite_ERAI_",lev[U],"hPa_",uproj,"_aust.nc",sep=""))
top_slp=ncvar_get(a,"ECL_SLP")[,,I]
top_ws=ncvar_get(a,"ECL_WS10")[,,I]

a=nc_open(paste(tdir,"raintemp_ERAI_",lev[U],"hPa_",uproj,"_aust.nc",sep=""))
top_gph500=ncvar_get(a,"ECL_GPH500")[,,I]
top_prcp=ncvar_get(a,"ECL_PRCP")[,,I]
top_tas=ncvar_get(a,"ECL_T2")[,,I]

a=nc_open(paste(tdir,"rain_TRMM_ERAI_",lev[U],"hPa_",uproj,"_aust_19982015.nc",sep=""))
p=ncvar_get(a,"ECL_PRCP")
top_prcpTRMM<-array(NaN,c(length(lat),length(lon),length(I)))
J=which(tmp$Year[I]>=1998 & tmp$Year[I]<=2015)
K=which(tmp$Year>=1998 & tmp$Year<=2015)
K2=which(tmp[K,1]%in%topcyc[,1])
top_prcpTRMM[,,J]=p[,,K2]

a=nc_open(paste(tdir,"lightning_ERAI_",lev[U],"hPa_",uproj,"_aust_20052015.nc",sep=""))
p=ncvar_get(a,"ECL_lightning")
top_lng<-array(NaN,c(length(lat),length(lon),length(I)))
J=which(tmp$Year[I]>=2005 & tmp$Year[I]<=2015)
K=which(tmp$Year>=2005 & tmp$Year<=2015)
K2=which(tmp[K,1]%in%topcyc[,1])
top_lng[,,J]=p[,,K2]

## Add some pertinent fields of impacts to both surf & top cyclones
## Make masks

dists=matrix(0,length(lat),length(lon))
for(i in 1:81)
  for(j in 1:81)
    dists[i,j]=(lat[i]^2 + lon[j]^2)^0.5
mask500=array(NaN,dim(dists))
for(i in 1:length(lon))
{
  I=which(dists[i,]<=5)
  mask500[i,I]=1
}

dists=matrix(0,length(elat),length(elon))
for(i in 1:29)
  for(j in 1:29)
    dists[i,j]=(elat[i]^2 + elon[j]^2)^0.5
emask500=array(NaN,dim(dists))
for(i in 1:length(elon))
{
  I=which(dists[i,]<=5)
  emask500[i,I]=1
}

## Surface

aa=dim(surfcyc)
surfcyc=cbind(surfcyc,array(NaN,c(length(surfcyc[,1]),10)))
colnames(surfcyc)=c(colnames(surfcyc[1:aa[2]]),"MEANWIND500","MAXWIND500",
                  "MEANT500","MAXT500","MINT500","TOTALLNG500",
                  "MEANRAIN500","MAXRAIN500","MEANRAIN500TRMM","MAXRAIN500TRMM")
for(i in 1:aa[1])
{
  surfcyc$MEANWIND500[i]=mean(surf_ws[,,i]*emask500,na.rm=T)
  surfcyc$MAXWIND500[i]=max(surf_ws[,,i]*emask500,na.rm=T)
  surfcyc$MEANT500[i]=mean(surf_tas[,,i]*emask500,na.rm=T)
  surfcyc$MAXT500[i]=max(surf_tas[,,i]*emask500,na.rm=T)
  surfcyc$MINT500[i]=min(surf_tas[,,i]*emask500,na.rm=T)
  surfcyc$MEANRAIN500[i]=mean(surf_prcp[,,i]*emask500,na.rm=T)
  surfcyc$MAXRAIN500[i]=max(surf_prcp[,,i]*emask500,na.rm=T)
  if(mean(!is.na(surf_prcpTRMM[21:61,21:61,i]))>=0.9)
  {
    surfcyc$MEANRAIN500TRMM[i]=mean(surf_prcpTRMM[,,i]*mask500,na.rm=T)
    surfcyc$MAXRAIN500TRMM[i]=max(surf_prcpTRMM[,,i]*mask500,na.rm=T)
  }
  if(mean(!is.na(surf_lng[21:61,21:61,i]))>=0.9)
  {
    surfcyc$TOTALLNG500[i]=sum(surf_lng[,,i]*mask500,na.rm=T)
  }
}

## And upper


aa=dim(topcyc)
topcyc=cbind(topcyc,array(NaN,c(length(topcyc[,1]),10)))
colnames(topcyc)=c(colnames(topcyc[1:aa[2]]),"MEANWIND500","MAXWIND500",
                    "MEANT500","MAXT500","MINT500","TOTALLNG500",
                    "MEANRAIN500","MAXRAIN500","MEANRAIN500TRMM","MAXRAIN500TRMM")
for(i in 1:aa[1])
{
  topcyc$MEANWIND500[i]=mean(top_ws[,,i]*emask500,na.rm=T)
  topcyc$MAXWIND500[i]=max(top_ws[,,i]*emask500,na.rm=T)
  topcyc$MEANT500[i]=mean(top_tas[,,i]*emask500,na.rm=T)
  topcyc$MAXT500[i]=max(top_tas[,,i]*emask500,na.rm=T)
  topcyc$MINT500[i]=min(top_tas[,,i]*emask500,na.rm=T)
  topcyc$MEANRAIN500[i]=mean(top_prcp[,,i]*emask500,na.rm=T)
  topcyc$MAXRAIN500[i]=max(top_prcp[,,i]*emask500,na.rm=T)
  if(mean(!is.na(top_prcpTRMM[21:61,21:61,i]))>=0.9)
  {
    topcyc$MEANRAIN500TRMM[i]=mean(top_prcpTRMM[,,i]*mask500,na.rm=T)
    topcyc$MAXRAIN500TRMM[i]=max(top_prcpTRMM[,,i]*mask500,na.rm=T)
  }
  if(mean(!is.na(top_lng[21:61,21:61,i]))>=0.9)
  {
    topcyc$TOTALLNG500[i]=sum(top_lng[,,i]*mask500,na.rm=T)
  }
}

### Everything is loaded and ready for analyses
### Huzzah, huzzah
