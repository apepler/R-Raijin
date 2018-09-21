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
#thresh=c(1.5,13,11.5,10.5,12.5,17) # Closed, any length - may need altering for diff domains
#thresh=c(1.4,12,11,10,12,16) # Closed, any length - alternative with higher freq
#thresh=c(1.5,11,10,10,11,14.5) # Open + 2 fixes
#thresh=c(1.5,13,11.5,11,14,19.5) # Open, any length
#thresh=c(1.4,12,10.5,10,13,18.5) # Open, any length, v2
thresh=c(1.4,11.5,11.5,11.5,11.5,11.5) # New constant freq


lenlim=F # Do cyclones need to last for 2 fixes?
closed=F # Do cyclones need to be closed?
dist=500 # Distance for identifying upper cyclones within

years=1979:2016
toplev=500 # Will usually be top level, but maybe want to track surface lows up to 300 sometimes, or track down from 700 or 850


lat=seq(-89.5,89.5,1)
lon=seq(0.5,359.5,1)

freqvheight_surf<-freqvheight_upper<-array(0,c(length(lon),length(lat),length(lev),length(years),12))

## Done with setup, let's run this thing
## Load the cyclones for all levels


for(y in 1:length(years))
{
  print(years[y])
  lows=list()
  for(u in 1:(length(lev)))
  {
    if(u==1) dir=paste("/short/eg3/asp561/cts.dir/gcyc_out/ERAI/",sproj2,"/",sep="") else
      dir=paste("/short/eg3/asp561/cts.dir/gcyc_out/ERAI/",lev[u],"hPa_z/",uproj2,"/",sep="")
    
    fname=paste(dir,"/tracks_",years[y],".dat",sep="")
    read.table(fname, sep="",skip=1)->fixes
    colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Depth","Radius","Up","Vp")
    fixes$Year=floor(fixes$Date/10000)
    fixes$Month=floor(fixes$Date/100)%%100
    if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
    fixes$Date2=as.POSIXct(paste(as.character(years[y]*10000+fixes$Date%%10000),substr(fixes$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    if(closed) fixes=fixes[fixes$Open%%10==0,]
    fixes$Lon2=fixes$Lon
    fixes$Lon2[fixes$Lon>180]=fixes$Lon[fixes$Lon>180]-360
    lows[[u]]=fixes[abs(fixes$CV)>=thresh[u],] 
  }

## Now need to do two sets of vertical tracking, from surface and from top
surfcyc=lows[[1]]

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

## Some extra variables for analysis

surfcyc$TopLevel=1000
for(i in 1:length(surfcyc[,1])) surfcyc$TopLevel[i]=lev[max(which(!is.na(surf_verttrack[i,,1])))]

## And same but going down

U=which(lev==toplev)
topcyc=lows[[U]]

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

topcyc$LowLevel=toplev
for(i in 1:length(topcyc[,1])) topcyc$LowLevel[i]=rev(lev)[max(which(!is.na(top_verttrack[i,,1])))]

### Now plot their vertical distribution

surfcyc$Lat2=floor(surfcyc$Lat)
surfcyc$Lon2=floor(surfcyc$Lon)%%360
freqvheight_surf[,,,y,]=table(factor(surfcyc$Lon2,levels=0:359),factor(surfcyc$Lat2,levels=-90:89),factor(surfcyc$TopLevel,levels=lev),surfcyc$Month)

topcyc$Lat2=floor(topcyc$Lat)
topcyc$Lon2=floor(topcyc$Lon)%%360
freqvheight_upper[,,,y,]=table(factor(topcyc$Lon2,levels=0:359),factor(topcyc$Lat2,levels=-90:89),factor(topcyc$LowLevel,levels=lev),topcyc$Month)
}

save(lat,lon,years,freqvheight_surf,lev,sproj,uproj,thresh,dist,freqvheight_upper,file="globalcyclones_ERAIproj240rad2_depthbylocation_constantthresh1.4_open.RData")


## 1. By season (top & bottom), plot of # of cyclones in these three categories

mlist=cbind(5:10,c(11:12,1:4))
snames=c("MJJASO","NDJFMA")
names=c("Surface-only cyclones","Mid-level cyclones",
        "Deep cyclones","Upper-only cyclones")

breaks1=c(0,0.05,seq(0.1,1,0.1),100)
col1=colorRampPalette(brewer.pal(9,"Blues"))(length(breaks1))
col1=col1[-1]

pdf(file="globalcyclones_ERAIproj240rad2_depthbylocation_constantthresh1.4_open_smooth2.pdf",width=10,height=12,pointsize=12)
layout(cbind(c(1:4),c(5:8),rep(9,4)),width=c(1,1,0.3))
par(mar=c(2,2,3,1),cex.lab=1.2,cex.axis=1.2)

for(m in 1:2)
{
  freq=abind(apply(freqvheight_surf[,,,,mlist[,m]],c(1,2,3),sum),apply(freqvheight_upper[,,,,mlist[,m]],c(1,2,3),sum),along=4)
  
  image(lon,lat,makesmooth(apply(freq[,,lev>=850,1],c(1,2),sum),winwid=2)/length(years),
        breaks=breaks1,col=col1,main=paste(names[1],"-",snames[m]),xlab="",ylab="",cex.main=1.5)
  contour(lon,lat,makesmooth(apply(freq[,,lev>=850,1],c(1,2),sum),winwid=5)/length(years),
          levels=breaks1[seq(2,length(breaks1),2)],add=T,lwd=1.5,col="black",drawlabels=F)
  map("world2",add=T)
  box()
  
  image(lon,lat,makesmooth(apply(freq[,,lev%in%c(600,700),1],c(1,2),sum),winwid=2)/length(years),
        breaks=breaks1,col=col1,main=paste(names[2],"-",snames[m]),xlab="",ylab="",cex.main=1.5)
  contour(lon,lat,makesmooth(apply(freq[,,lev%in%c(600,700),1],c(1,2),sum),winwid=5)/length(years),
          levels=breaks1[seq(2,length(breaks1),2)],add=T,lwd=1.5,col="black",drawlabels=F)
  map("world2",add=T)
  box()
  
  image(lon,lat,makesmooth(freq[,,lev==500,1],winwid=2)/length(years),
        breaks=breaks1,col=col1,main=paste(names[3],"-",snames[m]),xlab="",ylab="",cex.main=1.5)
  contour(lon,lat,makesmooth(freq[,,lev==500,1],winwid=5)/length(years),
          levels=breaks1[seq(2,length(breaks1),2)],add=T,lwd=1.5,col="black",drawlabels=F)
  map("world2",add=T)
  box()
  
  image(lon,lat,makesmooth(apply(freq[,,lev<=700,2],c(1,2),sum),winwid=2)/length(years),
        breaks=breaks1,col=col1,main=paste(names[4],"-",snames[m]),xlab="",ylab="",cex.main=1.5)
  contour(lon,lat,makesmooth(apply(freq[,,lev<=700,2],c(1,2),sum),winwid=5)/length(years),
          levels=breaks1[seq(2,length(breaks1),2)],add=T,lwd=1.5,col="black",drawlabels=F)
  map("world2",add=T)
  box()
}
ColorBar(brks=breaks1,cols=col1,vert=T)
dev.off()

### Shallow proportion, by season

mlist=cbind(5:10,c(11:12,1:4))
snames=c("MJJASO","NDJFMA")
names=c("Shallow % (surface)","Deep % (surface)","Shallow % (upper)")

breaks2=c(seq(0,50,5),100)
col2=colorRampPalette(brewer.pal(9,"Blues"))(length(breaks2)-1)

pdf(file="globalcyclones_ERAIproj240rad2_shallowpropbylocation_constantthresh1.4_open_smooth2.pdf",width=10,height=9,pointsize=12)
layout(cbind(c(1:3),c(4:6),rep(7,3)),width=c(1,1,0.3))
par(mar=c(2,2,3,1),cex.lab=1.2,cex.axis=1.2)

for(m in 1:2)
{
  freq=abind(apply(freqvheight_surf[,,,,mlist[,m]],c(1,2,3),sum),apply(freqvheight_upper[,,,,mlist[,m]],c(1,2,3),sum),along=4)
  
  tmp=100*makesmooth(apply(freq[,,lev>=850,1],c(1,2),sum),winwid=2)/makesmooth(apply(freq[,,,1],c(1,2),sum),winwid=2)
  tmp[makesmooth(apply(freq[,,,1],c(1,2),sum),winwid=2)<1]=NaN
  image(lon,lat,tmp,breaks=breaks2,col=col2,
        main=paste(names[1],"-",snames[m]),xlab="",ylab="",cex.main=1.5)
  tmp=100*makesmooth(apply(freq[,,lev>=850,1],c(1,2),sum),winwid=5)/makesmooth(apply(freq[,,,1],c(1,2),sum),winwid=5)
  tmp[makesmooth(apply(freq[,,,1],c(1,2),sum),winwid=5)<1]=NaN
  contour(lon,lat,tmp,
          levels=breaks2[seq(3,length(breaks2),2)],add=T,lwd=1.5,col="black",drawlabels=F)
  map("world2",add=T)
  box()
  
  tmp=100*makesmooth(freq[,,lev==500,1],winwid=2)/makesmooth(apply(freq[,,,1],c(1,2),sum),winwid=2)
  tmp[makesmooth(apply(freq[,,,1],c(1,2),sum),winwid=2)<1]=NaN
  image(lon,lat,tmp,breaks=breaks2,col=col2,
        main=paste(names[2],"-",snames[m]),xlab="",ylab="",cex.main=1.5)
  tmp=100*makesmooth(freq[,,lev==500,1],winwid=5)/makesmooth(apply(freq[,,,1],c(1,2),sum),winwid=5)
  tmp[makesmooth(apply(freq[,,,1],c(1,2),sum),winwid=5)<1]=NaN
  contour(lon,lat,tmp,
          levels=breaks2[seq(3,length(breaks2),2)],add=T,lwd=1.5,col="black",drawlabels=F)
  map("world2",add=T)
  box()
  
  tmp=100*makesmooth(apply(freq[,,lev<=700,2],c(1,2),sum),winwid=2)/makesmooth(apply(freq[,,,2],c(1,2),sum),winwid=2)
  tmp[makesmooth(apply(freq[,,,2],c(1,2),sum),winwid=2)<1]=NaN
  image(lon,lat,tmp,breaks=breaks2,col=col2,
        main=paste(names[3],"-",snames[m]),xlab="",ylab="",cex.main=1.5)
  tmp=100*makesmooth(apply(freq[,,lev<=700,2],c(1,2),sum),winwid=5)/makesmooth(apply(freq[,,,2],c(1,2),sum),winwid=5)
  tmp[makesmooth(apply(freq[,,,2],c(1,2),sum),winwid=5)<1]=NaN
  contour(lon,lat,tmp,
          levels=breaks2[seq(3,length(breaks2),2)],add=T,lwd=1.5,col="black",drawlabels=F)
  map("world2",add=T)
  box()
}
ColorBar(brks=breaks2,cols=col2)
dev.off()

