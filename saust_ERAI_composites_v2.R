library(ncdf4)
library(sp)
source("make_colours.R")

dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_highs_rad10cv0.075/"

snames=c("MAM","JJA","SON","DJF","MJJASO","NDJFMA","Annual")
mlist=list(3:5,6:8,9:11,c(1:2,12),5:10,c(1:4,11:12),1:12)
breaks=c(seq(0,3,0.25),100)

library(RColorBrewer)
cols=colorRampPalette(brewer.pal(9,"Blues"))(length(breaks)-1)

breaks2=c(-100,seq(-2,2,0.5),100)
cols2=colorRampPalette(rev(brewer.pal(11,"RdBu")))(length(breaks2)-1)

regnames=c("Aust","BigTas","GAB","SEIO")

## Get and save data

fixes=read.csv(paste0(dir,"UM_highs_ERAI_proj100_rad10cv0.075_bigaust_fixes.csv"),stringsAsFactors=F)
events=read.csv(paste0(dir,"UM_highs_ERAI_proj100_rad10cv0.075_bigaust_events.csv"),stringsAsFactors=F)
for(i in 1:length(events$ID))
{
  I=which(fixes$ID==events$ID[i])
  events$Move[i]=spDistsN1(cbind(fixes$Lon[min(I)],fixes$Lat[min(I)]),
                         cbind(fixes$Lon[max(I)],fixes$Lat[max(I)]),longlat=T)
}
  events=events[events$Move>=500,]
  fixes=fixes[fixes$Lat>=-50 & fixes$Lat<=-10 & fixes$Lon>=90 & fixes$Lon<=180,]
  fixes$Month=floor(fixes$Date/100)%%100

  include<-match(fixes$ID,events$ID)
  J<-which(is.na(include)==0 & (fixes$Open==0 | fixes$Open==10))

    fixes$BigTas<-fixes$GAB<-fixes$SEIO<-fixes$Aust<-0
    I<-which(fixes$Lon>=110 & fixes$Lon<=155 & fixes$Lat<(-10) & fixes$Lat>=-45)
    fixes$Aust[I]<-1
    I<-which(fixes$Lon>=120 & fixes$Lon<=145 & fixes$Lat<=(-27.5) & fixes$Lat>=-42.5)
    fixes$GAB[I]<-1
    I<-which(fixes$Lon>=90 & fixes$Lon<=110 & fixes$Lat<=(-27.5) & fixes$Lat>=-42.5)
    fixes$SEIO[I]<-1
    I<-which(fixes$Lon>=150 & fixes$Lon<=170 & fixes$Lat<=(-27.5) & fixes$Lat>=-42.5)
    fixes$BigTas[I]<-1

a=nc_open(paste0(dir,"UM_highs_ERAI_proj100_rad10cv0.075_bigaustrain.nc"))
lat=ncvar_get(a,"lat")
lon=ncvar_get(a,"lon")
cycrain=ncvar_get(a,"PRCP")
cycslp=ncvar_get(a,"SLP")
a=nc_open(paste0(dir,"UM_highs_ERAI_proj100_rad10cv0.075_bigaustT2.nc"))
cyctemp=ncvar_get(a,"T2")

for(r in 1:4)
for(s in 5:7)
{
S=which(fixes$Month[J]%in%mlist[[s]] & fixes[J,regnames[r]]==1)
avrain=apply(cycrain[,,J[S]],c(1,2),mean,na.rm=T)
avslp=apply(cycslp[,,J[S]],c(1,2),mean,na.rm=T)
avtemp=apply(cyctemp[,,J[S]],c(1,2),mean,na.rm=T)

## Make some plots

pdf(file=paste0("UM_highs_ERAI_proj100_rad10cv0.075_",regnames[r],"raintemp_500km_",snames[s],"_RdBu2.pdf"),height=4,width=10)
par(mar=c(2,2,3,1))
layout(cbind(1,3,2,4),width=c(1,0.25,1,0.25))
image(lon,lat,avrain[,length(lat):1],breaks=breaks,col=cols, # Have to reverse because ERAI is upside down
      xlab="",ylab="",main="6-hour rainfall (mm)")
contour(lon,lat,avslp[,length(lat):1],levels=seq(980,1040,4),lwd=2,col="black",add=T)
image(lon,lat,avtemp[,length(lat):1],breaks=breaks2,col=cols2, # Have to reverse because ERAI is upside down
      xlab="",ylab="",main="2m air temp anomaly (C)")
contour(lon,lat,avslp[,length(lat):1],levels=seq(980,1040,4),lwd=2,col="black",add=T)

ColorBar(breaks,cols,subsampleg=1,vert=T)
ColorBar(breaks2,cols2,subsampleg=1,vert=T)
dev.off()

}



