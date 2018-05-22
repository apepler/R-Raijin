source("init_cyclones.R")
library(fields)
reanal="BARRA"
dir=paste0("/short/eg3/asp561/cts.dir/gcyc_out/",reanal,"/proj240_lows_rad2cv1/")
ylim=c(2011,2015)
#lat<-lon<-seq(-10.5,10.5,length.out=29)
lat<-lon<-seq(-10.01,10.01,length.out=183)

ColorBar2 <- function(brks,cols,vert=T,subsampleg=1)
{
  if(vert) {
    par(mar = c(2, 1, 2, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
    image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols,
          xlab = '', ylab = '')
    box()
    axis(4, at = seq(0.5, length(brks) - 0.5, subsampleg), tick = TRUE,
         labels = brks[seq(1, length(brks)-1, subsampleg)])
  } else {
    par(mar = c(1.5, 1, 1, 1), mgp = c(1.5, 0.3, 0), las = 1, cex = 1)
    image(1:length(cols), 1, t(t(1:length(cols))), axes = FALSE, col = cols,
          xlab = '', ylab = '')
    box()
    axis(1, at = seq(0.5, length(brks) - 0.5, subsampleg),
         labels = brks[seq(1, length(brks)-1, subsampleg)])
  }
}


pdf(file=paste("cyclone_",reanal,"_wind_rain_panel_ECL.pdf",sep=""),width=9,height=8,pointsize=12)
layout(cbind(c(1,3),c(2,4)))
par(mar=c(3,3,3,6),cex.lab=1.2,cex.axis=1.2)
pal1 <- color.palette(c("darkred","red","yellow","white","cyan","blue","darkblue"),c(10,20,10,10,20,10))

fixes=read.csv(paste0(dir,"UM_lows_",reanal,"_proj240_rad2cv1_bigaust_fixes.csv"))
fixes=fixes[fixes$Date>=ylim[1]*10000 & fixes$Date<=(ylim[2]+1)*10000 & fixes$Lon>=110 & fixes$Lon<=160 & fixes$Lat>=-45 & fixes$Lat<=-10,]
print(dim(fixes))

fixes$Location=0
#fixes$Location[fixes$Lon>=110 & fixes$Lon<=155 & fixes$Lat>=-45 & fixes$Lat<=-10]=1
#fixes$Location[fixes$Lon>=135 & fixes$Lon<=148 & fixes$Lat>=-42.5 & fixes$Lat<=-32.5]=1

    I<-which(fixes$Lon>=149 & fixes$Lon<=161 & fixes$Lat<(-37) & fixes$Lat>=-41)
    fixes$Location[I]<-1
    I<-which(fixes$Lon>=(149+(37+fixes$Lat)/2) & fixes$Lon<=161 & fixes$Lat<(-31) & fixes$Lat>=-37)
    fixes$Location[I]<-1
    I<-which(fixes$Lon>=152 & fixes$Lon<=161 & fixes$Lat<=(-24) & fixes$Lat>=-31)
    fixes$Location[I]<-1

I=which(fixes$Location==1)

if(reanal=="BARRA") {
fixes$Year=floor(fixes$Date/10000)
for(y in seq(ylim[1],ylim[2]))
{
 I=which(fixes$Location[fixes$Year==y]==1)
 a=nc_open(paste0(dir,"ECLcomposite_",reanal,"_proj240_rad2cv1_aust_",y,".nc"))
 if(y==ylim[1])
 {
  tmp=ncvar_get(a,"ECL_U10")
  u=apply(tmp[,,I],c(1,2),sum,na.rm=T)
  tmp=ncvar_get(a,"ECL_V10")
  v=apply(tmp[,,I],c(1,2),sum,na.rm=T)
  tmp=ncvar_get(a,"ECL_WS10")
  ws=apply(tmp[,,I],c(1,2),sum,na.rm=T)
  tmp=ncvar_get(a,"ECL_SLP")
  slp=apply(tmp[,,I],c(1,2),sum,na.rm=T)
  tmp=ncvar_get(a,"ECL_PRCP")
  prcp=apply(tmp[,,I],c(1,2),sum,na.rm=T)/6
  count=apply(!is.na(tmp[,,I]),c(1,2),sum,na.rm=T)
 } else {
  tmp=ncvar_get(a,"ECL_U10")
  u=u+apply(tmp[,,I],c(1,2),sum,na.rm=T)
  tmp=ncvar_get(a,"ECL_V10")
  v=v+apply(tmp[,,I],c(1,2),sum,na.rm=T)
  tmp=ncvar_get(a,"ECL_WS10")
  ws=ws+apply(tmp[,,I],c(1,2),sum,na.rm=T)
  tmp=ncvar_get(a,"ECL_SLP")
  slp=slp+apply(tmp[,,I],c(1,2),sum,na.rm=T)
  tmp=ncvar_get(a,"ECL_PRCP")
  prcp=prcp+apply(tmp[,,I],c(1,2),sum,na.rm=T)/6
  count=count+apply(!is.na(tmp[,,I]),c(1,2),sum,na.rm=T)
 }
 }
 u=u/count
 v=v/count
 ws=ws/count
 slp=slp/count
 prcp=prcp/count
} else {
a=nc_open(paste0(dir,"ECLcomposite_",reanal,"_proj240_rad2cv1_aust_20112015.nc"))
tmp=ncvar_get(a,"ECL_U10")
u=apply(tmp[,,I],c(1,2),mean,na.rm=T)
tmp=ncvar_get(a,"ECL_V10")
v=apply(tmp[,,I],c(1,2),mean,na.rm=T)
tmp=ncvar_get(a,"ECL_WS10")
ws=apply(tmp[,,I],c(1,2),mean,na.rm=T)
tmp=ncvar_get(a,"ECL_SLP")
slp=apply(tmp[,,I],c(1,2),mean,na.rm=T)
tmp=ncvar_get(a,"ECL_PRCP")
prcp=apply(tmp[,,I],c(1,2),mean,na.rm=T)/6
}

u[u>8]=8
u[u<=-8]=-8
v[v>8]=8
v[v<=-8]=-8
prcp[prcp>=3]=3
ws[ws>=11]=11

image.plot(lat,lon,u,breaks=seq(-8,8,1),col=pal1(16),zlim=c(-8,8),
      main="a) Zonal wind (m/s)",xlab="",ylab="",border="black",cex.main=1.5,
      lab.breaks=c("-8","","-6","","-4","","-2","","0","","+2","","+4","","+6","","+8"))
contour(lat,lon,u,levels=seq(-8,-1,1),lty=2,add=T,labcex=1,
        labels=c("8","","6","","4","","2",""))
contour(lat,lon,u,levels=seq(0,8,1),add=T,labcex=1,
        labels=c("0","","2","","4","","6","","8"))
points(0,0,pch=4,cex=2,lwd=2)
box()
image.plot(lat,lon,v,breaks=seq(-8,8,1),col=pal1(16),zlim=c(-8,8),
      main="b) Meridional wind (m/s)",xlab="",ylab="",border="black",cex.main=1.5,
      lab.breaks=c("-8","","-6","","-4","","-2","","0","","+2","","+4","","+6","","+8"))
contour(lat,lon,v,levels=seq(-8,-1,1),lty=2,add=T,labcex=1,
        labels=c("8","","6","","4","","2",""))
contour(lat,lon,v,levels=seq(0,8,1),add=T,labcex=1,
        labels=c("0","","2","","4","","6","","8"))
points(0,0,pch=4,cex=2,lwd=2)
box()
image.plot(lat,lon,ws,breaks=seq(0,11,1),col=pal2(11),zlim=c(0,11),
      main="c) Wind speed (m/s)",xlab="",ylab="",border="black",cex.main=1.5)
contour(lat,lon,ws,levels=seq(0,11,1),add=T,labcex=1,
        labels=c("0","","2","","4","","6","","8","","10",""))
points(0,0,pch=4,cex=2,lwd=2)
box()
image.plot(lat,lon,prcp,breaks=seq(0,3,0.25),col=pal2(12),zlim=c(0,2),
      main="d) Rain rate (mm/hr)",xlab="",ylab="",border="black",cex.main=1.5)
points(0,0,pch=4,cex=2,lwd=2)
contour(lat,lon,slp,levels=seq(990,1020,2),add=T,labcex=1)
box()
dev.off()
