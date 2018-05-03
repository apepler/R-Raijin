source("init_cyclones.R")
library(ncdf4)
library(maps)
library(oz)
library(RColorBrewer)

years=1979:2015
dir=c("/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_highs_rad10cv0.075/",
      "/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_lows_rad5cv0.15/")
fname=c("austhighs_7deg_proj100_rad10cv0.075_500km_","austlows_6deg_proj100_rad5cv0.15_D2_")
ptype=c("Anticyclone","Cyclone")

for(f in 1:2)
for(y in 1:length(years))
{
a=nc_open(paste0(dir[f],fname[f],years[y],".nc"))

if(y==1 & f==1)
{
 time=ncvar_get(a,"time")
 lon=ncvar_get(a,"lon")
 lat=ncvar_get(a,"lat")
 data=array(0,c(length(lon),length(lat),2))
 data[,,f]=data[,,f]+apply(ncvar_get(a,"systems"),c(1,2),sum)

 } else {
 if(f==1) time=c(time,ncvar_get(a,"time"))
 data[,,f]=data[,,f]+apply(ncvar_get(a,"systems"),c(1,2),sum)
 }
}

## First : Number of hours

data2=data/length(years)
data2[data==0]=NaN
pdf(file="vicci_cycanticycfreq_19792015_aust2.pdf",width=7,height=3)
layout(cbind(1,2,3),width=c(1,1,0.3))
par(mar=c(2,2,4,1))
#breaks=seq(0,300,20)
breaks=c(seq(0,300,25),1000)
breaks2=seq(0,500,50)
#cols=col_val(length(breaks)-1)
cols=colorRampPalette(brewer.pal(9,"YlOrRd"))(length(breaks)-1)

for(n in 1:2)
{
 image(lon,lat,data2[,,n],breaks=breaks,col=cols,xlab="",ylab="",
          #xlim=c(140,151),ylim=c(-40,-32),
          xlim=c(110,155),ylim=c(-45,-10),
          main=paste0(ptype[n]," frequency (hours/year)"))
 contour(lon,lat,data2[,,n],levels=breaks2,add=T,xlim=c(140,151),ylim=c(-40,-32))
 map("world2",add=T)
 vic(add=T,lwd=2)
}
ColorBar(breaks,cols,subsampleg=1,vert=T)
dev.off()

## Second : percent of hours

data2=100*data/length(time)
data2[data==0]=NaN
pdf(file="vicci_cycanticycfreq_19792015_aust2_pc.pdf",width=7,height=3)
layout(cbind(1,2,3),width=c(1,1,0.3))
par(mar=c(2,2,4,1))
breaks=c(seq(0,20,2),100)
#cols=col_val(length(breaks)-1)
cols=colorRampPalette(brewer.pal(9,"YlOrRd"))(length(breaks)-1)

for(n in 1:2)
{
 image(lon,lat,data2[,,n],,breaks=breaks,col=cols,xlab="",ylab="",
          #xlim=c(140,151),ylim=c(-40,-32),
          xlim=c(110,155),ylim=c(-45,-10),
          main=paste0(ptype[n]," frequency (%)"))
 contour(lon,lat,data2[,,n],levels=breaks,add=T,xlim=c(140,151),ylim=c(-40,-32)) 
 map("world2",add=T)
 vic(add=T,lwd=2)
}
ColorBar(breaks,cols,subsampleg=1,vert=T)
dev.off()








