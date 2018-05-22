source("init_cyclones.R")
dirs=c("/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj240_lows_rad2cv1/","/short/eg3/asp561/cts.dir/gcyc_out/BARRA/proj240_lows_rad2cv1/")
reanals=c("ERAI","BARRA")
ylim=c(2011,2015)
pal2 <- color.palette(c("white","cyan","blue","black"), c(10,20,10))

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


pdf(file="TRMMrain_cyclone_panel_ERAI_BARRA_SEA_2.pdf",width=7,height=3,pointsize=12)
layout(cbind(1,2,3),width=c(1,1,0.35))
par(mar=c(2,2,3,2),cex.lab=1.2,cex.axis=1.2)
breaks=c(seq(0,3,0.25),10)
col=pal2(length(breaks)-1)

for(r in 1:2)
{
fixes=read.csv(paste0(dirs[r],"UM_lows_",reanals[r],"_proj240_rad2cv1_bigaust_fixes.csv"))
fixes=fixes[fixes$Date>=ylim[1]*10000 & fixes$Date<=(ylim[2]+1)*10000 & fixes$Lon>=110 & fixes$Lon<=160 & fixes$Lat>=-45 & fixes$Lat<=-10,]
print(dim(fixes))

fixes$Location=0
#fixes$Location[fixes$Lon>=110 & fixes$Lon<=155 & fixes$Lat>=-45 & fixes$Lat<=-10]=1


#    I<-which(fixes$Lon>=149 & fixes$Lon<=161 & fixes$Lat<(-37) & fixes$Lat>=-41)
#    fixes$Location[I]<-1
#    I<-which(fixes$Lon>=(149+(37+fixes$Lat)/2) & fixes$Lon<=161 & fixes$Lat<(-31) & fixes$Lat>=-37)
#    fixes$Location[I]<-1
#    I<-which(fixes$Lon>=152 & fixes$Lon<=161 & fixes$Lat<=(-24) & fixes$Lat>=-31)
#    fixes$Location[I]<-1
fixes$Location[fixes$Lon>=135 & fixes$Lon<=148 & fixes$Lat>=-42.5 & fixes$Lat<=-32.5]=1

I=which(fixes$Location==1)

a=nc_open(paste0(dirs[r],"ECLrain_TRMM_",reanals[r],"_proj240_rad2cv1_aust_20112015a.nc"))
lat<-lon<-seq(-10,10,length.out=81)
rain=ncvar_get(a,"ECL_PRCP")
print(dim(rain))
image(lon,lat,apply(rain[,,I],c(1,2),mean,na.rm=T),xlab="",ylab="",main=reanals[r],
   breaks=breaks,col=col)
points(0,0,pch=4,cex=2,lwd=2)
box()
}
ColorBar2(breaks,col,vert=T,subsampleg=2)
dev.off()

