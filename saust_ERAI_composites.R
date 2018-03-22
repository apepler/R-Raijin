library(ncdf4)
source("make_colours.R")

dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_highs_rad10cv0.075/"

snames=c("MAM","JJA","SON","DJF","MJJASO","NDJFMA","Annual")
mlist=list(3:5,6:8,9:11,c(1:2,12),5:10,c(1:4,11:12),1:12)
regnames=c("Aust","SAust","SEA","SWWA","Tasman")
breaks=c(seq(0,5,0.5),100)
cols=col_val(length(breaks)-1)

## Get and save data

fixes=read.csv(paste0(dir,"UM_highs_ERAI_proj100_rad10cv0.075_fixes_500km.csv"),stringsAsFactors=F)
fixes=fixes[fixes$Aust==1,]
fixes$Month=floor(fixs$Date/100)%%100

a=nc_open(paste0(dir,"UM_highs_ERAI_proj100_rad10cv0.075_austrain_500km.nc"))
lat=ncvar_get(a,"lat")
lon=ncvar_get(a,"lon")
cycrain=ncvar_get(a,"PRCP")
cycslp=ncvar_get(a,"SLP")

for(s in 1:7)
for(r in 1:length(regnames)) 
{
S=which(fixes$Month%in%mlist[[s]] & fixes[,regnames[r]]==1])
avrain=apply(cycrain[,,S],c(1,2),mean)
avslp=apply(cycslp[,,S],c(1,2),mean)

## Make some plots

pdf(file=paste0("UM_highs_ERAI_proj100_rad10cv0.075_",regnames[r],"rain_500km_",snames[s],".pdf"),height=4,width=5)
par(mar=c(2,2,1,1))
layout(cbind(1,2),width=c(1,0.25))
image(lon,lat,avrain[,length(lat):1,r],breaks=breaks,col=cols, # Have to reverse because ERAI is upside down
      xlab="",ylab="",main="")
contour(lon,lat,avslp[,length(lat):1,r],levels=seq(980,1040,4),lwd=2,col="black",add=T)
ColorBar(breaks,cols,subsampleg=1,vert=T)
dev.off()

}



