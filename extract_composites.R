library(ncdf4)
dir="/short/eg3/asp561/cts.dir/gcyc_out/ERAI/proj100_highs_rad10cv0.075/"

latrange=cbind(c(-40,-35),c(-45,-30),c(35,40),c(30,45),c(45,70))
latname=c("35-40S","30-45S","35-40N","30-45N","45-70N")

lat<-lon<-seq(-15,15,0.75)
years=seq(1980,2016)

statnames=c("Count","MSLP","CV","Depth","Radius")
stats=array(NaN,c(length(years),length(latname),length(statnames)))
dimnames(stats)[[1]]=years
dimnames(stats)[[2]]=latname
dimnames(stats)[[3]]=statnames

composites<-array(NaN,c(length(lon),length(lat),length(years),length(latname),2))
dimnames(composites)[[3]]=years
dimnames(composites)[[4]]=latname
dimnames(composites)[[5]]=c("slp","z500")


for(y in 1:length(years))
{
fixes=read.csv(paste0(dir,"ERAI_UM_highs_ERAI_proj100_rad10cv0.075_500km_",years[y],".csv"))
a=nc_open(paste0(dir,"ERAI_UM_highs_ERAI_proj100_rad10cv0.075_500km_",years[y],"_slp.nc"))
slp=ncvar_get(a,"psl")

for(r in 1:length(latrange[1,]))
{
I=which(fixes$Lat>=latrange[1,r] & fixes$Lat<=latrange[2,r])
stats[y,r,1]=length(I)
stats[y,r,2:5]=apply(fixes[I,9:12],2,mean)
composites[,,y,r,1]=apply(slp[,,I],c(1,2),mean)
}
rm(slp)

z=ncvar_get(a,"z500")
for(r in 1:length(latrange[1,]))
{
  I=which(fixes$Lat>=latrange[1,r] & fixes$Lat<=latrange[2,r])
  composites[,,y,r,2]=apply(z[,,I],c(1,2),mean)
}
rm(z)
}

save(lat,lon,years,latrange,stats,composites,file=paste0(dir,"ERAI_UM_highs_ERAI_proj100_rad10cv0.075_500km_annualcomposites.RData"))


