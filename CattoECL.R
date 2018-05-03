library(ncdf4)

indir="/g/data/eg3/ajd548/vicci/cyclone_data_JC/"
years=seq(1979,2015,1)

a=nc_open(paste0(indir,"cyclones_2015_0.75.nc"))
lat=ncvar_get(a,"latitude")
lon=ncvar_get(a,"longitude")

#ECLmask<-matrix(0,length(lon),length(lat))
#I=which(lon>=149 & lon<=161)
#J=which(lat>=-41 & lat<=-37)
#ECLmask[I,J]<-1
#I=which(lon>=152 & lon<=161)
#J=which(lat>=-31 & lat<=-24)
#ECLmask[I,J]<-1
#J=which(lat>=-37 & lat<=-31)
#for(j in 1:length(J))
#{
#  I<-which(lon>=(149+(37+lat[J[j]])/2) & lon<=161)
#  ECLmask[I,J[j]]=1
#}

ECLmask<-matrix(0,length(lon),length(lat))
I=which(lon>=152 & lon<=158)
J=which(lat>=-38 & lat<=-37)
ECLmask[I,J]<-1
I=which(lon>=155 & lon<=158)
J=which(lat>=-31 & lat<=-27)
ECLmask[I,J]<-1
J=which(lat>=-37 & lat<=-31)
for(j in 1:length(J))
{
  I<-which(lon>=(152+(37+lat[J[j]])/2) & lon<=158)
  ECLmask[I,J[j]]=1
}

for(year in years)
{
print(year)
print(min(years))

a=nc_open(paste0(indir,"cyclones_",year,"_0.75.nc"))
time=seq.POSIXt(as.POSIXct(paste0(year,"0101 00:00"),format="%Y%m%d %H:%M",tz="UTC"),
                 as.POSIXct(paste0(year,"1231 18:00"),format="%Y%m%d %H:%M",tz="UTC"),by="6 hours")

#time=as.POSIXlt(ncvar_get(a,"time")*60*60,origin="1900-01-01 00:00",tz="UTC")
cat_c=ncvar_get(a,"FLAG")

days=as.integer(format(time,"%Y%m%d"))
tmp=unique(days)
tmp2<-data.frame(Day=tmp,ECL=rep(0,length(tmp)))
for(i in 1:length(tmp2$Day))
{
  I=which(days==tmp2$Day[i])
  tmp2$ECL[i]=max(apply(cat_c[,,I],c(1,2),max)*ECLmask)
}

if(year==min(years)) ECLs=tmp2 else ECLs=rbind(ECLs,tmp2)

}

write.csv(ECLs,"Catto_ECLs2.csv")
