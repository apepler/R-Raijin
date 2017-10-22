for(cv in "075")
{
t1=Sys.time()
source("/short/eg3/asp561/R/init_cyclones.R")
library(abind)
library(verification)
years=1990:2012
members=paste("e",sprintf("%2.2d",1:11),sep="")

setwd("/short/eg3/asp561/cts.dir/gcyc_out/")
type="proj240_highs_rad10cv0.075"
type2="globalanticyclones_proj240_rad10cv0.075.nc"

a=nc_open(paste0("ERAI/daily_",type,"/ERAIdaily_UM_",type2))
erai=ncvar_get(a,"systems")
lat=ncvar_get(a,"lat")
lon=ncvar_get(a,"lon")

a=nc_open(paste("access-s1/",type,"/globalanticyclones_4monthlead_day1_cv0.",cv,".nc",sep=""))

erai2<-array(0,c(length(lon),length(lat),dim(erai)[3:4]))

for(i in 1:length(lon))
  for(j in 1:length(lat))
  {
    I=which((lon>=lon[i]-5 & lon<lon[i]+5) | lon>=lon[i]+355 | lon<lon[i]-355)
    J=which(lat>=lat[j]-5 & lat<lat[j]+5)
    erai2[i,j,,]=apply(erai[I,J,,],c(3,4),sum)
  }

skill<-array(NaN,c(length(lon),length(lat),12,4,4))
dimnames(skill)[[5]]=c("ACCESS.mean","AnomalyCorrelation","PercentConsistent","BrierSkill")

for(m in 1:12)
{
  print(m)
  tmp=ncvar_get(a,"cyclones",start=c(1,1,1,m,1,1),
                count=c(length(lon),length(lat),length(years),1,11,4))

  for(i in 1:length(lon))
    for(j in 1:length(lat))
    {
      I=which((lon>=lon[i]-5 & lon<lon[i]+5) | lon>=lon[i]+355 | lon<lon[i]-355)
      J=which(lat>=lat[j]-5 & lat<lat[j]+5)
      tmp2=apply(tmp[I,J,,,],c(3,4,5),sum)
      
      if(m<=9) tmpE=erai2[i,j,1:23,m:(m+3)] else tmpE=abind(erai2[i,j,1:23,m:12],erai2[i,j,2:24,1:(m%%9)])
      med=apply(tmp2,c(2,3),median)
      
      skill[i,j,m,,1]=apply(tmp2,3,mean)
      for(x in 1:4) 
        if(mean(med[,x])>=1 & median(tmpE[,x])>=1)
           {
             skill[i,j,m,x,2]=cor(tmpE[,x],apply(tmp2[,,x],1,mean))
             
             above=rep(0,length(years))
             for(y in 1:length(years)) above[y]=mean(tmp2[y,,x]>med[,x])
             
             skill[i,j,m,x,3]=length(which((above>0.5 & tmpE[,x]>median(tmpE[,x])) | (above<0.5 & tmpE[,x]<median(tmpE[,x]))))/length(years)
             
             b=brier(as.numeric(tmpE[,x]>median(tmpE[,x])),above)
             skill[i,j,m,x,4]=b$ss
        }

    }
  rm(tmp)
}

Sys.time()
Sys.time()-t1

save(skill,file=paste("ACCESS_skill_proj240_highs_rad10cv",cv,".RData",sep=""))

###### Want a panel plot - bias, anomaly correlation

names=c("Frequency Bias","Anomaly Correlation","Percent Consistent","Brier Skill Score")
levlist=list(seq(-10,10,2),seq(0.35,0.8,0.05),
             seq(45,90,5),seq(-0.05,0.5,0.05))
lead=paste(c("01","02","03","04","05","06","07","08","09",10:12),"01",sep="")

for(lt in 1:4)
for(m in 1:12)
  if(m!=4)
{
  pdf(file=paste("ACCESS_skill_anticyclonescv",cv,"_",lead[m],"_mo",lt,".pdf",sep=""),width=10,height=6,pointsize=12)
  layout(cbind(c(1,2),c(5,6),c(3,4),c(7,8)),width=c(1,0.2,1,0.2))
  par(mar=c(3,2,3,1))
  for(n in 1:4)
  {
    lev=levlist[[n]]
    if(n==1) 
    {
      m2=m+lt-1
      if(m2<=12) tmp=skill[,,m,lt,1]-apply(erai2[,,1:23,m2],c(1,2),mean) else
        tmp=skill[,,m,lt,1]-apply(erai2[,,2:24,m2%%12],c(1,2),mean)
      
      tmp[tmp<min(lev)]=min(lev)
      tmp[tmp>max(lev)]=max(lev)
      
      cols=rev(col_anom(length(lev)-1))
      image(lon,lat,tmp,breaks=lev,col=cols,main=names[n],xlab="",ylab="")
    } else {
      tmp=skill[,,m,lt,n]
      tmp[tmp>max(lev)]=max(lev)
      if(n==3) tmp=tmp*100
      cols=col_val(length(lev)-1)
      image(lon,lat,tmp,breaks=lev,col=cols,main=names[n],xlab="",ylab="")
      
    }
    map("world2",add=T,wrap=c(0,360))
  }
  for(n in 1:4) {
    lev=levlist[[n]]
    if(n==1) cols=rev(col_anom(length(lev)-1)) else cols=col_val(length(lev)-1)
    ColorBar(lev,cols)
  }

  dev.off()
}
}

