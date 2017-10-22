years=seq(1851,2014)
library(ncdf4)
library(sp)

dir="/short/eg3/asp561/cts.dir/gcyc_out/20CR/"
type="proj100_highs_rad10cv0.075/"
    
outfile="UM_20CR_anticyclonestats_proj100_500km_subregions_long.RData"
hilo="high"
cvthresh=0.15
slpthresh=1020
durthresh=2
move=500

regnames=c("Aust","SAust","SEA","SWWA","Tas","ESB")
stypes=c(paste(regnames,"Hours",sep="."),paste(regnames,".Hours.CV",cvthresh,sep=""),paste(regnames,".Hours.SLP",slpthresh,sep=""))

stats=array(NaN,c(length(years),12,length(stypes),57))
dimnames(stats)[[1]]=years
dimnames(stats)[[2]]=1:12
dimnames(stats)[[3]]=stypes
dimnames(stats)[[4]]=c(1:56,"EnsMean")


for(y in 1:length(years))
{
  for(n in 1:56)
  {
    print(paste(years[y],n))
    fname=paste(dir,type,"/tracks_",years[y],"_",n,".dat",sep="")
    read.table(fname, sep="",skip=1)->fixes
    colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Depth","Radius","Up","Vp")
    fixes=data.frame(fixes)
    fixes$Year=floor(fixes$Date/10000)
    fixes$CV=abs(fixes$CV)
    fixes$Depth=abs(fixes$Depth)
    if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
    fixes$Date=(fixes$Date%%10000) + years[y]*10000
    fixes$Month=floor(fixes$Date/100)%%100
    numcols=dim(fixes)[2]
 
fixes$ESB<-fixes$Tas<-fixes$SWWA<-fixes$SEA<-fixes$Saus<-fixes$Location<-0
    I<-which(fixes$Lon>=110 & fixes$Lon<=155 & fixes$Lat<(-10) & fixes$Lat>=-45)
    fixes$Location[I]<-1
I<-which(fixes$Lon>=110 & fixes$Lon<=155 & fixes$Lat<=(-30) & fixes$Lat>=-40)
fixes$Saus[I]<-1
I<-which(fixes$Lon>=135 & fixes$Lon<=147.5 & fixes$Lat<=(-32.5) & fixes$Lat>=-42.5)
fixes$SEA[I]<-1
I<-which(fixes$Lon>=110 & fixes$Lon<=122.5 & fixes$Lat<=(-27.5) & fixes$Lat>=-40)
fixes$SWWA[I]<-1
I<-which(fixes$Lon>=140 & fixes$Lon<=150 & fixes$Lat<=(-40) & fixes$Lat>=-50)
fixes$Tas[I]<-1
I<-which(fixes$Lon>=147.5 & fixes$Lon<=155 & fixes$Lat<(-25) & fixes$Lat>=-40)
fixes$ESB[I]<-1

    x<-rle(fixes$ID)
    events<-data.frame(ID=x$values,Length=x$lengths,Location=rep(0,length(x$values),Move=rep(0,length(x$values))))
    for(i in 1:length(events$ID)) 
     {
     events$Location[i]=sum(fixes$Location[fixes$ID==events$ID[i]])
     I=which(fixes$ID==events[i,1])
     events$Move[i]=spDistsN1(cbind(fixes$Lon[min(I)],fixes$Lat[min(I)]),
                           cbind(fixes$Lon[max(I)],fixes$Lat[max(I)]),longlat=T)
     }
    events=events[events$Move>=move & events$Location>=1,]
    include<-match(fixes[,1],events[,1])
    J<-which(is.na(include)==0)
    fixes=fixes[J,]
    
    for(m in 1:12)
      for(x in 1:length(regnames))
      {
        J=which(fixes[,numcols+x]==1 & fixes$Month==m)
        stats[y,m,x,n]=length(J)
        J=which(fixes[,numcols+x]==1 & fixes$CV>=cvthresh & fixes$Month==m)
        stats[y,m,x+length(regnames),n]=length(J)
        if(hilo=="high") J=which(fixes[,numcols+x]==1 & fixes$MSLP>=slpthresh & fixes$Month==m) else J=which(fixes[,numcols+x]==1 & fixes$MSLP<=slpthresh & fixes$Month==m)
        stats[y,m,x+(2*length(regnames)),n]=length(J)

      }
      
  }
  n=57
  fname=paste(dir,"EnsMean/",type,"/tracks_",years[y],".dat",sep="")
  read.table(fname, sep="",skip=1)->fixes
  colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Depth","Radius","Up","Vp")
  fixes=data.frame(fixes)
  fixes$Year=floor(fixes$Date/10000)
  fixes$CV=abs(fixes$CV)
  fixes$Depth=abs(fixes$Depth)
  if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
  fixes$Date=(fixes$Date%%10000) + years[y]*10000  
  fixes$Month=floor(fixes$Date/100)%%100

    numcols=dim(fixes)[2]

fixes$ESB<-fixes$Tas<-fixes$SWWA<-fixes$SEA<-fixes$Saus<-fixes$Location<-0
    I<-which(fixes$Lon>=110 & fixes$Lon<=155 & fixes$Lat<(-10) & fixes$Lat>=-45)
    fixes$Location[I]<-1
I<-which(fixes$Lon>=110 & fixes$Lon<=155 & fixes$Lat<=(-30) & fixes$Lat>=-40)
fixes$Saus[I]<-1
I<-which(fixes$Lon>=135 & fixes$Lon<=147.5 & fixes$Lat<=(-32.5) & fixes$Lat>=-42.5)
fixes$SEA[I]<-1
I<-which(fixes$Lon>=110 & fixes$Lon<=122.5 & fixes$Lat<=(-27.5) & fixes$Lat>=-40)
fixes$SWWA[I]<-1
I<-which(fixes$Lon>=140 & fixes$Lon<=150 & fixes$Lat<=(-40) & fixes$Lat>=-50)
fixes$Tas[I]<-1
I<-which(fixes$Lon>=147.5 & fixes$Lon<=155 & fixes$Lat<(-25) & fixes$Lat>=-40)
fixes$ESB[I]<-1

    x<-rle(fixes$ID)
    events<-data.frame(ID=x$values,Length=x$lengths,Location=rep(0,length(x$values),Location2=rep(0,length(x$values))))
    for(i in 1:length(events$ID)) events$Location[i]=sum(fixes$Location[fixes$ID==events$ID[i]])

    events=events[events$Length>=durthresh & events$Location>=1,]
    include<-match(fixes[,1],events[,1])
    J<-which(is.na(include)==0)
    fixes=fixes[J,]

    for(m in 1:12)
      for(x in 1:length(regnames))
      {
        J=which(fixes[,numcols+x]==1 & fixes$Month==m)
        stats[y,m,x,n]=length(J)
        J=which(fixes[,numcols+x]==1 & fixes$CV>=cvthresh & fixes$Month==m)
        stats[y,m,x+length(regnames),n]=length(J)
        if(hilo=="high") J=which(fixes[,numcols+x]==1 & fixes$MSLP>=slpthresh & fixes$Month==m) else J=which(fixes[,numcols+x]==1 & fixes$MSLP<=slpthresh & fixes$Month==m)
        stats[y,m,x+(2*length(regnames)),n]=length(J)

      }

   
    
}

##In case netcdf doesn't work
save(stats,file=paste(dir,outfile,sep=""))

