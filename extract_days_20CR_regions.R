years=seq(1851,2014)
library(ncdf4)
library(sp)

dir="/short/eg3/asp561/cts.dir/gcyc_out/"
type="proj100_highs_rad10cv0.075/"

regnames=c("Aust","SAust","SEA","SWWA","Tasman")
outfile="UM_20CR_anticyclonedays_austregions2_proj100_rad10cv0.075_500km"

daylist=seq.Date(as.Date(paste0(min(years),"0101"),format="%Y%m%d",tz="UTC"),
                 as.Date(paste0(max(years),"1231"),format="%Y%m%d",tz="UTC"),by="1 day")

cyctime=data.frame(Date=as.numeric(format(daylist,"%Y%m%d")),
                   AWAP.Date=as.numeric(format(daylist+1,"%Y%m%d")),
                   matrix(0,length(daylist),length(regnames)))
colnames(cyctime)=c("Date","AWAP.Date",regnames)

cycdays=array(0,c(length(daylist),length(regnames),56))
dimnames(cycdays)[[2]]=regnames
dimnames(cycdays)[[3]]=1:56

for(y in 1:length(years))
  for(n in 1:56)
  {
    print(paste(years[y],n))
    fname=paste(dir,"/20CR/",type,"/tracks_",years[y],"_",n,".dat",sep="")
    read.table(fname, sep="",skip=1)->fixes
    colnames(fixes)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV","Depth","Radius","Up","Vp")
    fixes=data.frame(fixes)
    fixes$Year=floor(fixes$Date/10000)
    fixes$CV=abs(fixes$CV)
    fixes$Depth=abs(fixes$Depth)
    if(length(unique(fixes$Year))>1) fixes=fixes[fixes$Year==unique(fixes$Year)[2],]
    fixes$Date=(fixes$Date%%10000) + years[y]*10000

    fixes$Closed=0
    fixes$Closed[fixes$Open==0 | fixes$Open==10]=1

    fixes$Tasman<-fixes$SWWA<-fixes$SEA<-fixes$SAust<-fixes$Aust<-0
    I<-which(fixes$Lon>=110 & fixes$Lon<=155 & fixes$Lat<(-10) & fixes$Lat>=-45)
    fixes$Aust[I]<-1
    I<-which(fixes$Lon>=110 & fixes$Lon<=155 & fixes$Lat<=(-25) & fixes$Lat>=-40)
    fixes$SAust[I]<-1
    I<-which(fixes$Lon>=135 & fixes$Lon<=147.5 & fixes$Lat<=(-32.5) & fixes$Lat>=-42.5)
    fixes$SEA[I]<-1
    I<-which(fixes$Lon>=110 & fixes$Lon<=122.5 & fixes$Lat<=(-27.5) & fixes$Lat>=-40)
    fixes$SWWA[I]<-1
    I<-which(fixes$Lon>=147.5 & fixes$Lon<=160 & fixes$Lat<(-25) & fixes$Lat>=-40)
    fixes$Tasman[I]<-1
 
    x<-rle(fixes$ID)
    events<-data.frame(ID=x$values,Length=x$lengths,Date1=rep(0,length(x$values)),Move=rep(0,length(x$values)),Closed=rep(0,length(x$values)))
    enumcols=dim(events)[2]

    for(i in 1:length(events$ID)) 
    {
    events$Date1[i]=min(fixes$Date[fixes$ID==events$ID[i]])
    I=which(fixes$ID==events[i,1])
    events$Move[i]=spDistsN1(cbind(fixes$Lon[min(I)],fixes$Lat[min(I)]),
                           cbind(fixes$Lon[max(I)],fixes$Lat[max(I)]),longlat=T)
    events$Closed[i]=max(fixes$Closed[I])
    }
 
    events=events[events$Move>=500 & events$Closed==1,]
    include<-match(fixes[,1],events[,1])
    J<-which(is.na(include)==0)
    fixes=fixes[J,]

   ## Now, make into a list of cyclone days

    for(r in 1:length(regnames))
    {
      I=which(fixes[,regnames[r]]==1)
      J=which(cyctime$Date%in%fixes$Date[I])
      cycdays[J,r,n]=1
    } 
  }

save(daylist,regnames,cycdays,cyctime,file=paste(dir,outfile,".RData",sep=""))
cyctime[,2+(1:length(regnames))]=apply(cycdays,c(1,2),sum) # Total number with an anticyclone
write.csv(cyctime,paste(dir,outfile,".csv",sep=""))
