### SOURCE THIS AT START
### Lots of useful functions & base libraries

library(sp)
#library(fields)
library(abind)
library(ncdf4)


source("/short/eg3/asp561/R/color.palette.R")
col_anom <- color.palette(c("darkred","red","white","blue","darkblue"),c(10,20,20,10))
col_val <- color.palette(c("white","blue","darkblue","black"),c(20,10,5))

ColorBar <- function(brks,cols,vert=T,subsampleg=1)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 0.8)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(1.5, length(brks) - 1.5, subsampleg), tick = TRUE, 
       labels = brks[seq(2, length(brks)-1, subsampleg)])
}

### Compares the PDF of two datasets
makePDF = function(data1,data2,xlabel="Intensity",labloc="topright",leg=c("Data1","Data2"),tit="") {
  a=density(data1,na.rm=T)
  b=density(data2,na.rm=T)
  
  lims=range(data1,data2,na.rm=T)
  if((lims[2]-lims[1])<10)
  {
    lims[1]=floor(lims[1])
    lims[2]=ceiling(lims[2])
  } else {
    lims[1]=floor(lims[1]/5)*5
    lims[2]=ceiling(lims[2]/5)*5
  }
  
  plot(a,col=rgb(0,0,1,1/4),xlim=lims,ylim=range(0,a$y,b$y),
       xlab=xlabel,ylab="Frequency",main=tit)
  polygon(a,col=rgb(0,0,1,1/4),density=-1)
  polygon(b,col=rgb(1,0,0,1/4),density=-1)
  legend(labloc,legend=leg,
         col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)),lwd=4,cex=1,bty="n")   
  
  return(ks.test(data1,data2))
}


## Simple match script
matchlows<-function(src,comp,dist=500)
{
  src$Match.ID<-src$Match<-0
  for(i in 1:length(src$ID))
  {
    I=which(comp$Date2==src$Date2[i])
    if(length(I)>0)
    {
      tmp=spDistsN1(cbind(comp$Lon[I],comp$Lat[I]),
                    c(src$Lon[i],src$Lat[i]),longlat=T)
      if(min(tmp)<dist)
      {
        k=which(tmp==min(tmp))
        src$Match[i]=I[k]
        src$Match.ID[i]=comp$ID[I[k]]
      }
      
    }
  }
  return(src)
}

### Script to load & format an ECL fix file

load_fixes<-function(fname,ECL=F,single=F,year=NaN)
{
  if(single) {
  read.table(fname, sep="",skip=1)->tmp
  tmp=tmp[,1:9]
  colnames(tmp)=c("ID","Fix","Date","Time","Open","Lon","Lat","MSLP","CV") 
  tmp$Year=floor(tmp$Date/10000)
  a=unique(tmp$Year)
  if(length(a)>1) tmp=tmp[tmp$Year==a[2],]
  if(!is.na(year)) tmp$Date=tmp$Date%%10000+year*10000  
  } else {
    read.csv(fname)->tmp
  }

  if(!is.na(year)) tmp$Date=tmp$Date%%10000+year*10000
  tmp$Year=floor(tmp$Date/10000)
  tmp$Month=floor(tmp$Date/100)%%100
  
  tmp$Location=0
  if(ECL)
  {
    I<-which(tmp$Lon>=149 & tmp$Lon<=161 & tmp$Lat<(-37) & tmp$Lat>=-41)
    tmp$Location[I]<-1
    I<-which(tmp$Lon>=(149+(37+tmp$Lat)/2) & tmp$Lon<=161 & tmp$Lat<(-31) & tmp$Lat>=-37)
    tmp$Location[I]<-1
    I<-which(tmp$Lon>=152 & tmp$Lon<=161 & tmp$Lat<=(-24) & tmp$Lat>=-31)
    tmp$Location[I]<-1
  } else 
    I=which(tmp$Lat>=-45 & tmp$Lat<=-10 & tmp$Lon>=110 & tmp$Lon<=155)
  
  tmp$Location[I]=1
  tmp$CV=abs(tmp$CV)
  tmp$Date2=as.POSIXct(paste(tmp$Date,tmp$Time,sep=""),
                       format="%Y%m%d %H:%M",tz="GMT")
  return(tmp)
}

make_events<-function(fixes)
{
  x<-rle(fixes$ID)
  
  ##Colnames = ID, length, max cv, max mslp, ever aust, max cv, max mslp
  events<-cbind(x$values,x$lengths,matrix(data=NaN,nrow=length(x$values),ncol=7))
  colnames(events)=c("ID","Length","Date1","Date2","MaxCV","MaxMSLP",
                   "Aust","Aust.CV","Aust.MSLP")
  
  I<-which(events[,2]==1)
  y<-match(events[I,1],fixes$ID)
  events[I,3]<-events[I,4]<-fixes$Date[y]
  events[I,5]<-fixes$CV[y]
  events[I,6]<-fixes$MSLP[y]
  events[I,7]<-fixes$Location[y]
  J=which(events[I,7]==1)
  events[I[J],8]=fixes$CV[y[J]]
  events[I[J],9]=fixes$MSLP[y[J]]
  
  ##Multi-fix events
  J<-which(events[,2]>1)
  
  for(i in 1:length(J)) 
  {
    I<-which(fixes[,1]==events[J[i],1])
    events[J[i],3]=min(fixes$Date[I]) ##Date1
    events[J[i],4]=max(fixes$Date[I]) ##Date2
    events[J[i],5]=max(fixes$CV[I]) ##Min central pressure
    events[J[i],6]=max(fixes$MSLP[I]) ##Max curvature
    if(sum(fixes$Location[I])>0) 
    {
      events[J[i],7]=1 ##Ever in ECL region
      events[J[i],8]=max(fixes$CV[I]*fixes$Location[I]) # All outside region = 0
      events[J[i],9]=max(fixes$MSLP[I]*fixes$Location[I]) # All outside region = 0
    } else events[J[i],7]=0
    
  }
  return(as.data.frame(events))
}

