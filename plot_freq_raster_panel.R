## This is code that takes the 11 ensemble members of a POAMA run
## Started on a particular day
## And outputs some interested results for further analysis

## I care about ECLs, so main result is a daily timeseries of # of members that
## have an ECL on that day

## Also outputs a figure for each of the months in the file
## That shows the number of cyclones centred in that month
## Option to set a higher intensity threshold
library(ncdf4)
library(raster)
library(maps)
library(abind)
library(sp)
color.palette <- function(steps, n.steps.between=NULL, ...){
  
  if(is.null(n.steps.between)) n.steps.between <- rep(0, (length(steps)-1))
  if(length(n.steps.between) != length(steps)-1) stop("Must have one less n.steps.between value than steps")
  
  fill.steps <- cumsum(rep(1, length(steps))+c(0,n.steps.between))
  RGB <- matrix(NA, nrow=3, ncol=fill.steps[length(fill.steps)])
  RGB[,fill.steps] <- col2rgb(steps)
  
  for(i in which(n.steps.between>0)){
    col.start=RGB[,fill.steps[i]]
    col.end=RGB[,fill.steps[i+1]]
    for(j in seq(3)){
      vals <- seq(col.start[j], col.end[j], length.out=n.steps.between[i]+2)[2:(2+n.steps.between[i]-1)]  
      RGB[j,(fill.steps[i]+1):(fill.steps[i+1]-1)] <- vals
    }
  }
  
  new.steps <- rgb(RGB[1,], RGB[2,], RGB[3,], maxColorValue = 255)
  pal <- colorRampPalette(new.steps, ...)
  return(pal)
}

col_anom <- color.palette(c("darkblue","blue","white","red","darkred"),c(10,20,20,10))
col_val <- color.palette(c("white","blue","darkblue","black"),c(20,10,5))

ColorBar <- function(brks,cols,vert=T,subsampleg=1)
{
  if(vert) {
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols,
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(1.5, length(brks) - 1.5, subsampleg), tick = TRUE,
       labels = brks[seq(2, length(brks)-1, subsampleg)])
  } else {
    par(mar = c(1.5, 1, 1, 1), mgp = c(1.5, 0.3, 0), las = 1, cex = 1)
    image(1:length(cols), 1, t(t(1:length(cols))), axes = FALSE, col = cols,
          xlab = '', ylab = '')
    box()
    axis(1, at = seq(1.5, length(brks) - 1.5, subsampleg),
         labels = brks[seq(2, length(brks)-1, subsampleg)])
  }
}

lat=seq(-89.5,89.5)
lon=seq(0,359.5)
makesmooth<-function(data,winwid=5)
{
  a=dim(data)
  m1=abind(data[(a[1]-winwid):a[1],a[2]:1],data[,a[2]:1],data[1:winwid,a[2]:1],along=1)
  m2=raster(t(m1),xmn=min(lon)-winwid,xmx=max(lon)+winwid,ymn=min(lat),ymx=max(lat))
  w2=focalWeight(m2,winwid,type="circle")
  m3=focal(m2,w2)
  tmp=t(as.matrix(m3)[a[2]:1,(winwid+1):(a[1]+winwid)])
  return(tmp)
}


plot_freq_panel<-function(year1,year2,seasons=rbind(c(5,10),c(11,4)),snames=c("MJJASO","NDJFMA"),
       dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf/",type2="_500km",
       type="cyclone",proj="proj100_rad5cv0.15",fout="output")
{

years=seq(year1,year2,1)

if(year1<1979) startN=2 else startN=1

reanals=c("ERAI","NCEP1","20CR")
varnames=c("systems","systems","EnsMean")
fend=c(".nc",".nc","_EnsMean.nc")

lat=seq(-89.5,89.5)
lon=seq(0,359.5)  ### Can always combine into bigger cells later

ss=length(snames)
breaks=c(0,0.05,seq(0.25,1,0.25),1.5,1000)
#breaks=c(0,0.05,seq(0.1,1,0.1),1000)
#col1=col_val(length(breaks)-1)
library(RColorBrewer)
col1=brewer.pal(length(breaks)-1,"Blues")

pnum=1
if(year1<1979)
{
pdf(file=fout,width=4*ss,height=6.3)
layout(rbind(seq(1,ss),seq(ss+1,ss*2),rep(ss*2+1,ss)),height=c(1,1,0.3))
} else {
pdf(file=paste0(fout,".pdf"),width=4*ss,height=9)
layout(rbind(seq(1,ss),seq(ss+1,ss*2),seq(ss*2+1,ss*3),rep(ss*3+1,ss)),height=c(1,1,1,0.3))
}

par(mar=c(2,2,4,1))

for(n in startN:3)
{
  a=nc_open(paste0(dir,"/",reanals[n],"_UM_global",type,"s_",proj,type2,fend[n]))
  systems=ncvar_get(a,varnames[n])
  years2=ncvar_get(a,"year")

  I=which(years2%in%years)
  J=which(years%in%years2)
  freq=array(NaN,c(length(lon),length(lat),length(I)))  

  for(s in 1:ss)
  {
   if(seasons[s,2]>=seasons[s,1])
   {
    for(y in 1:length(I)) freq[,,y]=apply(systems[,,I[y],seasons[s,1]:seasons[s,2]],c(1,2),sum,na.rm=T)    
   } else {
    tmp=abind(systems[,,I[-length(I)],seasons[s,1]:12],systems[,,I[-1],1:seasons[s,2]],along=4)

    for(y in 1:(length(I)-1)) freq[,,y]=apply(tmp[,,y,],c(1,2),sum,na.rm=T)
   }

   meanfreq=apply(freq,c(1,2),mean,na.rm=T)
   meanfreq2=makesmooth(meanfreq)

   image(lon,lat,meanfreq,breaks=breaks,col=col1,xlab="",ylab="",
          main=paste0(letters[pnum],") ",snames[s]," mean ",type," frequency: ",reanals[n]))
   map('world2',add=T)   
#   contour(lon,lat,meanfreq2,levels=breaks[seq(3,length(breaks),2)],add=T,lwd=2,col="black",drawlabels=F)
   contour(lon,lat,meanfreq2,levels=c(0.05,0.5,1,1.5),add=T,lwd=1.5,col="black",drawlabels=F)
   pnum=pnum+1
  }
}

ColorBar(breaks,col1,subsampleg=1,vert=F)
dev.off()
}
 
 
#plot_freq_panel(1980,2016,seasons=rbind(c(3,5),c(6,8),c(9,11),c(12,2)),snames=c("MAM","JJA","SON","DJF"),
#        dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf",
#        type="anticyclone",proj="proj100_rad10cv0.075",type2="_500km",
#        fout="paperfig_anticycfreq_3reanals_proj100_rad10cv0.075_500km_4seasons")

#plot_freq_panel(1980,2016,seasons=rbind(c(3,5),c(6,8),c(9,11),c(12,2)),snames=c("MAM","JJA","SON","DJF"),
#        dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf",
#        type="cyclone",proj="proj100_rad5cv0.15",type2="_500km",
#        fout="paperfig_cycfreq_3reanals_proj100_rad5cv0.15_500km_4seasons")

#plot_freq_panel(1980,2016,seasons=rbind(c(6,8),c(12,2)),snames=c("JJA","DJF"),
#        dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf",
#        type="anticyclone",proj="proj100_rad10cv0.075",type2="_500km",
#        fout="paperfig_anticycfreq_3reanals_proj100_rad10cv0.075_500km_jjadjf")

#plot_freq_panel(1980,2016,seasons=rbind(c(5,10),c(11,4)),snames=c("MJJASO","NDJFMA"),
#        dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf",
#        type="cyclone",proj="proj100_rad5cv0.15",type2="_500km",
#        fout="paperfig_cycfreq_3reanals_proj100_rad5cv0.15_500km_v2")


plot_freq_panel(1980,2016,seasons=rbind(c(5,10),c(11,4)),snames=c("MJJASO","NDJFMA"),
        dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf",
        type="anticyclone",proj="proj100_rad10cv0.075",type2="_500km",
        fout="paperfig_anticycfreq_3reanals_proj100_rad10cv0.075_500km_Bu")

#plot_freq_panel(1980,2016,seasons=rbind(c(5,10),c(11,4)),snames=c("MJJASO","NDJFMA"),
#        dir="/short/eg3/asp561/cts.dir/gcyc_out/netcdf",
#        type="cyclone",proj="proj100_rad5cv0.15",type2="_500km",
#        fout="paperfig_cycfreq_3reanals_proj100_rad5cv0.15_500km_v2")

