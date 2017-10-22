rm(list=ls())
source("Q:/R/init_cyclones.R")
library(ncdf4)

years=1851:2016

types=c("Anticyclones","Cyclones")
dir="\\nina\acaciap\Cyclones\"

names=c("NCEP1","ERAI","20CR.ensemble","20CR.mean")
yearS=c(1950,1980,1851,1851)
yearE=c(2016,2016,2014,2014)

fnames=cbind(c("NCEP1_UM_globalanticyclones_proj100_rad10cv0.075_D2.nc",
               "ERAI_UM_globalanticyclones_proj100_rad10cv0.075_D2.nc",
               "20CR_UM_anticyclonelats_proj100_rad10cv0.075_D2.nc"),
            c("NCEP1_UM_globalcyclones_proj100_rad5cv0.15_D2.nc",
               "ERAI_UM_globalcyclones_proj100_rad5cv0.15_D2.nc",
               "20CR_UM_cyclonelats_proj100_rad5cv0.15_D2.nc"))
