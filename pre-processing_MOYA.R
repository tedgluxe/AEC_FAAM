################################################################################
### Airborne Eddy Covariance files preparation from FAAM aircraft data ###

#Contributions from: Adam Vaughan, Freya Squires, Will Drysdale, Dominika Pasternak.

################################################################################
### Loading packages ###

library(dplyr)
library(magrittr)
library(rgdal)
library(sp)
library(ncdf4)
library(ggplot2)
library(ggmap)
library(viridis)
library(lubridate)
library(stringr)
library(stringi)
library(openair)
library(reshape2)



################################################################################
### Initial setup ###

#set working directory
setwd("G:/Shared drives/eddy4R AEC/pre_processing")

#list all files
core_32_files <-  list.files("./aimms_32hz",pattern = ".nc") # core 32 Hz data
core_01_files <- list.files("./aimms_1hz",pattern = ".nc") # core 1 Hz data (things that are wonky in 32 Hz)
FGGA_10_files <- list.files("./fgga_10hz",pattern = ".na") # FGGA 10 Hz data (methane)
FGGA_H2O_files <- list.files("./fgga_h2o",pattern = ".txt") # FGGA raw 10 Hz file (water vapour)
thermistor_files <- list.files("./thermistor_16Hz",pattern = ".csv") # Thermistor 16 Hz files (temperature)

#choose flight and UTM zone
flno <- 5
zn <- "+proj=utm +zone=35" #UTM zones: Uganda 35, Zambia 36, Finland 35
bl <- 800 #boundary layer in m (see inversions in profiles & choose minimal value)

#get full flight number
FLIGHT_num <- stringr::str_sub(FGGA_10_files, start= -7) %>% gsub(".na","",.)
fn <- FLIGHT_num[flno] 

#start date
origin_all <-  stringr::str_sub(FGGA_10_files, start= -19, end=-12) %>% gsub(".na","",.)
origin <-  origin_all[flno]
origin <-  paste0(origin, " 00:00")

#frequencies
core_freq <-  32
fgga_freq <- 10
rad_freq <-  1
temp_freq <- 32



################################################################################
### Loading core 32 Hz data ###

#open 32 Hz file
core32 <- ncdf4::nc_open(paste0("./aimms_32hz/",core_32_files[grep(fn,core_32_files,ignore.case=TRUE)]))

#create a list of variables needed (exc. ones that don't load well from 32 Hz files)
eddyvars <- c("ROLL_GIN", "U_C", "V_C", "W_C", "LAT_GIN", "LON_GIN","HDG_GIN","TAS", "PS_RVSM")

#turn NetCDF into data frame
for (i in 1:length(eddyvars)) {
  vname <- eddyvars[i]
  raw <- as.vector(ncdf4::ncvar_get(core32,vname,collapse_degen=FALSE))
  if(i==1){
    CORE_32Hz <- data.frame(raw)
    names(CORE_32Hz) <- vname
  } 
  else {
    CORE_32Hz <- cbind(CORE_32Hz,raw)
    names(CORE_32Hz)[ncol(CORE_32Hz)] <- vname
  }
}

#get time and adjust frequency 
core_time <- ncvar_get(core32, attributes(core32$dim)$names[1]) %>% as.vector()
date <- strptime(x = origin, format ="%Y%m%d %H:%M") + (core_time)
CORE_32Hz$date <- base::as.POSIXct(seq.POSIXt(from = min(date)+(1/core_freq),
                                               to = max(date)+1, 
                                               by = 1/core_freq,
                                               tz = "UTC"))

#tidy up
rm(core32, date, core_time, i, raw, vname)



################################################################################
### Loading core 1 Hz data ###

#open 1 Hz file
core1 <- ncdf4::nc_open(paste0("./aimms_1hz/",core_01_files[grep(fn,core_01_files,ignore.case=TRUE)]))

#create a list of variables needed
eddyvars <- c("HGT_RADR")

#turn NetCDF into data frame
for (i in 1:length(eddyvars)) {
  vname <- eddyvars[i]
  raw <- as.vector(ncdf4::ncvar_get(core1,vname,collapse_degen=FALSE))
  if(i==1){
    CORE_1Hz <- data.frame(raw)
    names(CORE_1Hz) <- vname
  } 
  else {
    CORE_1Hz <- cbind(CORE_1Hz,raw)
    names(CORE_1Hz)[ncol(CORE_1Hz)] <- vname
  }
}

#get time
rad_time <- ncvar_get(core1, attributes(core1$dim)$names[1]) %>% as.vector()
date <- strptime(x = origin, format ="%Y%m%d %H:%M") + (rad_time)
CORE_1Hz$date <- base::as.POSIXct(date, tz = "UTC")

#tidy up
rm(core1, date, rad_time, i, raw, vname)



################################################################################
### Loading 16 Hz Thermistor data ###

TEMP_16Hz <- read.csv(file=paste0("./thermistor_16Hz/",thermistor_files[grep(fn,thermistor_files,ignore.case=TRUE)]),
                     header = T, 
                     sep = ",",
                     stringsAsFactors = F) %>%
  dplyr::rename(date=datetime, temp=TAT_THERM_U) %>%
  dplyr::mutate(date=as.POSIXct(strptime(date, "%Y-%m-%d %H:%M:%OS", tz="UTC"))) %>%
  dplyr::select(date,temp)



################################################################################
### Loading calibrated 10 Hz FGGA data ###

#open file
FGGA_10Hz <- read.delim(file=paste0("./fgga_10hz/",FGGA_10_files[grep(fn,FGGA_10_files,ignore.case=TRUE)]),
                        skip = 62, 
                        header = F, 
                        sep=" ") %>%
  dplyr::rename(.,date=V1,CO2_ppm=V2,CH4_ppb=V4) %>%
  dplyr::select(.,c(date,CO2_ppm,CH4_ppb))

#format date
FGGA_10Hz$date <- FGGA_10Hz$date + strptime(x = origin, format ="%Y%m%d %H:%M")

#NA flagged data and pressure issues (sanity check)
FGGA_10Hz$CH4_ppb[FGGA_10Hz$CH4_ppb < 0] <- NA 
FGGA_10Hz$CO2_ppm[FGGA_10Hz$CO2_ppm < 0] <- NA
FGGA_10Hz$CH4_ppb[FGGA_10Hz$CH4_ppb >9500] <- NA 
FGGA_10Hz$CO2_ppm[FGGA_10Hz$CO2_ppm >950] <- NA



################################################################################
### Loading raw 10 Hz FGGA data (water) ###

#open file
H2O_10Hz <- read.delim(file=paste0("./fgga_h2o/",FGGA_H2O_files[grep(fn,FGGA_H2O_files,ignore.case=TRUE)]),
                       header = T, 
                       sep = ",",
                       stringsAsFactors = F,
                       skipNul = TRUE) %>%
  dplyr::rename(date=Time) %>%
  dplyr::mutate(date=as.POSIXct(strptime(date, "%d/%m/%Y %H:%M:%OS", tz="UTC"))) %>%
  dplyr::select(date,H2O_ppm)

#NA calibration data and pressure issues
H2O_10Hz$H2O_ppm[H2O_10Hz$H2O_ppm < 1] <- NA



################################################################################
### Merging ###

#crop all to calibrated FGGA file length
CORE_1Hz <- CORE_1Hz %>% subset(.,date>=min(FGGA_10Hz$date) & date<=max(FGGA_10Hz$date))
CORE_32Hz <- CORE_32Hz %>% subset(.,date>=min(FGGA_10Hz$date) & date<=max(FGGA_10Hz$date))
H2O_10Hz <- H2O_10Hz %>% subset(.,date>=min(FGGA_10Hz$date) & date<=max(FGGA_10Hz$date))
TEMP_16Hz <-  TEMP_16Hz %>% subset(.,date>=min(FGGA_10Hz$date) & date<=max(FGGA_10Hz$date))

#interpolate up to 30 Hz
FL_mrg <- cbind(CORE_32Hz,
                sapply(names(CORE_1Hz)[-which(names(CORE_1Hz)=="date")],
                       function(x)
                         x <- approx(CORE_1Hz$date,CORE_1Hz[[x]],
                                     CORE_32Hz$date)$y)) %>%
          cbind(.,
                sapply(names(TEMP_16Hz)[-which(names(TEMP_16Hz)=="date")],
                       function(x)
                         x <- approx(TEMP_16Hz$date,TEMP_16Hz[[x]],
                                     CORE_32Hz$date)$y)) %>%
          cbind(.,
                sapply(names(FGGA_10Hz)[-which(names(FGGA_10Hz)=="date")],
                       function(x)
                         x <- approx(FGGA_10Hz$date,FGGA_10Hz[[x]],
                                     CORE_32Hz$date)$y)) %>%
          cbind(.,
                sapply(names(H2O_10Hz)[-which(names(H2O_10Hz)=="date")],
                       function(x)
                         x <- approx(H2O_10Hz$date,H2O_10Hz[[x]],
                                     CORE_32Hz$date)$y)) %>%
          na.omit()

#tidy up
rm(CORE_32Hz, CORE_1Hz, FGGA_10Hz, H2O_10Hz, TEMP_16Hz)



################################################################################
### Snipping legs ###

#filter data based on altitude and roll angle
FL_brk <- FL_mrg %>%
  dplyr::filter(HGT_RADR<610) %>%
  dplyr::filter(abs(ROLL_GIN)<20)

#find time gaps longer than 5 seconds / 500 m
t_dif <- diff(FL_brk$date) %>% as.numeric(.)
t_gaps <- which(abs(t_dif)>5)

#calculate start and end positions for each leg
t_start <- FL_brk$date[c(1,(t_gaps[1:(length(t_gaps)-1)]+1))]
t_end <- FL_brk$date[t_gaps]

#save the legs for post-processing
legs <- data.frame(t_start,t_end)
saveRDS(legs, file=paste0("./processed/",fn,"legs.rds"))



################################################################################
### Preparing input files ###

for(bk in 1:length(t_start)){
  
  #clip data for each leg
  sub_mrg <- FL_mrg %>% subset(date>=t_start[bk] & date <= t_end[bk])
  
  #calculate UTM coordinates
  d <- data.frame(lon=sub_mrg$LON_GIN, lat=sub_mrg$LAT_GIN)
  coordinates(d) <- c("lon", "lat")
  proj4string(d) <- CRS("+proj=longlat")
  cov <- spTransform(d, CRS(zn))
  
  #calculate ground distance covered by the aircraft in m
  stretch <- c(0,sqrt((diff(cov@coords[,1]))^2 + (diff(cov@coords[,2]))^2))
  stretch <- cumsum(stretch)
  parcel <- cumsum(sub_mrg$TAS * 1/2)
  parcel <- parcel - min(parcel, na.rm=TRUE)
  
  #calculate the volume of air passed by the aircraft
  sub_mrg$d_xy_travel <- stretch
  sub_mrg$d_xy_flow <- parcel
  sub_mrg$d_x_utm <- cov@coords[,1]
  sub_mrg$d_y_utm <- cov@coords[,2]
  rm(stretch,parcel,cov)
  
  #filter out any data using FAAM's -9999 system
  for(gg in names(sub_mrg)){
    x <- which(sub_mrg[[gg]]==(-9999))
    if(length(x)>0){sub_mrg[[gg]][x] <- NA}
    rm(x)
  }
  rm(gg)
  
  #sanity check for high altitudes and rolls
  sub_mrg <- sub_mrg %>%
    dplyr::filter(HGT_RADR<610.01) %>%
    dplyr::filter(abs(ROLL_GIN)<20.01) %>%
    na.omit()
  
  #distance travelled check
  sub_mrg$d_xy_travel <- sub_mrg$d_xy_travel - sub_mrg$d_xy_travel[1]
  sub_mrg$d_xy_flow <- sub_mrg$d_xy_flow - sub_mrg$d_xy_flow[1]
  
  #flag whether fulfils criteria for AEC
  l_exp <- T #criteria:
    if(max(diff(sub_mrg$d_xy_travel))>500){l_exp <- F} #breaks sanity check
    if(range(sub_mrg$d_xy_travel)[2]<14999){l_exp <- F} #sufficient length
    if((max(sub_mrg$HGT_RADR)-min(sub_mrg$HGT_RADR))>400){l_exp <- F} #not a profile (differences might be big when flying at constant pressure not radar height; need double checking manually)
  
  #only export leg if criteria are met
  if(l_exp){
    
    #create data.frame for export
    eddy.data <- data.frame(date=sub_mrg$date,
                            t_utc=yorkFLUX::UTC_Time(sub_mrg$date), #yorkFLUX might not be up to date
                            t_doy_utc=yorkFLUX::DOY_Time(sub_mrg$date),
                            t_doy_local=yorkFLUX::DOY_Time(sub_mrg$date),
                            lat_a=sub_mrg$LAT_GIN,
                            lon_a=sub_mrg$LON_GIN,
                            d_x_utm=sub_mrg$d_x_utm,
                            d_y_utm=sub_mrg$d_y_utm,
                            d_z_m=sub_mrg$HGT_RADR,
                            d_z_ABL=rep(bl,length(sub_mrg$date)), 
                            d_xy_travel=sub_mrg$d_xy_travel,
                            d_xy_flow=sub_mrg$d_xy_flow,
                            PSI_aircraft=sub_mrg$HDG_GIN,
                            uvw_aircraft=sub_mrg$TAS,
                            u_met=sub_mrg$U,
                            v_met=sub_mrg$V,
                            uv_met=sqrt(sub_mrg$U^2 + sub_mrg$V^2),
                            w_met=sub_mrg$W,
                            T_air=sub_mrg$temp,
                            p_air=sub_mrg$press*100, #hPa to Pa
                            FD_mole_H2O=sub_mrg$H2O_ppm* 1e-6, 
                            FD_mole_CH4=sub_mrg$CH4_ppb* 1e-9) 
    
    saveRDS(eddy.data,
            file=paste0("./processed/",fn,"_leg_",bk,".rds"))
    
    #tidy up
    rm(eddy.data)
    
  }
  
  #tidy up
  rm(l_exp,sub_mrg)
  
}



################################################################################
### QA/QC of the snipped runs ###

#check the altitude
ggplot(data=dm, aes(date, d_z_m, colour=d_z_m)) + geom_point()+ scale_color_viridis()

#if need to trim
dm <-  dm %>% filter(date<ymd_hms("2019-02-02 10:28:20"))

#check for distance travelled
range(dm$d_xy_travel)

#save
saveRDS(dm,"./processed/c137_leg_27.rds")



################################################################################
      ### Now feed the legs to the AEC algorithm and hope for the best ###
################################################################################




