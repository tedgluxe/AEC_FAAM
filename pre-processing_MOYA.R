### Airborne Eddy Covariance files preparation from FAAM aircraft data ###

#Contributions from: Adam Vaughan, Freya Squires, Will Drysdale, Dominika Pasternak.


#load packages
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


##############################################################
#        CHANGE FLIGHT NUMBER, ZONE & BOUNDARY layer         #
##############################################################
# Uganda all 35, Zambia 36, Finland 35



################################################################################
#set working directory
setwd("G:/Shared drives/eddy4R AEC/pre_processing")

#set sensitivity for change in altitude
alt_change <- 1.5 #+- mean*1.5 times standard deviation 

#files
#AIMMS_20_files <- list.files("./aimms_20hz",pattern = ".nc") 
core_32_files <-  list.files("./aimms_32hz",pattern = ".nc")
AIMMS_01_files <- list.files("./aimms_1hz",pattern = ".nc")
FGGA_10_files <- list.files("./fgga_10hz",pattern = ".na")
FGGA_h2o_files <- list.files("./fgga_h2o",pattern = ".txt")
thermistor_files <- list.files("./thermistor_16Hz",pattern = ".csv")

#flight number
FLIGHT_num <- stringr::str_sub(FGGA_10_files, start= -7) %>% gsub(".na","",.)
fn <- FLIGHT_num[6] 

#start date
origin_all <-  stringr::str_sub(FGGA_10_files, start= -19, end=-12) %>% gsub(".na","",.)
origin <-  origin_all[6]
origin <-  paste0(origin, " 00:00")

#frequencies
#aimms_freq <- 20
core_freq <-  32
fgga_freq <- 10
rad_freq <-  1
temp_freq <- 32

################################################################################

# AIMMS (Jacobs code adopted) 
#https://www.r-bloggers.com/a-netcdf-4-in-r-cheatsheet/

# #open
# aimms_ncdf <- ncdf4::nc_open(paste0("./aimms_20hz/",AIMMS_20_files[grep(fn,AIMMS_20_files,ignore.case=TRUE)]))
# 
# #get variables
# #temp <- ncvar_get(aimms_ncdf, attributes(aimms_ncdf$var)$names[1]) %>% as.vector()
# U <- ncvar_get(aimms_ncdf, attributes(aimms_ncdf$var)$names[3]) %>% as.vector()
# V <- ncvar_get(aimms_ncdf, attributes(aimms_ncdf$var)$names[4]) %>% as.vector()
# W <- ncvar_get(aimms_ncdf, attributes(aimms_ncdf$var)$names[5]) %>% as.vector()
# lat <- ncvar_get(aimms_ncdf, attributes(aimms_ncdf$var)$names[6]) %>% as.vector()
# long <- ncvar_get(aimms_ncdf, attributes(aimms_ncdf$var)$names[7]) %>% as.vector()
# altitude <- ncvar_get(aimms_ncdf, attributes(aimms_ncdf$var)$names[8]) %>% as.vector()
# heading <- ncvar_get(aimms_ncdf, attributes(aimms_ncdf$var)$names[14]) %>% as.vector()
# tas <- ncvar_get(aimms_ncdf, attributes(aimms_ncdf$var)$names[15]) %>% as.vector()
# ws <- ncvar_get(aimms_ncdf, attributes(aimms_ncdf$var)$names[18]) %>% as.vector()
# wd <- ncvar_get(aimms_ncdf, attributes(aimms_ncdf$var)$names[19]) %>% as.vector()
# 
# #get time
# aimms_time <- ncvar_get(aimms_ncdf, attributes(aimms_ncdf$dim)$names[1]) %>% as.vector()
# date <- strptime(x = origin, format ="%Y%m%d %H:%M") + (aimms_time)
# date <- base::as.POSIXct(seq.POSIXt(from = min(date)+(1/aimms_freq),
#                                     to = max(date)+1, #20 points / 1s less than aimms has, who knows why
#                                     by = 1/aimms_freq,
#                                     tz = "UTC"))  
# 
# #put all together
# aimms <- data.frame(date, lat, long, altitude,  U,V,W, heading,tas, ws, wd) %>% na.omit
# 
# #tidy up
# rm(wd, ws, tas, heading, U, V, W, lat, long, altitude)



################################################################################

#CORE 32Hz

#open
core32 <- ncdf4::nc_open(paste0("./aimms_32hz/",
                                core_32_files[grep(fn,core_32_files,ignore.case=TRUE)]))

#create a list of variables that we want to select 
eddyvars <- c("ROLL_GIN","PTCH_GIN", "U_C", "V_C", "W_C", "LAT_GIN", "LON_GIN","HDG_GIN","TAS", "PS_RVSM")

#turn NetCDF into normal data frame
for (i in 1:length(eddyvars)) {
  vname <- eddyvars[i]
  raw <- as.vector(ncdf4::ncvar_get(core32,vname,collapse_degen=FALSE))
  if(i==1){
    AIMMS_32hz <- data.frame(raw)
    names(AIMMS_32hz) <- vname
  } else {
    AIMMS_32hz <- cbind(AIMMS_32hz,raw)
    names(AIMMS_32hz)[ncol(AIMMS_32hz)] <- vname
  }
}
rm(i,raw,vname)

#get time
core_time <- ncvar_get(core32, attributes(core32$dim)$names[1]) %>% as.vector()
date <- strptime(x = origin, format ="%Y%m%d %H:%M") + (core_time)
AIMMS_32hz$date <- base::as.POSIXct(seq.POSIXt(from = min(date)+(1/core_freq),
                                               to = max(date)+1, 
                                               by = 1/core_freq,
                                               tz = "UTC"))




# #get time
# Time <- as.vector(ncdf4::ncvar_get(core32,"Time"))
# start <- as.POSIXct(strptime(substr(core32$var[1][[names(core32$var)[1]]]$dim[[1]]$units,15,33),
#                              "%Y-%m-%d %H:%M:%S"),tz="UTC") + Time[1]
# end <- as.POSIXct(strptime(substr(AIMMS_32hz$var[1][[names(AIMMS_32hz$var)[1]]]$dim[[1]]$units,15,33),
#                            "%Y-%m-%d %H:%M:%S"),tz="UTC") + Time[1] + (nrow(AIMMS_32hz)/core_freq)
# AIMMS_32hz$date <- base::as.POSIXct(seq.POSIXt(from = start+(1/core_freq),
#                                          to = end,
#                                          by = 1/core_freq,
#                                          tz = "UTC"))

#tidy up
rm( core32, date)


################################################################################

#CORE 1Hz (based on the original code from now on)

#open
tmp_AIMMS <- ncdf4::nc_open(paste0("./aimms_1hz/",
                                   AIMMS_01_files[grep(fn,AIMMS_01_files,ignore.case=TRUE)]))

#get variables
HGT_RADR <- as.vector(ncdf4::ncvar_get(tmp_AIMMS,"HGT_RADR",collapse_degen=FALSE))
press <- as.vector(ncdf4::ncvar_get(tmp_AIMMS,"PS_RVSM",collapse_degen=FALSE))
ROLL_GIN <- as.vector(ncdf4::ncvar_get(tmp_AIMMS,"ROLL_GIN",collapse_degen=FALSE))
co <- as.vector(ncdf4::ncvar_get(tmp_AIMMS,"CO_AERO",collapse_degen=FALSE))
co_flag <- as.vector(ncdf4::ncvar_get(tmp_AIMMS,"CO_AERO_FLAG",collapse_degen=FALSE))
core <- data.frame(HGT_RADR, press, ROLL_GIN, co, co_flag)

#get time
Time <- as.vector(ncdf4::ncvar_get(tmp_AIMMS,"Time"))
start <- as.POSIXct(strptime(substr(tmp_AIMMS$var[1][[names(tmp_AIMMS$var)[1]]]$dim[[1]]$units,15,33),
                             "%Y-%m-%d %H:%M:%S"),tz="UTC") + Time[1]
end <- as.POSIXct(strptime(substr(tmp_AIMMS$var[1][[names(tmp_AIMMS$var)[1]]]$dim[[1]]$units,15,33),
                           "%Y-%m-%d %H:%M:%S"),tz="UTC") + Time[1] + (nrow(core)/rad_freq)
core$date <- base::as.POSIXct(seq.POSIXt(from = start+(1/rad_freq),
                                         to = end,
                                         by = 1/rad_freq,
                                         tz = "UTC"))

#tidy up
rm(Time,start,end, tmp_AIMMS, HGT_RADR, press, ROLL_GIN)


################################################################################
#temperature

#get data 
tmp_temp <- read.csv(file=paste0("./thermistor_16Hz/",
                                 thermistor_files[grep(fn,thermistor_files,ignore.case=TRUE)]),
                     header = T, sep = ",",
                     stringsAsFactors = F) %>%
  dplyr::rename(date=datetime, temp=TAT_THERM_U) %>%
  dplyr::mutate(date=as.POSIXct(strptime(date, "%Y-%m-%d %H:%M:%OS", tz="UTC"))) %>%
  dplyr::select(date,temp)

################################################################################

#FGGA 10 Hz

#read
FGGA_10hz <- read.delim(file=paste0("./fgga_10hz/",
                                    FGGA_10_files[grep(fn,FGGA_10_files,ignore.case=TRUE)]),
                        skip = 62, header = F, sep=" ") %>%
  dplyr::rename(.,date=V1,CO2_ppm=V2,CH4_ppb=V4) %>%
  dplyr::select(.,c(date,CO2_ppm,CH4_ppb))

#get date in right format
FGGA_10hz$date <- FGGA_10hz$date + strptime(x = origin, format ="%Y%m%d %H:%M")

#filter
FGGA_10hz$CH4_ppb[FGGA_10hz$CH4_ppb < 0] <- NA 
FGGA_10hz$CO2_ppm[FGGA_10hz$CO2_ppm < 0] <- NA
FGGA_10hz$CH4_ppb[FGGA_10hz$CH4_ppb >9500] <- NA 
FGGA_10hz$CO2_ppm[FGGA_10hz$CO2_ppm >950] <- NA




################################################################################

# FGGA H2O

#get data (ignore warning)
H2O_10Hz <- read.delim(file=paste0("./fgga_h2o/",
                                   FGGA_h2o_files[grep(fn,FGGA_h2o_files,ignore.case=TRUE)]),
                       header = T, sep = ",",
                       stringsAsFactors = F) %>%
  dplyr::rename(date=Time) %>%
  dplyr::mutate(date=as.POSIXct(strptime(date, "%d/%m/%Y %H:%M:%OS", tz="UTC"))) %>%
  dplyr::select(date,H2O_ppm)

#filter
H2O_10Hz$H2O_ppm[H2O_10Hz$H2O_ppm < 1] <- NA

FL_mrg$co[FL_mrg$co_flag>0] <- NA


################################################################################

# Crop all to FGGA 10Hz

core <- core %>% subset(.,date>=min(FGGA_10hz$date) & date<=max(FGGA_10hz$date))
AIMMS_32hz <- AIMMS_32hz %>% subset(.,date>=min(FGGA_10hz$date) & date<=max(FGGA_10hz$date))
H2O_10Hz <- H2O_10Hz %>% subset(.,date>=min(FGGA_10hz$date) & date<=max(FGGA_10hz$date))
tmp_temp <-  tmp_temp %>% subset(.,date>=min(FGGA_10hz$date) & date<=max(FGGA_10hz$date))


################################################################################

#interpolate up

tmp_date <- data.frame(date=AIMMS_32hz$date)
FL_mrg <- cbind(tmp_date,
                sapply(names(core)[-which(names(core)=="date")],
                       function(x)
                         x <- approx(core$date,core[[x]],
                                     tmp_date$date)$y))%>%
  cbind(.,
        sapply(names(tmp_temp)[-which(names(tmp_temp)=="date")],
               function(x)
                 x <- approx(tmp_temp$date,tmp_temp[[x]],
                             tmp_date$date)$y))%>%
  cbind(.,
        sapply(names(FGGA_10hz)[-which(names(FGGA_10hz)=="date")],
               function(x)
                 x <- approx(FGGA_10hz$date,FGGA_10hz[[x]],
                             tmp_date$date)$y)) %>%
  cbind(.,
        sapply(names(H2O_10Hz)[-which(names(H2O_10Hz)=="date")],
               function(x)
                 x <- approx(H2O_10Hz$date,H2O_10Hz[[x]],
                             tmp_date$date)$y)) %>%
  cbind(.,AIMMS_32hz[-which(names(H2O_10Hz)=="date")]) 

FL_mrg <- FL_mrg[,1:length(FL_mrg)-1]
FL_mrg <- na.omit(FL_mrg)
#tidy up
rm(AIMMS_32hz,core,FGGA_10hz,H2O_10Hz,tmp_date, tmp_temp)

#cut out fires
FL_mrg$CH4_ppb[FL_mrg$date > ymd_hms("2019-02-02 10:10:21") & FL_mrg$date < ymd_hms("2019-02-02 10:10:51")] <-  NA

FL_mrg <- na.omit(FL_mrg)




################################################################################

#snip legs

#C136
#FL_mrg <- FL_mrg %>%
#dplyr::filter(LON_GIN>20)%>%
#dplyr::filter(LAT_GIN>(-15)) %>%
#na.omit()


#filter data based on altitude and roll angle
FL_brk <- FL_mrg %>%
  dplyr::filter(HGT_RADR<610) %>%
  dplyr::filter(abs(ROLL_GIN)<20) %>%
  na.omit()

#map to check if sensible
# bbox = c(min(FL_mrg$LON_GIN-0.2),min(FL_mrg$LAT_GIN-0.1),max(FL_mrg$LON_GIN+0.2),max(FL_mrg$LAT_GIN+0.1))
# mymap = ggmap::get_stamenmap(bbox, zoom = 5)
# ggmap(mymap)+
#  geom_point(data = FL_brk,
#             aes(LON_GIN,LAT_GIN, colour=HGT_RADR),
#             size = 2) +
#  scale_color_viridis(option="magma") +
#  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=14), legend.title = element_blank())


#find heading breaks that are longer than 10 degrees
t_dif <- diff(FL_brk$date) %>% as.numeric(.)
t_gaps <- which(abs(t_dif)>5) #time difference of 5 seconds / 500 m

#Calculate start and end positions for each leg
t_start <- FL_brk$date[c(1,(t_gaps[1:(length(t_gaps)-1)]+1))]
t_end <- FL_brk$date[t_gaps]
t_leg <-  t_end - t_start 



################################################################################

#loop for each leg

for(bk in 1:length(t_start)){
  #bk <-  2
  
  #clip data to start and end times for each leg
  sub_mrg <- FL_mrg %>% subset(date>=t_start[bk] & date <= t_end[bk])
  
  #Calculate UTM coordinates
  d <- data.frame(lon=sub_mrg$LON_GIN, lat=sub_mrg$LAT_GIN)
  coordinates(d) <- c("lon", "lat")
  proj4string(d) <- CRS("+proj=longlat")
  
  ### CHANGE ZONE TO THE RIGHT ONE!!!!!!!!!!!! ###
  cov <- spTransform(d, CRS("+proj=utm +zone=35 ellps=WGS84"))
  
  #Calculate ground distance covered by the aircraft[m]
  stretch <- c(0,sqrt((diff(cov@coords[,1]))^2 + (diff(cov@coords[,2]))^2))
  stretch <- cumsum(stretch)
  parcel <- cumsum(sub_mrg$TAS * 1/2)
  parcel <- parcel - min(parcel, na.rm=TRUE)
  
  #Calculate the volume of air passed by the aircraft
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
  
  #final check for high altitudes and rolls
  sub_mrg <- sub_mrg %>%
    dplyr::filter(HGT_RADR<610.01) %>%
    dplyr::filter(abs(ROLL_GIN)<20.01) %>%
    na.omit()
  
  #distance traveled check
  sub_mrg$d_xy_travel <- sub_mrg$d_xy_travel - sub_mrg$d_xy_travel[1]
  sub_mrg$d_xy_flow <- sub_mrg$d_xy_flow - sub_mrg$d_xy_flow[1]
  
  l_exp <- T
  
  if(max(diff(sub_mrg$d_xy_travel))>500){l_exp <- F} #breaks sanity check
  if(range(sub_mrg$d_xy_travel)[2]<14999){l_exp <- F} #too short
  if((max(sub_mrg$HGT_RADR)-min(sub_mrg$HGT_RADR))>400){l_exp <- F} #profile
  
  #only export leg if conditions are made
  if(l_exp){
    
    #create data.frame for export
    eddy.data <- data.frame(date=sub_mrg$date,
                            t_utc=yorkFLUX::UTC_Time(sub_mrg$date), #yorkFLUX is outdated
                            t_doy_utc=yorkFLUX::DOY_Time(sub_mrg$date),
                            t_doy_local=yorkFLUX::DOY_Time(sub_mrg$date),
                            lat_a=sub_mrg$LAT_GIN,
                            lon_a=sub_mrg$LON_GIN,
                            d_x_utm=sub_mrg$d_x_utm,
                            d_y_utm=sub_mrg$d_y_utm,
                            d_z_m=sub_mrg$HGT_RADR,
                            d_z_ABL=rep(800,length(sub_mrg$date)), ##boundary layer (meters)
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
    
    #clean up
    rm(eddy.data)
    
  }
  
  #clean up
  rm(l_exp,sub_mrg)
  
}



##########################################################################

#check the runs
dm <- readRDS("G:/Shared drives/eddy4R AEC/pre_processing/processed/c137_leg_33.rds")
ggplot(data=dm, aes(date, d_z_m, colour=d_z_m)) + geom_point()+ scale_color_viridis()


# if need to trim
dm <-  dm %>% filter(date<ymd_hms("2019-02-02 10:28:20"))
dm$d_xy_travel <- dm$d_xy_travel - min(dm$d_xy_travel)
dm$d_xy_flow <- dm$d_xy_flow - min(dm$d_xy_flow)
max(diff(dm$d_xy_travel))
range(dm$d_z_m)
range(dm$d_xy_travel)

#save
saveRDS(dm,"./processed/c137_leg_27.rds")



FL_mrg2 <- FL_mrg %>% filter(date < ymd_hms("2019-01-25 12:30:00"))

ggplot()+geom_point(aes(FL_mrg2$temp, FL_mrg2$HGT_RADR), colour="blue") +theme_bw() 
#+labs(x="Time", y="Temperature / K") 


#############################

#crop data by date & time
dm <- FL_mrg %>%
  filter(between(date, 
                 ymd_hms("2019-01-25 12:04:00"),
                 ymd_hms("2019-01-25 14:23:00"))) 
dm$CH4_ppb[dm$CH4_ppb<1900] <- 1900
dm$CH4_ppb[dm$CH4_ppb>2400] <- 2400

#limits=c(1900, 2400)


bbox = c(min(dm$LON_GIN-0.1),min(dm$LAT_GIN-0.05),max(dm$LON_GIN+0.1),max(dm$LAT_GIN+0.05))
mymap = ggmap::get_stamenmap(bbox, zoom = 5)
ggmap(mymap)+
  geom_point(data = dm,
             aes(LON_GIN,LAT_GIN, colour=CH4_ppb),
             size = 2, alpha=.6) +
  scale_color_viridis(option="magma") +
  labs(x="Longitude", y="Latitude", title=bquote(''~CH[4]~ (ppb)*'')) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=14), legend.title = element_blank())















