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
library(rgdal)
library(rhdf5)
library(data.table)
library(openair)
library(Rmisc)
library(raster)
library(rgdal)
library(sp)
library(wsdmiscr)
library(sf)

###############################
### load and format files ### 

#load flight as 'dm' and its legs intervals as 'legs' 

#extract the eddy data
flux <-  data.frame(dm$INST$WAVE$date, dm$INST$WAVE$lat_a, dm$INST$WAVE$lon_a, dm$INST$WAVE$F_CO2_mass, dm$INST$error$ran.flux.F_CO2_mass, dm$INST$error$sys.flux.F_CO2_mass, dm$INST$itcs$ItcFlag_w_hor, dm$INST$itcs$Itc_w_hor, dm$INST$REYN$d_z_m, dm$INST$REYN$d_z_ABL) %>% 
  dplyr::rename(date = dm.INST.WAVE.date,
         lat=dm.INST.WAVE.lat_a,
         lon=dm.INST.WAVE.lon_a,
         flux=dm.INST.WAVE.F_CO2_mass,
         error_ran=dm.INST.error.ran.flux.F_CO2_mass,
         error_sys=dm.INST.error.sys.flux.F_CO2_mass,
         w_flag=dm.INST.itcs.ItcFlag_w_hor,
         w_hor=dm.INST.itcs.Itc_w_hor,
         alt = dm.INST.REYN.d_z_m,
         BLH = dm.INST.REYN.d_z_ABL)
flux$date <-  as.POSIXct(flux$date)

#get the legs times
t_start <- legs$t_start
t_end <- legs$t_end
rm(legs)

#cut flux into legs to add LODs
for(b in 1:length(t_start)){
  sub_mrg <- flux %>% subset(date>=t_start[b] & date <= t_end[b])
  assign(paste0("leg_",b), sub_mrg)
}
rm(sub_mrg, b) #tidy up

#remove empty data frames
legs <-  ls(pattern = "leg_.*") #list the legs
empty <- legs[sapply(legs, function(x) nrow(get(x)) == 0L)] #pick empty ones
rm(list = empty, empty, legs) #remove empty ones and tidy up

#remove unwanted legs and add LODs
leg_14 <- leg_45
leg_8 <-  leg_50
rm(leg_21,leg_31,leg_5)



#put back together
legs <- mget(ls(pattern = "leg_.*"))  #list the legs
flux <- bind_rows(legs, .id = "column_label") #bind all the legs
rm(legs) #tidy up




################################################################################
### flux correction for boundary layer depth ###

#correct the flux
C132$flux_c <-  C132$flux/(1-(C132$alt/(0.8*C132$BLH)))

#calculate error
C132$error_abs <- (C132$error_ran+C132$error_sys)*abs(C132$flux_c)/100

#filtering
C132 <- C132 %>% filter(w_flag<1)

fluxJ <- filter(flux2, lat > (-14.6) & lat < (-14.24) &lon > 27.61 & lon < 27.95)

flux4 <- openair::timeAverage(flux2, avg.time="30 sec")


flux_C137 <- flux

write.csv(flux4, "G:/My Drive/eddy_new/MOYA2/plotting/C137_30s_avg_CO2.csv")


# plot the bad boi

#correction
ggplot(flux2)+
  geom_point(aes(x=alt/BLH,y=(flux_c-flux)*100/flux), size=2, shape=1) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=14), 
        legend.title = element_blank()) +
  #xlim(min(flux$flux_c), max(flux$flux_c))+
  labs(x="Relative boundary layer depth", y="Percantage correction", title="C137")

#error
flux3 <-  flux2 %>% filter(flux<1)

ggplot(flux3)+
  geom_ribbon(aes(x=seq(1,nrow(flux3)), ymin=flux_c-(flux_c*error_ran/100), ymax=flux_c+(flux_c*error_ran/100)), fill="skyblue")+
  geom_hline(yintercept = 0, colour="skyblue4", size=2)+
  geom_point(aes(x=seq(1,nrow(flux3)),y=flux_c), size=2, shape=1) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x="Index", y=bquote(''~EC~(CH[4]~mg~m^-2~h^-1)*''), title="C138")



#flagging
ggplot(flux)+
  geom_point(aes(x=w_hor,y=flux_c, fill=w_flag), size=2, shape=21) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=14), 
        legend.title = element_blank()) +
  labs(x="Integral Turbulence Characteristics (%)" , y=bquote(''~EC~(CH[4]~mg~m^-2~h^-1)*''), title="C137, threshold 100%")+
  guides(fill = FALSE) 


#fire
ggplot(fire3)+
  geom_point(aes(x=CPC_CNTS,y=flux_c, colour=CO_AERO), size=2) +
  theme_bw()+
  scale_color_viridis(option="viridis") +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=14)) +
  labs(x="CPC counts" , y=bquote(''~EC~(CO[2]~mg~m^-2~h^-1)*''), colour="CO (ppb)", title="C128, C129 and C132")

#value vs alt + corr
ggplot(flux2)+
  geom_point(aes(x=alt,y=flux_c, colour=(flux_c-flux)*100/flux), size=2) +
  theme_bw()+
  scale_color_viridis(option="viridis") +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=14)) +
  labs(x="Radar height in m", y=bquote(''~EC~(CH[4]~mg~m^-2~h^-1)*''), colour="% correction", title="C137")

################################################################################
### statistics and visualisation ###

#average & standard deviation
d <- C132
mean(d$flux_c)
mean(d$error_abs)/sqrt(nrow(d))
sd(d$flux_c)
rm(d)

write.csv(dg2, "G:/My Drive/eddy_new/BL_corrected_for_paper_v4/C137_v10_foot_pt2.csv")

#basic map
data_map <-  flux4 %>% na.omit() #pick data 
bbox = c(min(data_map$lon-0.2),min(data_map$lat-0.1),max(data_map$lon+0.2),max(data_map$lat+0.1)) #pick area
mymap = ggmap::get_stamenmap(bbox, zoom = 7) #pick zoom

ggmap(mymap)+
  geom_point(data = data_map, 
             aes(x = lon,
                 y = lat, 
                 colour = flux_c, 
                 size = flux_c)) +
  scale_color_viridis(option="inferno", #limits=c(-10,150)
                      ) +
  labs(title=bquote(''~CH[4]~ (mg~m^-2~h^-1)*''))+
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=14), 
        legend.title = element_blank(), 
        axis.title = element_blank()) +
  guides(size = FALSE) 


#flux + error plots
ggplot() +
  geom_errorbar(data = flux, 
                aes(x = seq(1,length(flux)),
                    ymin = flux-error*flux/100, 
                    ymax = flux+error*flux/100), 
                width = 0.1, 
                colour = "dodgerblue4")+
  geom_point(data = flux,
             aes(x = seq(1,length(flux)), 
                 y = flux), 
             color = "dodgerblue", 
             size = 3) +
  theme_bw() +
  labs(x="Index", y=bquote(''~CH[4]~ (mg~m^-2~h^-1)*'')) +
  theme(plot.title = element_text(hjust = 0.5),  
        text = element_text(size=18))+
  guides(size = FALSE) 



################################################################################
### Footprint ###

#remove unwanted legs
#C137
feet <-  list(dm$FOOT$c137_leg_19.rds, dm$FOOT$c137_leg_24.rds, dm$FOOT$c137_leg_25.rds, dm$FOOT$c137_leg_27.rds,dm$FOOT$c137_leg_29.rds)
#feet <- feet %>% purrr::list_modify("c138_leg_8.rds"=NULL, "c138_leg_17.rds"=NULL, "c138_leg_14.rds"=NULL, "c138_leg_12.rds"=NULL)

#find the common extend
xmin <- min(sapply(feet, function(x) x@extent@xmin))
ymin <- min(sapply(feet, function(x) x@extent@ymin))
xmax <- max(sapply(feet, function(x) x@extent@xmax))
ymax <- max(sapply(feet, function(x) x@extent@ymax))
newextent=c(xmin, xmax, ymin, ymax)

#assign the new extend
legs <- lapply(feet, extend, newextent, value=0)

#sum all the footprints
all <-  stack(legs) #not sure why but otherwise doesn't work
all_s <-  calc(all, fun=sum) #sum

#convert coordinates
all_s <- projectRaster(all_s,crs = CRS("+proj=longlat +datum=WGS84 +no_defs"))

#convert into data frame to plot
dg <- as.data.frame(all_s, xy=TRUE)

df$x <- df$x-6

dg2 <- dg %>% filter(layer>wsdmiscr::f_cont(dg$layer, 90))


#plot footprint
ggplot(dg2)+
  geom_raster(aes(x,y, fill=layer)) +
  theme_bw() +
  scale_fill_viridis(option="inferno") +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=14), 
        legend.title = element_blank(), 
        axis.title = element_blank()) +
  guides(size = FALSE) 



################################################################################

fire <- read.csv("G:/My Drive/eddy_new/MOYA2/plotting/C137_C138_cpc_co.csv")
fire$date <- as.POSIXct(fire$date)


fire2 <- cbind(flux_all,
                 sapply(names(fire)[-which(names(fire)=="date")],
                       function(x)
                         x <- approx(fire$date,fire[[x]],
                                     flux_all$date)$y)) 


fire_Z <- rbind(fire2, fire8)

saveRDS(fire2, "G:/My Drive/eddy_new/MOYA2/Uganda_CO2_fires.RDS")

################################################################################
######## PS Will's magic ######## 

# make happy raster
parse_inter = function(x,gis,temp){
  if(length(x) == 0)
    return(temp)
  else{
    out = gis[x,] %>%
      dplyr::summarise_all(mean,na.rm = T)
  }
  return(out)
}
ra = all_s
happy_raster = raster::projectRaster(ra,crs = CRS("+proj=longlat +datum=WGS84 +no_defs"))
gis = happy_raster %>%
  gissr::ra_forify() %>%
  dplyr::rename(longitude = x,
                latitude = y)
sf = gis %>%
  gissr::sp_from_data_frame() %>%
  spTransform(CRS("+proj=utm +zone=36")) %>%
  st_as_sf()
grid = st_make_grid(sf,
                    square = T,
                    cellsize = c(800, 800)) %>%
  st_sf()
intersects = st_intersects(grid,sf)
temp = gis[1,] %>%
  dplyr::mutate_all(~NA)
inter_dat = intersects %>%
  as.list() %>%
  purrr::map_df(.x = .,
                .f = parse_inter,
                gis = gis,
                temp = temp)
inter_dat %>%
  mutate(longitude = round(longitude,2),
         latitude = round(latitude,2)) %>%
  ggplot()+
  geom_raster(aes(longitude,latitude,fill = value))+
  scale_fill_viridis_c(option = "A")+
  geom_point(data=flux4,aes(x=lon,y=lat),colour="white")



gmap = ggmap::get_stamenmap(bbox = c(min(inter_dat$longitude),
                                     min(inter_dat$latitude),
                                     max(inter_dat$longitude),
                                     max(inter_dat$latitude)))
plotDat = inter_dat %>%
  mutate(longitude = round(longitude,2),
         latitude = round(latitude,2),
         value = ifelse(value < wsdmiscr::f_cont(inter_dat$value, 90),NA,value))
ggmap(gmap)+
  coord_cartesian()+
  geom_raster(data = plotDat,aes(longitude,latitude,fill = value),alpha = 0.8)+
  scale_fill_viridis_c(option = "A",na.value = "NA")


