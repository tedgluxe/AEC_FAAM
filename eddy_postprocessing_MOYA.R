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

###############################
### load and format files ### 

#load flight as 'dm' and its legs intervals as 'legs' 

#extract the eddy data
flux <-  data.frame(dm$INST$WAVE$date, dm$INST$WAVE$lat_a, dm$INST$WAVE$lon_a, dm$INST$WAVE$F_CH4_mass, dm$INST$error$ran.flux.F_CH4_mass, dm$INST$WAVE$u_star) %>% 
  dplyr::rename(date = dm.INST.WAVE.date,
         lat=dm.INST.WAVE.lat_a,
         lon=dm.INST.WAVE.lon_a,
         flux=dm.INST.WAVE.F_CH4_mass,
         error=dm.INST.error.ran.flux.F_CH4_mass,
         u_star=dm.INST.WAVE.u_star)
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
#C137
rm(leg_12,leg_14,leg_16,leg_17,leg_2,leg_26,leg_28,leg_3,leg_32,leg_4,leg_8)

leg_18$lod <- 4.61
leg_19$lod <- 2.74
leg_20$lod <- 2.50
leg_21$lod <- 1.56
leg_24$lod <- 1.57
leg_25$lod <- 2.54
leg_27$lod <- 4.86
leg_29$lod <- 2.38
leg_30$lod <- 3.85
leg_31$lod <- 1.64
leg_5$lod  <- 4.98

#put back together
legs <- mget(ls(pattern = "leg_.*"))  #list the legs
flux <- bind_rows(legs, .id = "column_label") #bind all the legs
flux <- filter(flux, u_star >= 0.15 ) #filter by friction velocity
rm(legs) #tidy up

#calculate error
flux$error_abs <- flux$error*flux$flux/100



################################################################################
### statistics and visualisation ###

#average & standard deviation
d <- flux
mean(d$flux)
mean(d$error_abs)/sqrt(nrow(d))
sd(d$flux)
rm(d)


#basic map
data_map <-  flux %>% na.omit() #pick data 
bbox = c(min(data_map$lon-0.2),min(data_map$lat-0.1),max(data_map$lon+0.2),max(data_map$lat+0.1)) #pick area
mymap = ggmap::get_stamenmap(bbox, zoom = 7) #pick zoom

ggmap(mymap)+
  geom_point(data = data_map, 
             aes(x = lon,
                 y = lat, 
                 colour = flux, 
                 size = flux)) +
  scale_color_viridis(option="magma") +
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
feet <-  dm$FOOT
feet <- feet %>% purrr::list_modify("c137_leg_17.rds"=NULL, "c137_leg_21.rds"=NULL, "c137_leg_26.rds"=NULL, "c137_leg_28.rds"=NULL, "c137_leg_3.rds"=NULL, "c137_leg_4.rds"=NULL)

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
df <- as.data.frame(all_s, xy=TRUE)

#plot footprint
ggplot(df)+
  geom_raster(aes(x,y, fill=layer)) +
  theme_bw() +
  scale_fill_viridis(option="inferno") +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=14), 
        legend.title = element_blank(), 
        axis.title = element_blank()) +
  guides(size = FALSE) 



################################################################################
###Other LODs###

#c128
flux$lod <- 1.17

#C138
leg_8$lod <- 0.87
leg_12$lod <- 0.96
leg_14$lod <- 2.02
leg_17$lod <- 0.95

fluxJ <- filter(flux, lat > (-14.6) & lat < (-14.24) &lon > 27.61 & lon < 27.95)







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
ra = readRDS("sad_raster.RDS")
happy_raster = raster::projectRaster(ra,crs = CRS("+proj=longlat +datum=WGS84 +no_defs"))
gis = happy_raster %>%
  gissr::ra_forify() %>%
  dplyr::rename(longitude = x,
                latitude = y)
sf = gis %>%
  gissr::sp_from_data_frame() %>%
  spTransform(CRS("+proj=utm +zone=35")) %>%
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
  scale_fill_viridis_c(option = "A")


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


