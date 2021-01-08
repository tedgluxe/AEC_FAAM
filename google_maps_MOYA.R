################################################################################
### Google Maps plotting for MOYA AEC ###

#Contributions from: Will Drysdale, Dominika Pasternak.

################################################################################
### Loading packages ###

library(ggmap)
library(maps)
library(mapdata)


################################################################################
### Plotting ###

#Google map key, this one belongs to Jacob Shaw and can bes used only for MOYA paper!
ggmap::register_google(key = "AIzaSyAQxBVJMLIDDbIDkZ90MhgRDSwMSAOi3n8")

#adjust centre of the map from a standard stamenmap plot, zoom defines size of the map, scale - don't touch
map <- get_googlemap(center = c(lon = 27.1, lat = -15.9), zoom = 9, scale=2, maptype = "satellite")

#plot (ignore warings)
ggmap(map) + 
#  geom_rect(aes(xmin=27.61, xmax=27.95, ymin=(-14.6), ymax=(-14.24)),
#            fill=NA, color="lightblue", size=2) + #rectangle for C138
  geom_point(data = flux, 
             aes(lon,lat, colour=flux, size = flux)) +
  scale_color_viridis(option="inferno", name=bquote(''~CH[4]~ (mg~m^-2~h^-1)*'')) +
  xlim(c(min(flux$lon-0.01,na.rm = T),max(flux$lon+0.01,na.rm = T)))+ #forces reasonable size of the map
  ylim(c(min(flux$lat-0.005,na.rm = T),max(flux$lat+0.005,na.rm = T)))+
  theme(legend.text = element_text(size=14, colour="black"),
        axis.text = element_text(size=14, colour="black"),
        legend.title=element_text(size=14),
        axis.title = element_blank()
        ) +
  guides(size = FALSE, 
         color = guide_colourbar(barwidth = 1.5, 
                                 barheight = 24, 
                                 frame.colour = "black",
                                 ticks = F))

  

  