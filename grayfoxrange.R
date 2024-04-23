library(sf)

urocyon <- st_read("~/Downloads/redlist_species_data_40439b7f-ffb3-4ee3-8e20-b582b64ef368/data_0.shp")
world <- st_read("~/Downloads/World_Continents-7/World_Continents.shp")

if (length(urocyon$geometry[[1]]) >= 2) {
  # Remove the second list (polygon)
  # Keep the first and any after the second
  new_polygons <- urocyon$geometry[[1]][-c(2, 3)]  # -c(2, 3) removes the second and third elements
  
  # Update the geometry with the new list of polygons
  urocyon$geometry[[1]] <- st_multipolygon(new_polygons)
} else {
  cat("Not enough polygons to remove the second one.")
}




ggplot() + 
  geom_sf(data = world, size = 1, color = NA, fill = "gray") + 
  geom_sf(data = urocyon, size = 1, color = NA, fill = "goldenrod") + 
  geom_segment(aes(x = -100, y = 16, xend = -100, yend = 50))+
  coord_sf(xlim=c(-130,-60), ylim=c(0,51), expand = FALSE)+
  theme_classic()+
  theme(axis.text = element_text(color="black",size=9),
        legend.title = element_blank(),
        legend.position = "top",
        axis.title = element_blank())+
  scale_y_continuous(breaks=seq(0,50,10)) +
  scale_x_continuous(breaks=seq(-125,-60,15))

ggsave("~/Downloads/urocyon_gp.png", width = 5, height = 5, units = "in")
