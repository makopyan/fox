# Load required libraries
install.packages("sf")
install.packages("dplyr")
library(sf)
library(dplyr)

# Load the shapefile

data_points <- pop %>% select(long,lat)


data_points_sf <- st_as_sf(data_points, coords = c("long", "lat"), crs = 4326)

# Approximate buffer in degrees (very rough approximation)
# Example: 100 km is about 0.9 degrees at the equator
buffer_degree <- st_buffer(data_points_sf, dist = (100 / 111) + 1)

# Create convex hull around the combined buffer zones
combined_buffer <- st_union(buffer_degree)
convex_hull <- st_convex_hull(combined_buffer)

# Intersect convex hull with the species range, ensuring species_range is also in CRS 4326
species_range <- st_read("~/Downloads/redlist_species_data_40439b7f-ffb3-4ee3-8e20-b582b64ef368/data_0.shp")
species_range <- st_transform(species_range, crs = 4326)
convex_hull_within_range <- st_intersection(convex_hull, species_range)

ggplot() + 
  geom_sf(data = world, color = NA, fill = "sandybrown") + 
  geom_sf(data = species_range_utm, color = NA, fill = "royalblue") + 
  geom_sf(data = convex_hull_within_range, color = NA, fill = "brown3") +
  geom_sf(data = data_points_sf, size = 0.8, shape = 20) +
  coord_sf(xlim=c(-130,-62), ylim=c(2,51), expand = FALSE)+
  theme_classic()+
  theme(axis.text = element_text(color="black",size=9),
        legend.title = element_blank(),
        legend.position = "top",
        axis.title = element_blank())+
  scale_y_continuous(breaks=seq(0,50,10)) +
  scale_x_continuous(breaks=seq(-125,-60,15))

ggsave("~/Downloads/urocyon_outer_bounds.png", width = 5, height = 5, units = "in")



