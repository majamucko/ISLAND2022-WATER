library(ggplot2)
install.packages("marmap")
install.packages("terra", type = "binary")
library(terra)
install.packages("ggOceanMaps")
library(marmap)
library(patchwork)
library(ggOceanMaps)
library(ggrepel)
library(RColorBrewer)
library(scales)


par(mar=c(3,4,2,2))
display.brewer.all()



# Define bounding boxes for the maps
south_adriatic_bbox <- data.frame(lon = c(16.5, 18.5), lat = c(41.5, 43.5))
lastovo_bbox <- data.frame(lon = c(16.5, 17.2), lat = c(42.6, 42.8))

# Fetch bathymetry data (higher resolution)
bathy_south_adriatic <- getNOAA.bathy(lon1 = 16.5, lon2 = 18.5, lat1 = 41.5, lat2 = 43.5, resolution = 0.05)
bathy_lastovo_korcula <- getNOAA.bathy(lon1 = 16.6, lon2 = 17.1, lat1 = 42.65, lat2 = 42.85, resolution = 0.04)

# Convert bathymetry data to a dataframe
bat_south <- as.data.frame(as.xyz(bathy_south_adriatic))
bat_lastovo <- as.data.frame(as.xyz(bathy_lastovo_korcula))

# Rename columns to match ggplot
colnames(bat_south) <- c("lon", "lat", "depth")
colnames(bat_lastovo) <- c("lon", "lat", "depth")

# Define study locations
locations <- data.frame(
  lon = c(16.88753),
  lat = c(42.72289),
  site = c("S1")
)

# Reverse the RdYlBu palette for ocean-blue to land-red gradient
min(bathy_south_adriatic)
max(bathy_south_adriatic)

color_palette <- scale_fill_gradientn(
  colors = rev(brewer.pal(9, "YlGnBu")),  # Reverse the RdYlBu palette
  values = rescale(c(-1300, -1000,  -800,   -600, - 300,  -100,  0)),  # Rescale values to match depth range
  name = "(m)",
  limits = c(-1300, 0)
)

plot_bathymetry <- function(data, title, sites = NULL) {
  p <- ggplot() +
    # Background layer: all tiles in gray (land and non-ocean areas)
    geom_tile(data = data, aes(x = lon, y = lat), fill = "gray", show.legend = FALSE) +
    # Overlay: only ocean areas (depth < 0) with gradient fill
    geom_tile(data = subset(data, depth < 0), aes(x = lon, y = lat, fill = depth)) +
    # Contour lines (applied to full data)
    geom_contour(data = data, aes(x = lon, y = lat, z = depth), 
                 color = "gray60", bins = 20) +
    # Apply the color palette (assumed to be defined externally)
    color_palette +
    xlab("Longitude") +
    ylab("Latitude") +
    theme_minimal(base_size = 20) +
    coord_fixed() +
    theme(
      legend.position = "bottom",         # Place the legend at the bottom
      legend.direction = "horizontal",      # Display the legend items horizontally
      legend.box = "horizontal",            # Arrange legend items in a horizontal box
      legend.key.width = unit(2, "cm"),       # Adjust the width of the legend keys
      legend.box.spacing = unit(0.5, "cm")    # Add spacing between the legend and plot
    )
  
  # Optionally add site locations
  if (!is.null(sites)) {
    p <- p + 
      geom_point(data = sites, aes(x = lon, y = lat), color = "black", size = 3) + 
      geom_text_repel(data = sites, aes(x = lon, y = lat, label = site), color = "black", size = 5)
  }
  
  return(p)
}


# Generate the maps
south_adriatic_map <- plot_bathymetry(bat_south)
south_adriatic_map


min(bathy_lastovo_korcula)
max(bathy_lastovo_korcula)

#Set color palette
color_palette2 <- scale_fill_gradientn(
  colors = rev(brewer.pal(11, "YlGnBu")),
  values = rescale(c(-200, -100, -50, 0)),  # Adjust these values based on your data range
  name = "(m)",
  limits = c(-200, 0)  # Adjust to reflect only ocean depths
)


plot_bathymetry <- function(data, title, sites = NULL) {
  p <- ggplot() +
    # First layer: fill all areas with gray (non-ocean areas and background)
    geom_tile(data = data, aes(x = lon, y = lat),
              fill = "gray", show.legend = FALSE) +
    # Second layer: overlay ocean data (depth < 0) with gradient fill
    geom_tile(data = subset(data, depth < 0), aes(x = lon, y = lat, fill = depth)) +
    # Optional contours: you might choose to keep or remove these
    geom_contour(data = data, aes(x = lon, y = lat, z = depth),
                 color = "black", bins = 2) +
    color_palette2 +
    xlab("Longitude") +
    ylab("Latitude") +
    coord_fixed() +
    theme_minimal(base_size = 20) +
    theme(
      legend.position = "bottom",        # Legend at the bottom
      legend.direction = "horizontal",     # Legend is horizontal
      legend.box = "horizontal",           # Lay out legend items horizontally
      legend.key.width = unit(2, "cm"),      # Adjust key width
      legend.box.spacing = unit(0.5, "cm")   # Spacing between legend and plot
    )
  
  if (!is.null(sites)) {
    p <- p + 
      geom_point(data = sites, aes(x = lon, y = lat), color = "black", size = 3) +
      geom_text_repel(data = sites, aes(x = lon, y = lat, label = site),
                      color = "black", size = 5)
  }
  
  return(p)
}

# Generate the map
lastovo_korcula_map <- plot_bathymetry(bat_lastovo, sites = locations)
print(lastovo_korcula_map)


#save 
ggsave("south_adriatic.pdf", south_adriatic_map, width = 10, height = 7, units = "in")
ggsave("lastovo_2.pdf", lastovo_korcula_map, width = 10, height = 7, units = "in")

library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(scales)
library(ggspatial)

# Custom PuBu palette (reversed for deep blue ocean to lighter shallows)
# Get the full PuBu palette and trim to start from 3rd color
pubu_colors <- rev(brewer.pal(9, "PuBu")[3:9])

blue_gradient <- scale_fill_gradientn(
  colors = pubu_colors,
  values = rescale(c(-1300, -1000, -800, -600, -300, -100, 0)),
  name = "Depth (m)",
  limits = c(-1300, 0)
)


# Clean base theme
theme_map_clean <- theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = unit(2, "cm"),
    legend.box.spacing = unit(0.5, "cm"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# Final bathymetry plotting function
plot_bathymetry <- function(data, title = NULL, sites = NULL) {
  ggplot() +
    # Land/background layer
    geom_tile(data = data, aes(x = lon, y = lat), fill = "gray90", show.legend = FALSE) +
    
    # Ocean tiles with gradient fill
    geom_tile(data = subset(data, depth < 0), aes(x = lon, y = lat, fill = depth)) +
    
    # Contour lines
    geom_contour(data = data, aes(x = lon, y = lat, z = depth),
                 color = "black", bins = 5, size = 0.3) +
    
    # Color palette
    blue_gradient +
    
    # Study sites
    {
      if (!is.null(sites)) {
        list(
          geom_point(data = sites, aes(x = lon, y = lat), color = "black", size = 3),
          geom_text_repel(data = sites, aes(x = lon, y = lat, label = site),
                          color = "black", size = 5)
        )
      }
    } +
    
    # Theme and layout
    xlab("Longitude") +
    ylab("Latitude") +
    coord_fixed() +
    theme_map_clean +
    
    # Scale and north arrow
    annotation_scale(location = "bl", width_hint = 0.3) +
    annotation_north_arrow(location = "bl", which_north = "true",
                           style = north_arrow_fancy_orienteering) +
    
    # Optional title
    {
      if (!is.null(title)) ggtitle(title) else NULL
    }
}

south_adriatic_map <- plot_bathymetry(bat_south, title = "South Adriatic", sites = locations)
lastovo_korcula_map <- plot_bathymetry(bat_lastovo, title = "Lastovo-Korčula", sites = locations)

print(south_adriatic_map)
print(lastovo_korcula_map)
#save 
ggsave("south_adriatic_updated.pdf", south_adriatic_map, width = 10, height = 7, units = "in")
ggsave("lastovo_2_updated.pdf", lastovo_korcula_map, width = 10, height = 7, units = "in")
