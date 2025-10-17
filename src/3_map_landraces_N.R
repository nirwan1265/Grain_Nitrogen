####  Libraries
library(ggplot2)
library(maps)
library(ggspatial)
library(dplyr)
library(raster)
library(sp)
library(RColorBrewer) 
library(GGally)
library(scatterplot3d)
library(plotly)
library(fields)
library(viridis)
library(tidyverse)
library(mapdata)
library(patchwork)
library(vroom)
library(grid)     


# Load and filter data
dir_maize <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/data/Phenotype/"
maize_geo <- read.csv(paste0(dir_maize, "taxa_geoloc_pheno.csv")) %>%
  dplyr::filter(
    sp == "Zea mays",
    GEO3 %in% c("Caribbean", "Meso-America", "South America")
  ) %>%
  dplyr::select(2, 6, 7)
names(maize_geo) <- c("Genotypes", "long", "lat")

# Nitrate values
data_values <- vroom("/Users/nirwantandukar/Documents/Research/data/Phenotypes/maize_N_NHx.csv")
colnames(data_values)[1] <- "Genotypes"

# Combine
maize_geo <- maize_geo %>%
  left_join(data_values, by = "Genotypes") %>%
  dplyr::select(Genotypes, long, lat, TN_maize)



# Convert to spatial points
coordinates(maize_geo) <- ~long + lat
proj4string(maize_geo) <- CRS("+init=epsg:4326")

# Convert back to data frame for ggplot
maize_df <- as.data.frame(maize_geo)


# Log-transform NHx values
#maize_df$log_nhx <- log2(maize_df$NHx_maize)
maize_df$log_nhx <- maize_df$TN_maize
str(maize_df)


# helper: margin() if available, otherwise fall back to grid::unit
marg <- function(top = 0, right = 0, bottom = 0, left = 0, unit = "pt") {
  if (exists("margin", envir = asNamespace("ggplot2"))) {
    ggplot2::margin(top, right, bottom, left)
  } else {
    grid::unit(c(top, right, bottom, left), unit)
  }
}

## ─────────────────────────────────────────────────────────
##  A) MAIN MAP
## ─────────────────────────────────────────────────────────
main_map <- ggplot() +
  borders(
    "world",
    xlim   = c(-120, -30),
    ylim   = c(-50,  40),
    fill   = "gray95",
    colour = "gray70",
    size   = 0.2
  ) +
  geom_point(
    data  = maize_df,
    aes(long, lat, colour = log_nhx),
    size  = 0.8,
    alpha = 0.8
  ) +
  coord_fixed(ratio = 1.3, xlim = c(-120, -30), ylim = c(-50, 40)) +
  scale_colour_viridis_c(
    option = "D",  # "D" = viridis, "C" = magma, "B" = cividis
    name   = expression("Soil N (cg/kg)"),
    guide  = guide_colourbar(
      direction    = "vertical",
      barwidth     = unit(0.4, "cm"),
      barheight    = unit(4, "cm"),
      ticks.colour = "black",
      title.position = "top",
      title.hjust  = 0.5
    )
  ) +
  labs(
    x        = "Longitude",
    y        = "Latitude",
    title    = "Geographic distribution of soil nitrogen",
    subtitle = "Landrace collection sites across the Americas colored by soil nitrogen"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(size = 0.5, colour = "black"),
    axis.ticks       = element_line(size = 0.5, colour = "black"),
    legend.position  = c(0.9, 0.2),
    legend.justification = "left",
    legend.title     = element_text(size = 8, face = "bold"),
    legend.text      = element_text(size = 7),
    plot.title       = element_text(face = "bold", size = 11, hjust = 0),
    plot.subtitle    = element_text(size = 9, hjust = 0, margin = marg(bottom = 6)),
    plot.margin      = marg(6, 12, 6, 6)
  )

## ─────────────────────────────────────────────────────────
##  B) HISTOGRAM (larger inset, no grids)
## ─────────────────────────────────────────────────────────

# Calculate quantiles
quantiles <- quantile(maize_df$log_nhx, probs = c(0.25, 0.75), na.rm = TRUE)

hist_plot <- ggplot(maize_df, aes(log_nhx, fill = ..x..)) +
  geom_histogram(bins = 30, color = "black", size = 0.3) +
  # Add vertical lines for quantiles
  geom_vline(xintercept = quantiles, 
             color = "red", 
             linetype = "dashed", 
             size = 0.8,
             alpha = 0.7) +
  # Add text annotations
  annotate("text", 
           x = quantiles[1], 
           y = max(ggplot_build(hist_plot)$data[[1]]$count) * 0.95,
           label = "Bottom \n 25%",
           color = "red",
           hjust = 1.1,
           size = 2) +
  annotate("text", 
           x = quantiles[2], 
           y = max(ggplot_build(hist_plot)$data[[1]]$count) * 0.95,
           label = "Top \n 25%",
           color = "red",
           hjust = -0.1,
           size = 2) +
  scale_fill_viridis_c(option = "D", guide = "none") +
  labs(x = "Soil N (cg/kg)", y = "Count") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black", size = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  )

## ─────────────────────────────────────────────────────────
##  C) COMBINE & SAVE
## ─────────────────────────────────────────────────────────
combined_plot <- main_map +
  inset_element(
    hist_plot,
    left   = 0.02,
    bottom = 0.02,
    right  = 0.42,   # bigger
    top    = 0.44
  )

quartz()
combined_plot
ggsave(
  "Nitrogen_Map_Landraces_Maize_Romerro.png",
  plot   = combined_plot,
  width  = 8,
  height = 6,
  dpi    = 600,
  bg     = "white"
)
getwd()



GAPIT.Manhattan()





