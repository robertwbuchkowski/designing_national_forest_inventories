# Script to create a map of the plots:
# A: Robert W. Buchkowski

# Load packages:
library(lubridate)
library(lme4)
library(sf)
library(raster)
library(terra)
library(tidyverse)

# Load in the necessary data:
source("Scripts/00_load_data.R")

# Load in the map data:

# Function to convert each row to WGS84
convert_row_to_wgs84 <- function(row) {
  utm_crs <- paste0("+proj=utm +zone=", row$utm_zone, " +datum=WGS84 +units=m +no_defs")
  sf_point <- st_as_sf(row, coords = c("utm_e", "utm_n"), crs = utm_crs)
  st_transform(sf_point, crs = 4326)
}

# Apply the function to each row and combine results
# NOTE: Take a few seconds to convert the coordinates!
sf_points_wgs84 <- sitedata %>%
  rowwise() %>%
  group_split() %>%
  lapply(convert_row_to_wgs84) %>%
  do.call(rbind, .)

# Identify the NFI plots with data:
sf_points_wgs84 = sf_points_wgs84 %>%
  mutate(Mineral = ifelse(nfi_plot %in% unique(Mgpdata$nfi_plot), "Yes","No"),
         Organic = ifelse(nfi_plot %in% unique(Ogpdata$nfi_plot), "Yes","No")) %>%
  filter(!(Mineral == "No" & Organic == "No"))

rmda = gpremeas_calc_all %>%
  select(nfi_plot) %>%
  distinct() %>% pull()

sf_points_wgs84 %>%
  select(nfi_plot, Mineral, Organic) %>%
  distinct() %>%
  mutate(Remeasure = ifelse(nfi_plot %in% rmda, "Yes", "No")) %>%
  select(nfi_plot, Mineral, Organic, Remeasure) %>%
  write_sf("Data/plot_shapefile/plots.shp")
