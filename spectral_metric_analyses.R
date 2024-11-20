
# load packages -----------------------------------------------------------
setwd("C:/Users/adele/Documents/fg_spectral_diversity")

library(raster)
library(tidyverse) # dplyr masks raster::extract + raster::select
library(tools)
library(vegan) # for calculate_sv and chv
library(geometry) # for calculate_chv
library(sf)
library(pROC) # for optimum thresholds
library(terra)
#library(devtools)
#install_github("ternaustralia/ausplotsR", build_vignettes = TRUE, dependencies = TRUE)
library(ausplotsR) # for field div
source('funx.R')

# combining spectral band images into one tif file  -------------------------------------------------------

## folders in mosaics_dir should be carry AUSPLOT plot_id (e.g. NSABHC0009)
# file directory
sixteen_dir <- "data/2016_mosaics"

twentyfour_dir <- 'data/2024_mosaics'

create_multiband_image(sixteen_dir, c("green", "red", "red_edge", "nir"))

create_multiband_image(twentyfour_dir, c("blue", "green", "red", "red_edge", "nir"))

# calculating optimum thresholds for masking ------------------------------
## 2024
# ndvi 

ndvi_values_24 <- read_csv("data/ndvi_2024_values.csv")

ndvi_threshold_df_24 <- find_optimum_thresholds(ndvi_values_24, class = 'class', value = 'ndvi', site = 'location', class_value = 'ground')

# nir
nir_values_24 <- read_csv("data/nir_less_05_table.csv")

nir_threshold_df_24 <- find_optimum_thresholds(nir_values_24, class = 'class', value = 'nir', class_value = 'shadow', site = 'location')

# havent done random points for emu grazed, 0.045 looks about right haha
nir_threshold_df_24 <- rbind(nir_threshold_df_24, data.frame(site = "NSABHC0010", threshold = 0.045))

red_threshold_df_24 <- data.frame(site = c('NSABHC0009', 'NSABHC0010', 'NSABHC0011', 'NSABHC0012'),
                                  threshold = c(0.07, 0.08, 0.07, 0.09))
# 2016 
# ndvi
ndvi_values_16 <- read_csv("data/ndvi_2016_values.csv")

ndvi_threshold_df_16 <- find_optimum_thresholds(ndvi_values_16, class = 'class', value = 'ndvi', site = 'site')

# NIR - havent done random points arbitrarily chosen these values (for now)
nir_threshold_df_16 <- data.frame(
  site = c("NSABHC0009", "NSABHC0012"),
  threshold = c(0.15, 0.20))

# creating mask for raster --------------------------------------------------------------

# create 2024 masked rasters
create_masked_raster(input = 'data_out/combined_rasters/2024',
                     output_dir = 'data_out/combined_rasters/masked/2024',
                     NDVI_Thresh_df = ndvi_threshold_df_24,
                     NIR_Thresh_df = nir_threshold_df_24,
                     Red_Thresh = 0.1)  

# create 2016 masked rasters
create_masked_raster(input = 'data_out/combined_rasters/2016',
                     output_dir = 'data_out/combined_rasters/masked/2016',
                     band_wavelengths = c(550,660,735,790),
                     NDVI_Thresh_df = ndvi_threshold_df_16,
                     NIR_Thresh_df = nir_threshold_df_16,
                     Red_band = 660,
                     NIR_band = 790)  


# extract pixel values for all rasters by subplot  ------------------------------------------------------

subplot_files <- list.files('data/fishnets', pattern = '_fishnet.shp$', full.names = TRUE)

# FOR MASKED IMAGES
# 2016
raster_files_masked_16 <- list.files('data_out/combined_rasters/masked/2016', pattern = '_combined_image_masked.tif$', full.names = TRUE)

masked_pixel_values_16 <- extract_pixel_values(raster_files_masked_16, subplot_files, c('green', 'red', 'red_edge', 'nir'))

#2024
raster_files_masked_24 <- list.files('data_out/combined_rasters/masked/2024', pattern = '_combined_image_masked.tif$', full.names = TRUE)

masked_pixel_values_24 <- extract_pixel_values(raster_files_masked_24, subplot_files, c('blue', 'green', 'red', 'red_edge', 'nir'))

# FOR UNMASKED IMAGES
#2016
raster_files_unmasked_16 <- list.files('data_out/combined_rasters/2016', pattern = '_combined_image.tif$', full.names = TRUE)

unmasked_pixel_values_16 <- extract_pixel_values(raster_files_unmasked_16, subplot_files, c('green', 'red', 'red_edge', 'nir'))

#2024
raster_files_unmasked_24 <- list.files('data_out/combined_rasters/2024', pattern = '_combined_image.tif$', full.names = TRUE)

unmasked_pixel_values_24 <- extract_pixel_values(raster_files_unmasked_24, subplot_files, c('blue', 'green', 'red', 'red_edge', 'nir'))


# calculate 'min_points' (the minimum number of pixels in a given subplot)
pixel_per_plot <- masked_pixel_values_16 %>%
  na.omit() %>%
  group_by(identifier, subplot_id) %>%
  summarise(count = n()) 

min(pixel_per_plot$count)
#8691 for 2016 masked 9836 for unmasked 2016
#90361 for 2024 masked 1271009 unmasked


# apply spectral metrics to pixel value df --------------------------------

masked_metrics_24 <- calculate_spectral_metrics(masked_pixel_values_24, 
                                                masked = TRUE, 
                                                c('blue', 'green', 'red', 'red_edge', 'nir'))

unmasked_metrics_24 <- calculate_spectral_metrics(unmasked_pixel_values_24, 
                                         masked = FALSE, 
                                         c('blue', 'green', 'red', 'red_edge', 'nir'))

metrics_24 <- bind_rows(unmasked_metrics_24, masked_metrics_24)

masked_metrics_16 <- calculate_spectral_metrics(masked_pixel_values_16,
                                                masked = T,
                                                c('green', 'red', 'red_edge', 'nir'),
                                                min_points = 8691)

unmasked_metrics_16 <- calculate_spectral_metrics(unmasked_pixel_values_16,
                                                masked = F,
                                                c('green', 'red', 'red_edge', 'nir'),
                                                min_points = 9836)

metrics_16 <- bind_rows(unmasked_metrics_16, masked_metrics_16)


# field observations for every plot ---------------------------------------
library(devtools)
install_github("ternaustralia/ausplotsR", build_vignettes = TRUE, dependencies = TRUE)
library(ausplotsR)
library(vegan)
library(tidyverse)

source('funx.R')
# read survey data
twentyfour_survey_data <- read.csv('data/ausplots_march_24.csv')

# for 2016 data

plots_oi <- c('NSABHC0009', 'NSABHC0010', 'NSABHC0011', 'NSABHC0012')

veg <- get_ausplots(plots_oi, veg.vouchers = T, veg.PI = T)

sixteen_survey_data <- veg$veg.PI %>%
  left_join(veg$site.info %>% dplyr::select(site_unique, visit_start_date), by = "site_unique") %>%
  filter(substr(visit_start_date, 1, 4) == '2016')

sixteen_field_diversity <- calculate_field_diversity(sixteen_survey_data)


# combined field and spectral metrics -------------------------------------

#metrics_24 <- read_csv('data_out/2024_spec_div_values.csv')
#metrics_16 <- read_csv('data_out/2016_spec_div_values.csv')


tax_and_spec_diversity_values_24 <- left_join(twentyfour_field_diversity$final_results, metrics_24, by = c('site' = 'identifier', 'subplot_id'))
tax_and_spec_diversity_values_16 <- left_join(sixteen_field_diversity$final_results, metrics_16, by = c('site' = 'identifier', 'subplot_id'))


tax_and_spec_diversity_values_24$year <- 2024 
tax_and_spec_diversity_values_16$year <- 2016

tax_and_spec_div_values <- rbind(tax_and_spec_diversity_values_16, tax_and_spec_diversity_values_24)
tax_and_spec_div_values <- tax_and_spec_div_values %>%
  na.omit() # remove rows where diversity couldnt be calculated as sr <2 and
            # where spec div cant be calc (2016 cons cropped bit)


#write.csv(tax_and_spec_div_values, 'data_out/tax_and_spec_diversity_values.csv', row.names = F)
