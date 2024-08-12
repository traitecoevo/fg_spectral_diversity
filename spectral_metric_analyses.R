# refactor 1

# load packages -----------------------------------------------------------
library(raster)
library(tidyverse) # dplyr masks raster::extract + raster::select
library(vegan) # for calculate_sv
library(geometry) # for calculate_chv
library(sf)
source('funx.R')

setwd("C:/Users/adele/Documents/fg_spectral_diversity")

# combining spectral band images into one tif file  -------------------------------------------------------

## folders in mosaics_dir should be carry AUSPLOT plot_id (e.g. NSABHC0009)
# file directory
sixteen_dir <- "data/2016_mosaics"

twentyfour_dir <- 'data/2024_mosaics'

create_multiband_image(sixteen_dir, c("green", "red", "red_edge", "nir"))

create_multiband_image(twentyfour_dir, c("blue", "green", "red", "red_edge", "nir"))


# creating mask for raster using biodivmapR --------------------------------------------------------------

remotes::install_github('cran/dissUtils')
remotes::install_github('jbferet/biodivMapR', force = T)
library(biodivMapR)
library(terra)
library(tools)
library(raster)
library(stars)

envi_dir <- 'data_out/combined_rasters'

mask_dir <- 'C:/Users/adele/Documents/fg_spectral_diversity/biodivmapR/RESULTS/'

masked_raster_dir <- 'C:/Users/adele/Documents/fg_spectral_diversity/data_out/combined_rasters/masked/'


create_masked_raster(envi_dir, mask_dir, masked_raster_dir, '2024',
                     c("blue", "green", "red", "red_edge", "nir"), c(450,560,650,730,840),
                     NDVI_Thresh = 0.05, NIR_Thresh = 0.02,
                     Blue = 450, Red = 650, NIR = 840)


# 2016 raster
create_masked_raster(envi_dir, mask_dir, masked_raster_dir, '2016',
                     c("green", "red", "red_edge", "nir"), c(550,660,735,790),
                     NDVI_Thresh = 0.02, NIR_Thresh = 0.02,
                     Blue = 450, Red = 660, NIR = 790)


# extract pixel values for all rasters by subplot  ------------------------------------------------------

subplot_files <- list.files('data/fishnets', pattern = '_fishnet.shp$', full.names = TRUE)

# FOR MASKED IMAGES
raster_files_masked_16 <- list.files('data_out/combined_rasters/masked/2016', pattern = '_combined_image_masked.tif$', full.names = TRUE)

masked_pixel_values_list_16 <- extract_pixels_values(raster_files_masked_16, subplot_files, c('green', 'red', 'red_edge', 'nir'))
masked_pixel_values_16 <- bind_rows(masked_pixel_values_list_16, .id = 'identifier')

raster_files_masked_24 <- list.files('data_out/combined_rasters/masked/2024', pattern = '_combined_image_masked.tif$', full.names = TRUE)

masked_pixel_values_list_24 <- extract_pixels_values(raster_files_masked_24, subplot_files, c('blue', 'green', 'red', 'red_edge', 'nir'))
masked_pixel_values_24 <- bind_rows(masked_pixel_values_list_24, .id = 'identifier')

# FOR UNMASKED IMAGES
raster_files_unmasked_16 <- list.files('data_out/combined_rasters/2016', pattern = '_combined_image.tif$', full.names = TRUE)

unmasked_pixel_values_list_16 <- extract_pixels_values(raster_files_unmasked_16, subplot_files, c('green', 'red', 'red_edge', 'nir'))
unmasked_pixel_values_16 <- bind_rows(unmasked_pixel_values_list_16, .id = 'identifier')

raster_files_unmasked_24 <- list.files('data_out/combined_rasters/2024', pattern = '_combined_image.tif$', full.names = TRUE)

unmasked_pixel_values_list_24 <- extract_pixels_values(raster_files_unmasked_24, subplot_files, c('blue', 'green', 'red', 'red_edge', 'nir'))
unmasked_pixel_values_24 <- bind_rows(unmasked_pixel_values_list_24, .id = 'identifier')
head(unmasked_pixel_values_24)

# apply spectral metrics to pixel value df --------------------------------

source('funx.R')

masked_metrics_24 <- calculate_metrics(masked_pixel_values_24, masked = TRUE, c('blue', 'green', 'red', 'red_edge', 'nir'))
unmasked_metrics_24 <- calculate_metrics(unmasked_pixel_values_24, masked = FALSE, c('blue', 'green', 'red', 'red_edge', 'nir'))

metrics_24 <- bind_rows(unmasked_metrics_24, masked_metrics_24)

masked_metrics_16 <- calculate_metrics(masked_pixel_values_16, masked = TRUE, c('green', 'red', 'red_edge', 'nir'))
unmasked_metrics_16 <- calculate_metrics(unmasked_pixel_values_16, masked = FALSE, c('green', 'red', 'red_edge', 'nir'))

metrics_16 <- bind_rows(unmasked_metrics_16, masked_metrics_16)

#write.csv(metrics_24, 'data_out/2024_spec_div_values.csv')
#write.csv(metrics_16, 'data_out/2016_spec_div_values.csv')

metrics_24 <- read_csv('data_out/2024_spec_div_values.csv')

# field observations for every plot ---------------------------------------
library(devtools)
install_github("ternaustralia/ausplotsR", build_vignettes = TRUE, dependencies = TRUE)
library(ausplotsR)

# read survey data
twentyfour_survey_data <- read.csv('data/ausplots_march_24.csv')

OR

# for 2016 data

plots_oi <- c('NSABHC0009', 'NSABHC0010', 'NSABHC0011', 'NSABHC0012')

veg <- get_ausplots(plots_oi, veg.vouchers = T, veg.PI = T)

sixteen_survey_data <- veg$veg.PI %>%
  left_join(veg$site.info %>% select(site_unique, visit_start_date), by = "site_unique") %>%
  filter(substr(visit_start_date, 1, 4) == '2016')

twentyfour_field_diversity <- calculate_field_diversity(twentyfour_survey_data)

sixteen_field_diversity <- calculate_field_diversity(sixteen_survey_data)

# sanity check one of these 
head(final_results)



# erghhhhhhhhhhhhhhhhhhhh ... need to still clean up the survey data and think about stuff like... 
# dead shrub stuff, things IDd to genus level etc etc....
#community_matrices[["NSABHC012"]] |> View()


# combined field and spectral metrics -------------------------------------

library(gridExtra)

tax_and_spec_diversity_values_24 <- left_join(twentyfour_field_diversity, metrics_24, by = c('site' = 'identifier', 'subplot_id'))
tax_and_spec_diversity_values_16 <- left_join(sixteen_field_diversity, metrics_16, by = c('site' = 'identifier', 'subplot_id'))

OR
#tax_and_spec_diversity_values_24 <- read_csv("data_out/tax_and_spec_div_values_24.csv")
#tax_and_spec_diversity_values_16 <- read_csv('data_out/tax_and_spec_div_values_16.csv')

# pivot longer 
tax_spec_div_long <- tax_and_spec_diversity_values_24 %>%
  pivot_longer(cols = c(species_richness, shannon_diversity, simpson_diversity, pielou_evenness),
               names_to = "taxonomic_metric",
               values_to = "taxonomic_value") %>%
  pivot_longer(cols = c(CV, SV, CHV),
               names_to = "spectral_metric",
               values_to = "spectral_value")

taxonomic_labels <- c(
  species_richness = "Species\nRichness",
  shannon_diversity = "Shannon\nDiversity",
  simpson_diversity = "Simpson\nDiversity",
  pielou_evenness = "Pielou's\nEvenness"
)
# plot, showing spectral ~ taxonomic relos. add , color = masked and unmasked
unmasked_24_plot <- tax_spec_div_long %>%
  filter(image_type == 'unmasked') %>%
ggplot(aes(x = spectral_value, y = taxonomic_value, color = site)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = 'red') +
  facet_grid(taxonomic_metric ~ spectral_metric, scales = "free", labeller = labeller(taxonomic_metric = taxonomic_labels)) +
  labs(
    title = "unmasked 2024",
    x = "spectral diversity value",
    y = "taxonomic diversity value") +
  theme_minimal() +
  scale_color_manual(values = c('lightpink1', 'violet', 'skyblue1', 'springgreen3')) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.y = element_text(angle = 0, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

masked_24_plot <- tax_spec_div_long %>%
  filter(image_type == 'masked') %>%
  ggplot(aes(x = spectral_value, y = taxonomic_value, color = site)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = 'red') +
  facet_grid(taxonomic_metric ~ spectral_metric, scales = "free", labeller = labeller(taxonomic_metric = taxonomic_labels)) +
  labs(
    title = "masked 2024",
    x = "spectral diversity value",
    y = "taxonomic diversity value") +
  theme_minimal() +
  scale_color_manual(values = c('lightpink1', 'violet', 'skyblue1', 'springgreen3')) +
  
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.y = element_text(angle = 0, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 


grid.arrange(unmasked_24_plot, masked_24_plot, ncol = 2)


# for 2016

# pivot longer 
tax_spec_div_long_16 <- tax_and_spec_diversity_values_16 %>%
  pivot_longer(cols = c(species_richness, shannon_diversity, simpson_diversity, pielou_evenness),
               names_to = "taxonomic_metric",
               values_to = "taxonomic_value") %>%
  pivot_longer(cols = c(CV, SV, CHV),
               names_to = "spectral_metric",
               values_to = "spectral_value")


masked_16_plot <- tax_spec_div_long_16 %>%
  filter(image_type == 'masked') %>%
  ggplot(aes(x = spectral_value, y = taxonomic_value, color = site)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = 'red') +
  facet_grid(taxonomic_metric ~ spectral_metric, scales = "free", labeller = labeller(taxonomic_metric = taxonomic_labels)) +
  labs(
    title = "masked 2016",
    x = "spectral diversity value",
    y = "taxonomic diversity value") +
  theme_minimal() +
  scale_color_manual(values = c('lightpink1', 'springgreen3')) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.y = element_text(angle = 0, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

unmasked_16_plot <- tax_spec_div_long_16 %>%
  filter(image_type == 'unmasked') %>%
  ggplot(aes(x = spectral_value, y = taxonomic_value, color = site)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = 'red') +
  facet_grid(taxonomic_metric ~ spectral_metric, scales = "free", labeller = labeller(taxonomic_metric = taxonomic_labels)) +
  labs(
    title = "unmasked 2016",
    x = "spectral diversity value",
    y = "taxonomic diversity value") +
  theme_minimal() +
  scale_color_manual(values = c('lightpink1', 'springgreen3')) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.y = element_text(angle = 0, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 


combined_plot <- grid.arrange(unmasked_24_plot, masked_24_plot, unmasked_16_plot, masked_16_plot, ncol = 2)


