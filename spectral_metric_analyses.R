# refactor 1

# load packages -----------------------------------------------------------
library(raster)
library(tidyverse) # dplyr masks raster::extract + raster::select
library(vegan) # for calculate_sv
library(geometry) # for calculate_chv
library(sf)
source('funx.R')

# combining spectral band images into one tif file  -------------------------------------------------------

# file directory
mosaics_dir <- "data/2024_mosaics"

# sources function from funx.R
create_multiband_image(mosaics_dir)

# creating mask for raster using biodivmapR --------------------------------------------------------------

remotes::install_github('cran/dissUtils')
remotes::install_github('jbferet/biodivMapR', force = T)
library(biodivMapR)
library(terra)
library(tools)
library(raster)
library(stars)

raster_files <- list.files('data_out/combined_rasters', pattern = '_reflectance_combined_image.tif$', full.names = T)

envi_dir <- 'C:/Users/adele/Documents/fg_spectral_diversity/ENVI/'

mask_dir <- 'C:/Users/adele/Documents/fg_spectral_diversity/biodivmapR/RESULTS/'

masked_raster_dir <- 'C:/Users/adele/Documents/fg_spectral_diversity/data_out/combined_rasters/masked/'

# function sourced from funx.R 
create_masked_raster(raster_files, envi_dir, mask_dir, masked_raster_dir)

# extract pixel values for all rasters by subplot  ------------------------------------------------------

subplot_files <- list.files('data/fishnets', pattern = '_fishnet.shp$', full.names = TRUE)

# FOR MASKED IMAGES
raster_files_masked <- list.files('data_out/combined_rasters/masked', pattern = '_reflectance_combined_image_masked.tif$', full.names = TRUE)

masked_pixel_values_list <- extract_pixels_values(raster_files_masked, subplot_files)
masked_pixel_values <- bind_rows(masked_pixel_values_list, .id = 'identifier')

# FOR UNMASKED IMAGES
raster_files_unmasked <- list.files('data_out/combined_rasters', pattern = '_reflectance_combined_image.tif$', full.names = TRUE)

unmasked_pixel_values_list <- extract_pixels_values(raster_files_unmasked, subplot_files)
unmasked_pixel_values <- bind_rows(final_pixel_values_list, .id = 'identifier')


# apply spectral metrics to pixel value df --------------------------------

masked_metrics <- calculate_metrics(masked_pixel_values, masked = TRUE)
unmasked_metrics <- calculate_metrics(unmasked_pixel_values, masked = FALSE)

metrics <- bind_rows(unmasked_metrics, masked_metrics) 

# field observations for every plot ---------------------------------------

# read survey data
survey_data <- read.csv('data/ausplots_march_24.csv')

# get unique site names
unique_sites <- unique(survey_data$site_unique) 
unique_sites <- unique_sites[unique_sites != ""]

# list to store results for all lists
all_site_results <- list()

# list to store community matrices- this is a temporary step, to check that community matrices are correct :)
community_matrices <- list()

# Loop through each unique site
for (site in unique_sites) {
  # Filter data for the current site
  site_survey_data <- survey_data %>%
    filter(site_unique == site)
  
  # Extract only direction of the transect (no numbers)
  site_survey_data$transect_direction <- gsub('[[:digit:]]+', '', site_survey_data$transect)
  
  # Extract only number of the transect (no direction)
  site_survey_data$transect_number <- as.numeric(gsub(".*?([0-9]+).*", "\\1", site_survey_data$transect))
  
  # Create variable for fixed transect direction (to order them all transects in the same direction)
  site_survey_data$transect_direction2 <- NA 
  
  # Create variable for fixed point number (inverse in some cases as if they had been collected in the same direction)
  site_survey_data$point_number2 <- NA 
  
  # Create XY empty variables for plot XY coordinates
  site_survey_data$X_plot <- NA 
  site_survey_data$Y_plot <- NA
  
  # For loop to homogenize transects and numbers. It converts all W-E to E-W and all N-S to S-N 
  for (i in 1:nrow(site_survey_data)){
    if (site_survey_data[i, "transect_direction"] == "W-E") {
      site_survey_data[i, "point_number2"] <- 100 - site_survey_data[i, "point_number"] # If transect E-W, transect fixed is W-E and inverse numbers
      site_survey_data[i, "transect_direction2"] <- "E-W"
    }
    if (site_survey_data[i, "transect_direction"] == "E-W") {
      site_survey_data[i, "point_number2"] <- site_survey_data[i, "point_number"] # If transect W-E, all stays the same
      site_survey_data[i, "transect_direction2"] <- "E-W"
    }
    if (site_survey_data[i, "transect_direction"] == "S-N") {
      site_survey_data[i, "point_number2"] <- site_survey_data[i, "point_number"] # If transect N-S, all stays the same
      site_survey_data[i, "transect_direction2"] <- "S-N"
    }
    if (site_survey_data[i, "transect_direction"] == "N-S") {
      site_survey_data[i, "point_number2"] <- 100 - site_survey_data[i, "point_number"] # If transect S-N, transect fixed is N-S and inverse numbers
      site_survey_data[i, "transect_direction2"] <- "S-N"
    }
  }
  
  # For loop to assign plotXY coordinates to each point intercept
  for (i in 1:nrow(site_survey_data)){
    if (site_survey_data[i, "transect_direction2"] == "E-W") {
      if (site_survey_data[i, "transect_number"] == 1){
        site_survey_data[i, "Y_plot"] <- 10
        site_survey_data[i, "X_plot"] <- site_survey_data[i, "point_number2"]
      }
      if (site_survey_data[i, "transect_number"] == 2){
        site_survey_data[i, "Y_plot"] <- 30
        site_survey_data[i, "X_plot"] <- site_survey_data[i, "point_number2"]
      }
      if (site_survey_data[i, "transect_number"] == 3){
        site_survey_data[i, "Y_plot"] <- 50
        site_survey_data[i, "X_plot"] <- site_survey_data[i, "point_number2"]
      }
      if (site_survey_data[i, "transect_number"] == 4){
        site_survey_data[i, "Y_plot"] <- 70
        site_survey_data[i, "X_plot"] <- site_survey_data[i, "point_number2"]
      }
      if (site_survey_data[i, "transect_number"] == 5){
        site_survey_data[i, "Y_plot"] <- 90
        site_survey_data[i, "X_plot"] <- site_survey_data[i, "point_number2"]
      }
    }
    if (site_survey_data[i, "transect_direction2"] == "S-N") {
      if (site_survey_data[i, "transect_number"] == 1){
        site_survey_data[i, "X_plot"] <- 10
        site_survey_data[i, "Y_plot"] <- site_survey_data[i, "point_number2"]
      }
      if (site_survey_data[i, "transect_number"] == 2){
        site_survey_data[i, "X_plot"] <- 30
        site_survey_data[i, "Y_plot"] <- site_survey_data[i, "point_number2"]
      }
      if (site_survey_data[i, "transect_number"] == 3){
        site_survey_data[i, "X_plot"] <- 50
        site_survey_data[i, "Y_plot"] <- site_survey_data[i, "point_number2"]
      }
      if (site_survey_data[i, "transect_number"] == 4){
        site_survey_data[i, "X_plot"] <- 70
        site_survey_data[i, "Y_plot"] <- site_survey_data[i, "point_number2"]
      }
      if (site_survey_data[i, "transect_number"] == 5){
        site_survey_data[i, "X_plot"] <- 90
        site_survey_data[i, "Y_plot"] <- site_survey_data[i, "point_number2"]
      }
    }
  }
  
  # subplot rows and columns - +1 ensures 0 point values fall into correct subplot, 
  # pmin ensures 100 point values falls in correct subplot given +1
  site_survey_data$subplot_row <- pmin(ceiling((site_survey_data$Y_plot + 1) / 20), 5)
  site_survey_data$subplot_col <- pmin(ceiling((site_survey_data$X_plot + 1) / 20), 5)
  
  # single ID for subplot row and column
  site_survey_data$subplot_id <- paste(site_survey_data$subplot_row, site_survey_data$subplot_col, sep = "_")
  
  subplot_diversity <- site_survey_data %>%
    drop_na(standardised_name) %>%
    group_by(subplot_id) %>%
    summarise(species_richness = n_distinct(standardised_name))
  
  community_matrix <- site_survey_data %>%
    drop_na(standardised_name) %>%
    count(subplot_id, standardised_name) %>%
    spread(standardised_name, n, fill = 0)
  
  # store the community matrix in the list - this is a temp step to check!!!
  community_matrices[[site]] <- community_matrix
  
  # remove unwanted column if it exists -- WHY DOES THIS COLUMN EXIST~!!!!>???>
  if ("V1" %in% colnames(community_matrix)) {
    community_matrix <- community_matrix %>% 
      dplyr::select(-V1)
  }
  
  # calculate diversity indices
  shannon_diversity <- diversity(community_matrix[, -1], index = "shannon")
  simpson_diversity <- diversity(community_matrix[, -1], index = "simpson")
  
  subplot_diversity <- subplot_diversity %>%
    mutate(shannon_diversity = shannon_diversity,
      simpson_diversity = simpson_diversity,
      site = site)
  
  # store  result for  current site
  all_site_results[[site]] <- subplot_diversity
}

# combine into one df
final_results <- bind_rows(all_site_results, .id = "site")

# add identifier based on site name so that this can be joined with spectral results 
# in refactor two you could remove this step by just calling the sites by their ausplot name...
# would require changing the file names of OG rasters
final_results <- final_results %>%
  mutate(identifier = recode(site,
                       "NSABHC009" = "emu",
                       "NSABHC010" = "emgrazed",
                       "NSABHC011" = "sandstone",
                       "NSABHC012" = "cons"))


# sanity check one of these 
head(final_results)


# erghhhhhhhhhhhhhhhhhhhh ... need to still clean up the survey data and think about stuff like... 
# dead shrub stuff, things IDd to genus level etc etc....
#community_matrices[["NSABHC012"]] |> View()

# combined field and spectral metrics -------------------------------------

tax_and_spec_diversity_values <- left_join(final_results, metrics, by = c('identifier', 'subplot_id'))

# pivot longer 
tax_spec_div_long <- tax_and_spec_diversity_values %>%
  pivot_longer(cols = c(species_richness, shannon_diversity, simpson_diversity),
               names_to = "taxonomic_metric",
               values_to = "taxonomic_value") %>%
  pivot_longer(cols = c(CV, SV, CHV),
               names_to = "spectral_metric",
               values_to = "spectral_value")


# plot, showing spectral ~ taxonomic relos. add , color = masked and unmasked
ggplot(tax_spec_div_long, aes(x = spectral_value, y = taxonomic_value, color = image_type)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, size = 1.5, linetype = "solid", aes(fill = image_type), color = "black") +
  facet_grid(taxonomic_metric ~ spectral_metric, scales = "free") +
  labs(
    title = "comparison of taxonomic + spectral diversity metrics ",
    x = "spectral diversity value",
    y = "taxonomic diversity value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
