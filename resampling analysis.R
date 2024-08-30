
# RESAMPLING Yo! ----------------------------------------------------------

library(terra)
library(pROC)

library(terra)
# List and load raster files
r_2024 <- list.files('data_out/combined_rasters/2024', full.names = TRUE)
rast_2024 <- lapply(r_2024, rast)

# Initialize a list to store the resampled rasters
rasters_resampled <- list()

# Define the resolutions for each raster
resolutions <- list(
  c(0.07125, 0.07125),  # Resolution for the first raster
  c(0.07125, 0.07125),  # Resolution for the second raster
  c(0.08334, 0.08334),  # Resolution for the third raster
  c(0.08334, 0.08334)   # Resolution for the fourth raster
)

# Loop through each raster
for (i in seq_along(rast_2024)) {
  # Get the current raster
  current_raster <- rast_2024[[i]]
  
  # Get the resolution for the current raster
  new_res <- resolutions[[i]]
  
  # Create a template raster for the current raster with the new resolution
  r_template <- rast(ncol = ceiling((ext(current_raster)[2] - ext(current_raster)[1]) / new_res[1]),
                     nrow = ceiling((ext(current_raster)[4] - ext(current_raster)[3]) / new_res[2]),
                     ext = ext(current_raster),
                     crs = crs(current_raster))
  
  # Set the resolution of the template raster
  res(r_template) <- new_res
  
  # Resample the current raster
  r_resampled <- resample(current_raster, r_template, method = 'bilinear')
  
  # Store the resampled raster in the list
  rasters_resampled[[i]] <- r_resampled
}

output_files <- gsub("\\.tif$", "_RESAMPLED2.tif", basename(r_2024))
output_files <- file.path('data_out/combined_rasters/resampled_2024', output_files)
for (i in seq_along(rasters_resampled)) {
  writeRaster(rasters_resampled[[i]], filename = output_files[i], overwrite = TRUE)
}


# create 2024 resampled masked rasters
create_masked_raster(input = 'data_out/combined_rasters/resampled_2024',
                     output_dir = 'data_out/combined_rasters/masked/resampled_2024',
                     band_wavelengths = c(450,560,650,730,840),
                     NDVI_Thresh_df = ndvi_threshold_df_24,
                     NIR_Thresh_df = nir_threshold_df_24,
                     Red_band = 650,
                     NIR_band = 840)  

#masked_r <- list.files('data_out/combined_rasters/masked/resampled_2024', full.names = T)
#masked_r <- lapply(masked_r, rast) 
#plot(masked_r[[1]])

#r <- list.files('data_out/combined_rasters/masked/2024', pattern = '\\.tif$', full.names = T)
#r <- lapply(r, rast)
#plot(r[[1]])

#remove blue layer 
remove_blue <- function(raster_file){
  
  #read in stack
  r <- rast(raster_file)
  
  #remove first layer
  r_modified <- r[[2:nlyr(r)]]
  
  new_raster_path <- gsub(".tif", "_no_blue.tif", raster_file)
  
  writeRaster(r_modified, filename = new_raster_path, overwrite = TRUE)
}

lapply(raster_files_masked_24_RS, remove_blue)
lapply(raster_files_unmasked_24_RS,remove_blue)

subplot_files <- list.files('data/fishnets', pattern = '_fishnet.shp$', full.names = TRUE)

raster_files_masked_24_RS <- list.files('data_out/combined_rasters/masked/resampled_2024', pattern = 'no_blue.tif$', full.names = TRUE)
raster_files_unmasked_24_RS <- list.files('data_out/combined_rasters/resampled_2024', pattern = 'no_blue.tif$', full.names = TRUE)
masked_pv_24_rs <- extract_pixel_values(raster_files_masked_24_RS, subplot_files, c('green', 'red', 'red_edge', 'nir'))
unmasked_pv_24_rs <- extract_pixel_values(raster_files_unmasked_24_RS, subplot_files, c('green', 'red', 'red_edge', 'nir'))

tic()
rs_masked_metrics_24 <- calculate_metrics(masked_pv_24_rs, masked = TRUE, c('green', 'red', 'red_edge', 'nir'))
toc()

tic()
rs_unmasked_metrics_24 <- calculate_metrics(unmasked_pv_24_rs, masked = FALSE, c('green', 'red', 'red_edge', 'nir'))
toc()

rs_metrics_24 <- bind_rows(rs_unmasked_metrics_24, rs_masked_metrics_24)

rs_tax_and_spec_div_values <- left_join(twentyfour_field_diversity$final_results, rs_metrics_24, by = c('site' = 'identifier', 'subplot_id'))

write.csv(rs_tax_and_spec_div_values, 'data_out/rs_24_tax_and_spec_div_values.csv')



# plotting ----------------------------------------------------------------

tax_spec_div_long_rs <- rs_tax_and_spec_div_values %>%
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
unmasked_24_plot_rs <- tax_spec_div_long_rs %>%
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

masked_24_plot_rs <- tax_spec_div_long_rs %>%
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


# for 2016

# pivot longer 
tax_spec_div_long_16 <- tax_and_spec_diversity_values_16 %>%
  pivot_longer(cols = c(species_richness, shannon_diversity, simpson_diversity, pielou_evenness),
               names_to = "taxonomic_metric",
               values_to = "taxonomic_value") %>%
  pivot_longer(cols = c(CV, SV, CHV, CHV_nopca),
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


rs_combined_plot <- grid.arrange(unmasked_24_plot_rs, masked_24_plot_rs, unmasked_16_plot, masked_16_plot, ncol = 2)


