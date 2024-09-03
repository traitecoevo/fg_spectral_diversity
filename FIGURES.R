# FIGURES, BABY! 
setwd("C:/Users/adele/Documents/fg_spectral_diversity")

library(reshape2)
library(ggplot2)
source('funx.R')
library(tools)
library(tidyverse)
library(sf)


# main figure ------------------------------------------------------------

data_frames <- list(
  tax_and_spec_div_values = tax_and_spec_div_values,
  alive_tax_and_spec_div_values = alive_tax_and_spec_div_values,
  rs_tax_and_spec_div_values = rs_tax_and_spec_div_values
)

create_plots <- function(df) {
  # taxonomic labels in correct order
  taxonomic_labels <- c(
    species_richness = "Species\nRichness",
    shannon_diversity = "Shannon\nDiversity",
    exp_shannon = "Exponential\nShannons",
    simpson_diversity = "Simpson\nDiversity",
    inv_simpson = "Inverse\nSimpsons",
    pielou_evenness = "Pielou's\nEvenness"
  )
  
  # pivot longer
  df_long <- df %>%
    pivot_longer(cols = c(species_richness, shannon_diversity, simpson_diversity, pielou_evenness, exp_shannon, inv_simpson),
                 names_to = "taxonomic_metric",
                 values_to = "taxonomic_value") %>%
    pivot_longer(cols = c(CV, SV, CHV),
                 names_to = "spectral_metric",
                 values_to = "spectral_value") %>%
    mutate(taxonomic_metric = factor(taxonomic_metric, levels = names(taxonomic_labels)))
  
  # 2024 plots
  unmasked_plot_24 <- df_long %>%
    filter(image_type == 'unmasked', year == 2024) %>%
    ggplot(aes(x = spectral_value, y = taxonomic_value, color = site)) +
    geom_point(alpha = 0.7, show.legend = F) +
    geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = 'black') +
    facet_grid(taxonomic_metric ~ spectral_metric, scales = "free", labeller = label_value) +
    labs(
      title = "Unmasked Rasters",
      y = '2024') +
    theme_minimal() +
    scale_color_manual(values = c('darkgreen', 'chocolate', 'navajowhite2', 'darkseagreen')) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.y.right = element_blank())
  
  masked_plot_24 <- df_long %>%
    filter(image_type == 'masked', year == 2024) %>%
    ggplot(aes(x = spectral_value, y = taxonomic_value, color = site)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = 'black') +
    facet_grid(taxonomic_metric ~ spectral_metric, scales = "free", labeller = labeller(taxonomic_metric = taxonomic_labels)) +
    labs(
      title = "Masked Rasters",
      color = 'Site') +
    theme_minimal() +
    scale_color_manual(values = c('darkgreen', 'chocolate', 'navajowhite2', 'darkseagreen')) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text.y = element_text(angle = 0, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank())
  
  # For 2016
  masked_16_plot <- df_long %>%
    filter(image_type == 'masked', year == 2016) %>%
    ggplot(aes(x = spectral_value, y = taxonomic_value, color = site)) +
    geom_point(show.legend = F, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = 'black') +
    facet_grid(taxonomic_metric ~ spectral_metric, scales = "free", labeller = labeller(taxonomic_metric = taxonomic_labels)) +
    theme_minimal() +
    scale_color_manual(values = c('darkgreen', 'darkseagreen')) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text.y = element_text(angle = 0, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank())
  
  unmasked_16_plot <- df_long %>%
    filter(image_type == 'unmasked', year == 2016) %>%
    ggplot(aes(x = spectral_value, y = taxonomic_value, color = site)) +
    geom_point(show.legend = F, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = 'black') +
    facet_grid(taxonomic_metric ~ spectral_metric, scales = "free", labeller = label_value) +
    labs(y = '2016') +
    theme_minimal() +
    scale_color_manual(values = c('darkgreen', 'darkseagreen')) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text.y.right = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # Combine plots for 2024
  combined_24_plot <- (unmasked_plot_24 + masked_plot_24) +
    plot_layout(ncol = 2, guides = "collect") & 
    theme(
      legend.position = "bottom",
      axis.title.x = element_blank()
    )
  
  # Combine plots for 2016
  combined_16_plot <- (unmasked_16_plot + masked_16_plot) +
    plot_layout(ncol = 2, guides = "collect") & 
    theme(
      legend.position = "bottom",
      axis.title.x = element_blank()
    )
  
  # Combine all plots for 2024 and 2016
  combined_all_plots <- (combined_24_plot / combined_16_plot) +
    plot_layout(ncol = 1) &
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 90, hjust = 1),
      strip.text.y = element_text(angle = 0, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  return(combined_all_plots)
}

create_plots(tax_and_spec_div_values)
create_plots(alive_tax_and_spec_div_values)
create_plots(rs_tax_and_spec_div_values)


# reflectance values across sites - subsample -----------------------------

raster_files_unmasked_24 <- list.files('data_out/combined_rasters/2024', pattern = '_combined_image.tif$', full.names = TRUE)
subplot_files <- list.files('data/fishnets', pattern = '_fishnet.shp$', full.names = TRUE)

unmasked_pixel_values_24 <- extract_pixel_values(raster_files_unmasked_24, subplot_files, c('blue', 'green', 'red', 'red_edge', 'nir'))

unmasked_24_sample <- unmasked_pixel_values_24 %>%
  filter(subplot_id == '1_1')

long_24_pixels_um <- melt(unmasked_24_sample, 
                id.vars = 'identifier', 
                measure.vars = c('red', 'blue', 'green', 'red_edge', 'nir'),
                variable.name = 'Wavelength', 
                value.name = 'PixelValue')

ggplot(long_24_pixels_um, aes(x = Wavelength, y = PixelValue, fill = identifier)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  labs(x = 'Wavelength',
       y = 'Standardised Pixel Value',
       fill = 'Site') +
  theme_minimal() +
  scale_x_discrete(labels= c('Red', 'Blue', 'Green', 'Red Edge', 'Near Infrared')) + 
  scale_fill_manual(values = c('darkgreen','chocolate','beige','darkseagreen')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

head(unmasked_24_sample)

ave_wavelengths_sample <- unmasked_24_sample %>%
  group_by(identifier) %>%
  summarise(across(blue:nir, mean, na.rm = TRUE)) 

av_wavelength_sample <- melt(ave_wavelengths_sample ,
                             id.varse = 'identifier',
                             measure.vars = c('red', 'blue', 'green', 'red_edge', 'nir'),
                             value.name = 'PixelValue',
                             variable.name = 'Wavelength')

av_wavelength_sample$Wavelength <- recode(av_wavelength_sample$Wavelength,
                                          blue = 450,
                                          green = 560,
                                          red = 650,
                                          red_edge = 730,
                                          nir = 840
                                          )
ggplot(av_wavelength_sample, aes(x = Wavelength, y = PixelValue, color = identifier)) +
  geom_line() +
  theme_minimal()


# presence/absense vs ndvi ------------------------------------------------
source('funx.R')
ndvi_files <- list.files(path = 'data/ndvi', pattern = '\\.tif$', full.names = T)
subplot_files <- list.files('data/fishnets', pattern = '_fishnet.shp$', full.names = TRUE)

ndvi_values <- extract_pixel_values(ndvi_files, subplot_files, 'ndvi')

mean_ndvi_values <- ndvi_values %>%
  group_by(identifier, subplot_id) %>%
  summarise(mean_ndvi = mean(ndvi, na.rm = T))

presence_proportion <- do.call(rbind, lapply(names(twentyfour_field_diversity$presence_absence_matrices), 
                                             function(site) {
  site_data <- twentyfour_field_diversity$presence_absence_matrices[[site]] %>%
    dplyr::select(subplot_id, presence_proportion) %>%
    mutate(identifier = site) 
  return(site_data)
}))

abundance_values <- left_join(mean_ndvi_values, presence_proportion, by = c('identifier', 'subplot_id'))
abundance_values <- abundance_values %>%
  mutate(mean_ndvi_normalized = (mean_ndvi - min(mean_ndvi, na.rm = TRUE)) / 
           (max(mean_ndvi, na.rm = TRUE) - min(mean_ndvi, na.rm = TRUE)))

ggplot(abundance_values, aes(x = presence_proportion, y = mean_ndvi_normalized, color = identifier)) +
  geom_point() +
  geom_smooth() +
  theme_minimal()

# ndvi and diversity ------------------------------------------------------

diversity_24 <- as.data.frame(twentyfour_field_diversity$final_results)

diversity_24 <- left_join(diversity_24, abundance_values, by = c('site' = 'identifier','subplot_id'))

ggplot(diversity_24, aes(x = mean_ndvi_normalized, y = shannon_diversity)) +
  geom_point() +
  geom_smooth() +
  theme_minimal()

