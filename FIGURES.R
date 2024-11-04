# FIGURES, BABY! 
setwd("C:/Users/adele/Documents/fg_spectral_diversity")

library(reshape2)
library(ggplot2)
source('funx.R')
library(tools)
library(tidyverse)
library(sf)
library(patchwork)
install.packages('ggh4x')
library(ggh4x)
# TRY TIDY MODELS FOR YOUR MODEL!

# main figure ------------------------------------------------------------


tax_and_spec_div_values <- read.csv('data_out/tax_and_spec_diversity_values.csv')


df_long <- tax_and_spec_div_values %>%
  pivot_longer(cols = c(species_richness, exp_shannon, inv_simpson, pielou_evenness),
               names_to = "taxonomic_metric",
               values_to = "taxonomic_value") %>%
  pivot_longer(cols = c(CV, SV, log.CHV_nopca),
               names_to = "spectral_metric",
               values_to = "spectral_value") %>%
  mutate(taxonomic_metric = factor(taxonomic_metric, levels = names(taxonomic_labels)))

df_long$predicted_values <- NA

for (year in unique(tax_and_spec_div_values$year)) {
  for (img_type in unique(tax_and_spec_div_values$image_type)) {
    # Filter for year and image type from the original data
    df_filtered <- tax_and_spec_div_values %>%
      filter(year == !!year, image_type == img_type)
    
    for (tax in unique(df_long$taxonomic_metric)) {
      for (spec in unique(df_long$spectral_metric)) {
        
        # Fit the model for the specific tax and spectral metrics
        model <- try(glmmTMB(
          formula = as.formula(paste(tax, "~ scale(", spec, ") + (1 | site)")),
          data = df_filtered
        ), silent = TRUE)
        
        # Check if the model fitting was successful
        if (!inherits(model, "try-error")) {
          # Predict values using the model
          predictions <- predict(model, newdata = df_filtered, type = "response")
          
          # Add predictions to the relevant rows in df_long
          df_long$predicted_values[df_long$year == year & df_long$image_type == img_type & 
                                     df_long$taxonomic_metric == tax & 
                                     df_long$spectral_metric == spec & 
                                     df_long$subplot_id %in% df_filtered$subplot_id & 
                                     df_long$site %in% df_filtered$site] <- predictions
        }
      }
    }
  }
}


create_plots <- function(df) {
  # taxonomic labels in correct order
  taxonomic_labels <- c(
    species_richness = "Species\nRichness",
    exp_shannon = "Exponential\nShannon's",
    inv_simpson = "Inverse\nSimpson's",
    pielou_evenness = "Pielou's\nEvenness"
  )
  
  spectral_labels <- c(
    CV = 'CV',
    log.CHV_nopca = 'log(CHV)',
    SV = 'SV'
  )
  
  # pivot longer
  #df_long <- df %>%
  #  pivot_longer(cols = c(species_richness, exp_shannon, inv_simpson, pielou_evenness),
  #               names_to = "taxonomic_metric",
  #               values_to = "taxonomic_value") %>%
  #  pivot_longer(cols = c(CV, SV, log.CHV_nopca),
  #               names_to = "spectral_metric",
  #               values_to = "spectral_value") %>%
  #  mutate(taxonomic_metric = factor(taxonomic_metric, levels = names(taxonomic_labels)))
  
  spectral_range <- range(df_long$spectral_value, na.rm = TRUE)
  taxonomic_range <- range(df_long$taxonomic_value, na.rm = T)

  results_df <- results_df %>% 
    mutate(taxonomic_metric = factor(taxonomic_metric, levels = names(taxonomic_labels)))
  
  # 2024 plots
  unmasked_plot_24 <- df_long %>%
    filter(image_type == 'unmasked', year == 2024) %>%
    ggplot(aes(x = spectral_value, y = taxonomic_value, color = site)) +
    geom_point(alpha = 0.7, show.legend = F) +
    #geom_abline(aes(intercept = intercept, slope = beta), 
    #            data = results_df %>% filter(year == 2024, image_type == 'unmasked'), 
    #            color = "black", linewidth = 1) +
    geom_line(aes(x = spectral_value, y = predicted_values), linetype = "solid", color = "black") + 
    #geom_smooth(method = 'lm', se = FALSE, linetype = "solid", color = 'black') +
    facet_grid(taxonomic_metric ~ spectral_metric, scales = "free", labeller = labeller(spectral_metric = spectral_labels)) +
    labs(
      title = "Unmasked Rasters",
      y = '2024') +
    theme_minimal() +
    scale_color_manual(values = c('darkgreen', 'chocolate', 'navajowhite2', 'darkseagreen')) +
    ggh4x::facetted_pos_scales(y = list(
      taxonomic_metric == 'species_richness' ~ scale_y_continuous(breaks = c(5, 10, 15), limits = c(0, 20)),
      taxonomic_metric == 'exp_shannon' ~ scale_y_continuous(breaks = c(5, 10, 15), limits = c(0.7, 16)),
      taxonomic_metric == 'inv_simpson' ~ scale_y_continuous(breaks = c(4, 8, 12), limits = c(0, 12.5)),
      taxonomic_metric == 'pielou_evenness' ~ scale_y_continuous(breaks = c(0.8, 0.9, 1), limits = c(0.65, 1)))) +
    ggh4x::facetted_pos_scales(x = list(
      spectral_metric == 'CV' ~ scale_x_continuous(breaks = c(0.2, 0.3, 0.4, 0.5))
    )) +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.y.right = element_blank())
  
  masked_plot_24 <- df_long %>%
    filter(image_type == 'masked', year == 2024) %>%
    ggplot(aes(x = spectral_value, y = taxonomic_value, color = site)) +
    geom_point(alpha = 0.7) +
    #geom_line(aes(y = predicted_values), linetype = "solid", color = "black") + 
    #geom_smooth(method = 'lm', se = FALSE, linetype = "solid", color = 'black') +
    facet_grid(taxonomic_metric ~ spectral_metric, scales = "free", labeller = labeller(taxonomic_metric = taxonomic_labels,
                                                                                        spectral_metric = spectral_labels)) +
    labs(
      title = "Masked Rasters",
      color = 'Site') +
    theme_minimal() +
    scale_color_manual(values = c('darkgreen', 'chocolate', 'navajowhite2', 'darkseagreen')) +
    ggh4x::facetted_pos_scales(y = list(
      taxonomic_metric == 'species_richness' ~ scale_y_continuous(breaks = c(5, 10, 15), limits = c(0, 20)),
      taxonomic_metric == 'exp_shannon' ~ scale_y_continuous(breaks = c(5, 10, 15), limits = c(0.7, 16)),
      taxonomic_metric == 'inv_simpson' ~ scale_y_continuous(breaks = c(4, 8, 12), limits = c(0, 12.5)),
      taxonomic_metric == 'pielou_evenness' ~ scale_y_continuous(breaks = c(0.8, 0.9, 1), limits = c(0.65, 1)))) +    
    ggh4x::facetted_pos_scales(x = list(
      spectral_metric == 'CV' ~ scale_x_continuous(breaks = c(0.2, 0.3, 0.4, 0.5))
    )) +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
          strip.text.y = element_text(angle = 0, hjust = 1),
          strip.text.y.right = element_text(lineheight = 0.8, size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank())
  
  # For 2016
  masked_16_plot <- df_long %>%
    filter(image_type == 'masked', year == 2016) %>%
    ggplot(aes(x = spectral_value, y = taxonomic_value, color = site)) +
    geom_point(show.legend = F, alpha = 0.7) +
    #geom_line(aes(y = predicted_values), linetype = "solid", color = "black") + 
    #geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = 'black') +
    facet_grid(taxonomic_metric ~ spectral_metric, scales = "free", labeller = labeller(taxonomic_metric = taxonomic_labels,
                                                                                        spectral_metric = spectral_labels)) +
    theme_minimal() +
    scale_color_manual(values = c('darkgreen', 'darkseagreen')) +
    ggh4x::facetted_pos_scales(y = list(
      taxonomic_metric == 'species_richness' ~ scale_y_continuous(breaks = c(5, 10, 15), limits = c(0, 20)),
      taxonomic_metric == 'exp_shannon' ~ scale_y_continuous(breaks = c(5, 10, 15), limits = c(0.7, 16)),
      taxonomic_metric == 'inv_simpson' ~ scale_y_continuous(breaks = c(4, 8, 12), limits = c(0, 12.5)),
      taxonomic_metric == 'pielou_evenness' ~ scale_y_continuous(breaks = c(0.8, 0.9, 1), limits = c(0.7, 1)))) +    
    ggh4x::facetted_pos_scales(x = list(
      spectral_metric == 'CV' ~ scale_x_continuous(breaks = c(0.2, 0.3, 0.4, 0.5))
    )) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
          strip.text.y = element_text(angle = 0, hjust = 1),
          strip.text.y.right = element_text(lineheight = 0.8, size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank())
  
  unmasked_16_plot <- df_long %>%
    filter(image_type == 'unmasked', year == 2016) %>%
    ggplot(aes(x = spectral_value, y = taxonomic_value, color = site)) +
    geom_point(show.legend = F, alpha = 0.7) +
    #geom_line(aes(y = predicted_values), linetype = "solid", color = "black") + 
    #geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = 'black') +
    facet_grid(taxonomic_metric ~ spectral_metric, scales = "free", labeller = labeller(spectral_metric = spectral_labels)) +
    labs(y = '2016') +
    theme_minimal() +
    scale_color_manual(values = c('darkgreen', 'darkseagreen')) +
    ggh4x::facetted_pos_scales(y = list(
      taxonomic_metric == 'species_richness' ~ scale_y_continuous(breaks = c(5, 10, 15), limits = c(0, 20)),
      taxonomic_metric == 'exp_shannon' ~ scale_y_continuous(breaks = c(5, 10, 15), limits = c(0.7, 16)),
      taxonomic_metric == 'inv_simpson' ~ scale_y_continuous(breaks = c(4, 8, 12), limits = c(0, 12.5)),
      taxonomic_metric == 'pielou_evenness' ~ scale_y_continuous(breaks = c(0.8, 0.9, 1), limits = c(0.7, 1)))) +    
    ggh4x::facetted_pos_scales(x = list(
      spectral_metric == 'CV' ~ scale_x_continuous(breaks = c(0.2, 0.3, 0.4, 0.5))
    )) +
theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
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
      axis.text.x = element_text(angle = 90, hjust = 1, , size = 10),
      axis.text.y = element_text(size = 10),
      strip.text.y = element_text(angle = 0, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  return(combined_all_plots)
}

create_plots(tax_and_spec_div_values)
create_plots(tall_tax_and_spec_diversity_values_24)
create_plots(alive_tax_and_spec_div_values)
create_plots(rs_tax_and_spec_div_values)





ggsave('maps_graphs/main_fig.png', height = 8, width = 10, dpi = 600)

library(terra)
sandstone_masked <- rast('data_out/combined_rasters/masked/2024/NSABHC0011_combined_image_masked.tif')
plot(sandstone_masked[[1]])
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
       y = 'Reflectance',
       fill = 'Site') +
  theme_minimal() +
  scale_x_discrete(labels= c('Red', 'Blue', 'Green', 'Red Edge', 'Near Infrared')) + 
  scale_fill_manual(values = c('darkgreen','chocolate','beige','darkseagreen')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 16))

ggsave('maps_graphs/reflectance.png', width = 4, height = 3, dpi = 300)

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

#need to amend stack() in extract_pixel_values to raster() for this to work :)
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

# alive and normy r2 ------------------------------------------------------
library(readr)

results_df <- read_csv('data_out/modelling_results2.csv')
alive_results_df <- read_csv('data_out/alive_modelling_results.csv')
 
results_df$analysis <- 'all plants'
alive_results_df$analysis <- 'alive plants only'

merged_results$r2_marginal_results <- as.numeric(merged_results$r2_marginal_results)
merged_results$r2_marginal_alive <- as.numeric(merged_results$r2_marginal_alive)
merged_results$r2_marginal_tall <- as.numeric(merged_results$r2_marginal_tall)

tall_results_df$year <- as.numeric(tall_results_df$year)
results_df$year <- as.numeric(results_df$year)

tall_all_results <- results_df %>%
  inner_join(tall_results_df, by = c("year", "image_type", "spectral_metric", "taxonomic_metric"), suffix = c("_results", "_tall"))
tall_all_results$r2_marginal_results <- as.numeric(tall_all_results$r2_marginal_results)
tall_all_results$r2_marginal_tall <- as.numeric(tall_all_results$r2_marginal_tall)
library(showtext)
font_add_google('Roboto','roboto')
showtext_auto()

tall_all_results %>%
  distinct %>%
  filter(taxonomic_metric != 'pielou_evenness') %>%
 filter(year == 2024) %>%
ggplot(aes(x = r2_marginal_results, y = r2_marginal_tall)) +
  geom_point(aes(shape = taxonomic_metric, 
                 color = spectral_metric,
                 alpha = image_type), 
             size = 4,
             position = position_jitter(width = 0.005, height = 0.005)) +
  geom_abline(intercept = 0, slope = 1, color = "grey50", linetype = "dashed") +
  scale_x_continuous(limits = c(0,0.4),
                    breaks = seq(0, 0.4, by = 0.1)) +
  scale_y_continuous(limits = c(0,0.4),
                     breaks = seq(0, 0.4, by = 0.1)) +
  scale_color_manual(values = c('#4B8F29', '#E6AB02', '#D95F02')) +
  scale_shape_manual(values = c(15:17),
                     labels = c('Exponential Shannons',
                                'Inverse Simpsons',
                                'Species Richness')) +
  scale_alpha_manual(values = c(0.4, 1)) +
  labs(x = 'R2 (all plant model)', 
       y = 'R2 (tall plant model)',
       color = 'Spectral Metric',
       shape = 'Taxonomic Metric',
       alpha = 'Image Type') + 
  theme_minimal(base_family = 'roboto') +
  theme(panel.grid.minor = element_blank(),  # Remove minor grid lines for clarity
        panel.grid.major.x = element_line(color = "grey90", size = 0.5),  # Subtle major grid lines
        panel.grid.major.y = element_line(color = "grey90", size = 0.5),
        legend.position = "right",  # Adjust legend position
        legend.box = "vertical",  # Stack legend items vertically
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Center and bold title
        axis.title.x = element_text(face = "bold"),  # Bold axis titles
        axis.title.y = element_text(face = "bold"),
        text = element_text(size = 22),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))

ggsave('maps_graphs/tall_all_model.png', height = 4, width = 5, dpi = 300, bg = 'white')

alive_all_results <- results_df %>%
  inner_join(alive_results_df, by = c("year", "image_type", "spectral_metric", "taxonomic_metric"), suffix = c("_results", "_alive"))

alive_all_results %>%
filter(taxonomic_metric != 'pielou_evenness') %>%
  ggplot(aes(x = r2_marginal_results, y = r2_marginal_alive)) +
  geom_point(aes(shape = taxonomic_metric, color = spectral_metric), 
             size = 3,
             position = position_jitter(width = 0.005, height = 0.005)) +
  geom_abline(intercept = 0, slope = 1, color = "grey50", linetype = "dashed") +
  scale_x_continuous(limits = c(0,0.4),
                     breaks = seq(0, 0.4, by = 0.1)) +
  scale_y_continuous(limits = c(0,0.4),
                     breaks = seq(0, 0.4, by = 0.1)) +
  scale_color_manual(values = c('#4B8F29', '#E6AB02', '#D95F02')) +
  scale_shape_manual(values = c(15:17),
                     labels = c('Exponential Shannons',
                                'Inverse Simpsons',
                                'Species Richness')) +
  labs(x = 'R2 (all plant model)', 
       y = 'R2 (alive plant model)',
       color = 'Spectral Metric',
       shape = 'Taxonomic Metric') + 
  theme_minimal(base_family = 'roboto') +
  theme(panel.grid.minor = element_blank(),  # Remove minor grid lines for clarity
        panel.grid.major.x = element_line(color = "grey80", size = 0.5),  # Subtle major grid lines
        panel.grid.major.y = element_line(color = "grey80", size = 0.5),
        legend.position = "right",  # Adjust legend position
        legend.box = "vertical",  # Stack legend items vertically
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Center and bold title
        axis.title.x = element_text(face = "bold"),  # Bold axis titles
        axis.title.y = element_text(face = "bold"))



# comparing masked and unmasked -------------------------------------------
results_df$r2_marginal <- as.numeric(results_df$r2_marginal)
unmasked_df <- results_df %>%
  filter(image_type == "unmasked") %>%
  dplyr::select(year, taxonomic_metric, spectral_metric, r2_marginal) %>%
  rename(r2_unmasked = r2_marginal)

masked_df <- results_df %>%
  filter(image_type == "masked") %>%
  dplyr::select(year, taxonomic_metric, spectral_metric, r2_marginal) %>%
  rename(r2_masked = r2_marginal)

# Merge the two data frames on year, taxonomic_metric, and spectral_metric
merged_df <- left_join(unmasked_df, masked_df, by = c("year", "taxonomic_metric", "spectral_metric"))

merged_df %>% 
  filter(taxonomic_metric != 'pielou_evenness') %>%
ggplot(aes(x = r2_unmasked, y = r2_masked, shape = taxonomic_metric, color = spectral_metric)) +
  geom_point(aes(shape = taxonomic_metric, color = spectral_metric, alpha = year), 
             size = 4, 
             position = position_jitter(width = 0.005, height = 0.005)) +
  geom_abline(intercept = 0, slope = 1, color = "grey50", linetype = "dashed") +
  scale_x_continuous(limits = c(0,0.4),
                     breaks = seq(0, 0.4, by = 0.1)) +
  scale_y_continuous(limits = c(0,0.4),
                     breaks = seq(0, 0.4, by = 0.1)) +
  scale_color_manual(values = c('#4B8F29', '#E6AB02', '#D95F02')) +
  scale_shape_manual(values = c(15:17),
                     labels = c('Exponential Shannons',
                                'Inverse Simpsons',
                                'Species Richness')) +
  scale_alpha_manual(values = c(0.4, 1)) +
  labs(x = "R²  (Unmasked Images)",
       y = "R²  (Masked Images)",
       color = 'Spectral metric',
       shape = 'Taxonomic metric',
       alpha = 'Year') +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),  # Remove minor grid lines for clarity
         panel.grid.major.x = element_line(color = "grey80", size = 0.5),  # Subtle major grid lines
         panel.grid.major.y = element_line(color = "grey80", size = 0.5),
         legend.position = "right",  # Adjust legend position
         legend.box = "vertical",  # Stack legend items vertically
         plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Center and bold title
         axis.title.x = element_text(face = "bold"),  # Bold axis titles
         axis.title.y = element_text(face = "bold"))







library(ggplot2)

results_df <- results_df %>%
  mutate(year = factor(year, levels = c(2024, 2016))) 

results_df %>%
  mutate(taxonomic_metric = factor(taxonomic_metric, 
                                   levels = c('species_richness', 'exp_shannon', 'inv_simpson', 'pielou_evenness'))) %>%
  ggplot() +
  geom_errorbarh(aes(xmin = beta_ci_lower, xmax = beta_ci_upper, 
                     y = spectral_metric, color = spectral_metric), height = 0.1, size = 0.8) +
  geom_point(aes(x = beta, y = spectral_metric, shape = image_type, 
                 color = spectral_metric), size = 3.5, stroke = 1, fill = 'white') +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey50') +
  scale_color_manual(values = c('#4B8F29', '#E6AB02', '#D95F02')) +
  scale_shape_manual(values = 21:22) + 
  facet_grid(year ~ taxonomic_metric, scales = "free_x",
             labeller = as_labeller(c(
               species_richness = "Species Richness",
               inv_simpson = "Inverse Simpson",
               exp_shannon = "Exponential Shannon",
               pielou_evenness = "Pielou's Evenness",
               '2016' = '2016',
               '2024' = '2024'))) +
  scale_y_discrete(labels = c(
    CV_rf = "CV", 
    SV = "SV", 
    log.CHV_nopca = "log(CHV)"
  )) + 
  labs(y = '',
       x = 'Standardised Slope',
       shape = 'Image type') +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 100),
        legend.title = element_text(size = 98),
        legend.text = element_text(size = 96),
        strip.text = element_text(size = 98),
        strip.text.y.right = element_text(face = 'bold')) +
        guides(color = 'none') 


ggsave('maps_graphs/beta_metrics.png', height = 6, width= 14, dpi = 600, bg = 'white')

# Remove the image_type part since your columns don't have masked/unmasked suffixes
library(tidyr)

results_long <- results_df %>%
  dplyr::select(year, spectral_metric, taxonomic_metric, beta, beta_ci_lower, beta_ci_upper, image_type) %>%
  pivot_wider(
    names_from = image_type,
    values_from = c(beta, beta_ci_lower, beta_ci_upper),
    names_sep = "_"
  )

grep("^beta|beta_ci_lower|beta_ci_upper", names(results_df), value = TRUE)





masked_ss <- rast('data_out/combined_rasters/masked/2024/NSABHC0011_combined_image_masked.tif')
emu <- rast('data_out/combined_rasters/2024/NSABHC0009_combined_image.tif')

plotRGB(masked_ss, r = 3, g = 2, b = 1, stretch = "lin", main = "RGB Image")

plotRGB(emu, r = 3, g = 2, b = 1, stretch = "lin", main = "RGB Image")

png('maps_graphs/emu_RGB.png', res = 300)

red <- raster::subset(masked_emu, 3)  # Red band
green <- raster::subset(masked_emu, 2)  # Green band
blue <- raster::subset(masked_emu, 1)  # Blue band

rgb_stack <- stack(red, green, blue)


# % cover vs diversity ----------------------------------------------------

pixel_per_plot_24_masked <- masked_pixel_values_24 %>%
  na.omit() %>%
  group_by(identifier, subplot_id) %>%
  summarise(count = n()) 

pixel_per_plot_24_unmasked <- unmasked_pixel_values_24 %>%
  na.omit() %>%
  group_by(identifier, subplot_id) %>%
  summarise(count = n()) 

proportion_pixel <- data.frame(identifier = pixel_per_plot_24_masked$identifier,
                               subplot_id = pixel_per_plot_24_masked$subplot_id,
                               prop = pixel_per_plot_24_masked$count / pixel_per_plot_24_unmasked$count)

metrics_24 <- read_csv('data_out/tax_and_spec_diversity_values_24.csv')

proportion_pixel <- left_join(proportion_pixel, metrics_24 %>%
                                dplyr::select(site, subplot_id, species_richness,
                                       exp_shannon, inv_simpson, pielou_evenness),
                              by = c('identifier'='site', 'subplot_id'))

proportion_pixel2 <- comparison_df %>%
  dplyr::select(c('subplot_id','site','CV_rf_masked','CV_rf_unmasked', 'CHV_masked', 'CHV_unmasked', 'SV_unmasked','SV_masked')) %>%
  left_join(proportion_pixel, by = c('site'='identifier','subplot_id'))

prop_pix_long2 <- pivot_longer(proportion_pixel2, 
                              cols = c(CV_rf_masked, CV_rf_unmasked, CHV_masked, CHV_unmasked, SV_unmasked, SV_masked),
                              names_to = 'metric_masktype',
                              values_to = 'spectral_diversity')


prop_pix_long <- pivot_longer(proportion_pixel, 
                              cols = c(species_richness, exp_shannon, inv_simpson, pielou_evenness),
                              names_to = 'metric',
                              values_to = 'taxonomic_diversity')

prop_pix_long$metric <- factor(prop_pix_long$metric, 
                               levels = c('species_richness', 'exp_shannon', 'inv_simpson', 'pielou_evenness'))


ggplot(prop_pix_long, aes(x = prop, y = taxonomic_diversity, color = metric, group = metric)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_smooth(method ='lm', se = F, alpha = 0.6) +
  labs(x = "% vegetation cover", 
       y = "taxonomic diversity") +
  scale_color_manual(values = c('#4B8F29', '#E6AB02', '#D95F02', '#1F78B4'),
                     labels = c('species richness', 'exponential shannons', 'inverse simpsons', 'pielous evenness')) +
  labs(color = 'taxonomic metric') + 
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        legend.position = 'bottom')

ggplot(prop_pix_long2, aes(x = prop, y = spectral_diversity, color = metric_masktype)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_smooth(method ='lm', se = F, alpha = 0.6) +
  labs(x = "% vegetation cover") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        legend.position = 'bottom') 

# ausplots map and scatterplot diversity ----------------------------------
library(ausplotsR)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)
library(ggthemes)
ausplots <- get_ausplots(veg.vouchers = T, veg.PI = T)
veg <- ausplots$veg.vouch
veg.PI <- ausplots$veg.PI

richness <- veg %>%
  group_by(site_unique) %>%
  summarise(count = n())

richness_PI <- veg.PI %>%
  group_by(site_unique) %>%
  drop_na(herbarium_determination) %>%
  summarise(count = n_distinct(herbarium_determination))

locations <- ausplots$site.info %>%
  dplyr::select(site_unique, latitude, longitude) #%>%
  #st_as_sf(coords = c('longitude', 'latitude'), crs = 4326)
  
ausplot_richness <- left_join(richness, locations, by = 'site_unique')

fg_points <- ausplot_richness %>%
  subset(grepl("NSABHC0009|NSABHC0010|NSABHC0011|NSABHC0012", site_unique))

ggplot(ausplot_richness) +
  geom_point(aes(x = count, y = latitude), color = 'grey70') +
  geom_point(data = fg_points, aes(x = count, y = latitude), color = 'red', size = 3) + 
  labs(x = 'species richness') +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

aus <- ne_states(country = 'australia', returnclass = 'sf')

fg_points$fowlers_gap <- 'Fowlers Gap'

aus <- st_transform(aus, crs = 3577)

ausplot_richness <- ausplot_richness[order(ausplot_richness$count), ]

library(showtext)
font_add_google('Open Sans','open sans')
font_add_google('Tahoma','tahoma')
font_add_google('Roboto','roboto')
showtext_auto()

fg_points <- st_as_sf(fg_points, coords = c("longitude", "latitude"), crs = 4326)
fg_points <- st_transform(fg_points, crs = 3577)

ausplot_richness <- st_as_sf(ausplot_richness, coords = c("longitude", "latitude"), crs = 4326)
ausplot_richness <- st_transform(ausplot_richness, crs = 3577)


ggplot() +
  geom_sf(data = aus, fill = 'white', color = 'black') + 
  geom_sf(data = ausplot_richness, aes(color = count), alpha = 0.6, size = 3, height = 0.15, width = 0.15) +
  geom_sf(data = fg_points, aes(fill = fowlers_gap), color = 'red', shape = 0, size = 5, stroke = 1) + 
  scale_color_viridis_c(option = 'H') +
  labs(color = 'species richness') +  
  theme_minimal(base_family = 'roboto') +
  guides(fill = guide_legend(title = NULL)) +
  ylim(-5000000, -700000) +  # Adjust the limits for the projected system
  xlim(-2050000, 2200000) + 
  theme(text = element_text(size = 60),
        legend.title = element_text(margin = margin(b = 15)),
        panel.grid = element_blank()) +
coord_sf(crs = 3577)

ggsave('maps_graphs/map_of_ausplot_sites_with_richness.png', width = 10, height = 7, dpi = 300, bg = 'white')

name_stats <- twentyfour_survey_data %>%
  filter(!is.na(standardised_name)) %>%
  group_by(site_location_name, standardised_name) %>%
  summarise(
    count = n(),
    average_height = mean(height, na.rm = TRUE) 
  )


ggplot(name_stats) +
  geom_point(aes(x = count, y = average_height, color = site_location_name)) +
  ylim(NA, 0.1) +
  theme_minimal() +
  theme(panel.grid = element_blank())
