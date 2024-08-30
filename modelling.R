library(tidyverse)
library(glmmTMB)
library(easystats)
library(performance)
library(lme4)


tax_and_spec_diversity_values <- read_csv('data_out/tax_and_spec_diversity_values.csv')

taxonomic_metrics <- c("species_richness", "shannon_diversity", "simpson_diversity", "pielou_evenness")
spectral_metrics <- c("CV", "SV", "CHV", "CHV_nopca")
image_types <- c("masked", "unmasked")
years <- c('2016', '2024')

# Initialize the results data frame
results_df <- data.frame(
  year = integer(),
  image_type = character(),
  taxonomic_metric = character(),
  spectral_metric = character(),
  r2_marginal = numeric(),
  r2_conditional = numeric(),
  correlation_coefficient = numeric(),
  stringsAsFactors = FALSE
)

for (year in years) {
  for (img_type in image_types) {
    # Filter data for the current year and image type
    df_filtered <- tax_and_spec_diversity_values %>%
      filter(year == !!year, image_type == img_type)
    
    for (tax in taxonomic_metrics) {
      for (spec in spectral_metrics) {
        
        # Fit the model
        model <- try(glmmTMB(
          formula = as.formula(paste(tax, "~ scale(", spec, ") + (1 | site)")),
          data = df_filtered
        ), silent = TRUE)
        
        # Handle any errors in fitting the model
        if (!inherits(model, "try-error")) {
          # Calculate R2 values using r2_nakagawa
          r2_values <- r2_nakagawa(model)
          
          r2_marginal <- format(r2_values$R2_marginal, scientific = FALSE)
          r2_conditional <- format(r2_values$R2_conditional, scientific = FALSE)
          
          # Calculate predicted values from the model
          predicted_values <- predict(model)
          
          # Calculate Pearson correlation between observed and predicted values
          correlation_coefficient <- cor(df_filtered[[tax]], predicted_values, use = "complete.obs")
          
          # Store the results
          results_df <- rbind(results_df, data.frame(
            year = year,
            image_type = img_type,
            taxonomic_metric = tax,
            spectral_metric = spec,
            r2_marginal = r2_marginal,
            r2_conditional = r2_conditional,
            correlation_coefficient = correlation_coefficient
          ))
        }
      }
    }
  }
}



# Print the final results
print(results_df)




# alive -------------------------------------------------------------------


alive_tax_and_spec_diversity_values <- read_csv('data_out/alive_tax_and_spec_diversity_values.csv')
library(glmmTMB)
library(performance)
library(dplyr)
library(readr)

# Load the data
alive_tax_and_spec_diversity_values <- read_csv('data_out/alive_tax_and_spec_diversity_values.csv')

# Define metrics, image types, and years
taxonomic_metrics <- c("species_richness", "shannon_diversity", "simpson_diversity", "pielou_evenness")
spectral_metrics <- c("CV", "SV", "CHV", "CHV_nopca")
image_types <- c("masked", "unmasked")
years <- c('2016', '2024')

# Initialize the results data frame
alive_results_df <- data.frame(
  year = integer(),
  image_type = character(),
  taxonomic_metric = character(),
  spectral_metric = character(),
  r2_marginal = numeric(),
  r2_conditional = numeric(),
  correlation_coefficient = numeric(),
  stringsAsFactors = FALSE
)

# Loop through years, image types, and metrics
for (year in years) {
  for (img_type in image_types) {
    # Filter data for the current year and image type
    df_filtered <- alive_tax_and_spec_diversity_values %>%
      filter(year == !!year, image_type == img_type)
    
    # Check if the filtered data is empty
    if (nrow(df_filtered) == 0) {
      cat("No data for year:", year, "and image type:", img_type, "\n")
      next
    }
    
    for (tax in taxonomic_metrics) {
      for (spec in spectral_metrics) {
        
        # Fit the model
        model <- try(glmmTMB(
          formula = as.formula(paste(tax, "~ scale(", spec, ") + (1 | site)")),
          data = df_filtered
        ), silent = TRUE)
        
        # Handle any errors in fitting the model
        if (!inherits(model, "try-error")) {
          # Calculate R2 values using r2_nakagawa
          r2_values <- r2_nakagawa(model)
          
          r2_marginal <- format(r2_values$R2_marginal, scientific = FALSE)
          r2_conditional <- format(r2_values$R2_conditional, scientific = FALSE)
          
          # Calculate Pearson correlation between observed and predicted values
          predicted_values <- predict(model)
          correlation_coefficient <- cor(df_filtered[[tax]], predicted_values, use = "complete.obs")
          
          # Store the results
          alive_results_df <- rbind(alive_results_df, data.frame(
            year = year,
            image_type = img_type,
            taxonomic_metric = tax,
            spectral_metric = spec,
            r2_marginal = r2_marginal,
            r2_conditional = r2_conditional,
            correlation_coefficient = correlation_coefficient
          ))
        } else {
          cat("Error fitting model for year:", year, ", image type:", img_type, 
              ", tax metric:", tax, ", spectral metric:", spec, "\n")
        }
      }
    }
  }
}

# Check the results
print(alive_results_df)


# resampled data ----------------------------------------------------------

#okay so first we will do year on year
#then combine and do one model for all years but use 'year' as a random effect??


rs_tax_and_spec_div_values_1624 <- rbind(rs_tax_and_spec_div_values, tax_and_spec_diversity_values_16)

taxonomic_metrics <- c("species_richness", "shannon_diversity", "simpson_diversity", "pielou_evenness")
spectral_metrics <- c("CV", "SV", "CHV")
image_types <- c("masked", "unmasked")
years <- c('2016', '2024')

# Initialize the results data frame
rs_results <- data.frame(
  year = integer(),
  image_type = character(),
  taxonomic_metric = character(),
  spectral_metric = character(),
  r2_marginal = numeric(),
  r2_conditional = numeric(),
  correlation_coefficient = numeric(),
  stringsAsFactors = FALSE
)

# Loop through years, image types, and metrics
for (year in years) {
  for (img_type in image_types) {
    # Filter data for the current year and image type
    df_filtered <- rs_tax_and_spec_div_values_1624 %>%
      filter(year == !!year, image_type == img_type)
    
    # Check if the filtered data is empty
    if (nrow(df_filtered) == 0) {
      cat("No data for year:", year, "and image type:", img_type, "\n")
      next
    }
    
    for (tax in taxonomic_metrics) {
      for (spec in spectral_metrics) {
        
        # Fit the model
        model <- try(glmmTMB(
          formula = as.formula(paste(tax, "~ scale(", spec, ") + (1 | site)")),
          data = df_filtered
        ), silent = TRUE)
        
        # Handle any errors in fitting the model
        if (!inherits(model, "try-error")) {
          # Calculate R2 values using r2_nakagawa
          r2_values <- r2_nakagawa(model)
          
          r2_marginal <- format(r2_values$R2_marginal, scientific = FALSE)
          r2_conditional <- format(r2_values$R2_conditional, scientific = FALSE)
          
          # Calculate Pearson correlation between observed and predicted values
          predicted_values <- predict(model)
          correlation_coefficient <- cor(df_filtered[[tax]], predicted_values, use = "complete.obs")
          
          # Store the results
          rs_results <- rbind(rs_results, data.frame(
            year = year,
            image_type = img_type,
            taxonomic_metric = tax,
            spectral_metric = spec,
            r2_marginal = r2_marginal,
            r2_conditional = r2_conditional,
            correlation_coefficient = correlation_coefficient
          ))
        } else {
          cat("Error fitting model for year:", year, ", image type:", img_type, 
              ", tax metric:", tax, ", spectral metric:", spec, "\n")
        }
      }
    }
  }
}


# add year as random effect to your model and combine. would year or subplot
# be the random effect?

taxonomic_metrics <- c("species_richness", "shannon_diversity", "simpson_diversity", "pielou_evenness")
spectral_metrics <- c("CV", "SV", "CHV")
image_types <- c("masked", "unmasked")
years <- c('2016', '2024')


rs_results_2 <- data.frame(
  year = integer(),
  image_type = character(),
  taxonomic_metric = character(),
  spectral_metric = character(),
  r2_marginal = numeric(),
  r2_conditional = numeric(),
  correlation_coefficient = numeric(),
  stringsAsFactors = FALSE
)

# Update model loop to include year as a fixed effect and subplots nested within site as random
for (img_type in image_types) {
  # Filter combined data for the current image type
  df_filtered <- rs_tax_and_spec_div_values_1624 %>%
    filter(image_type == img_type)
  
  for (tax in taxonomic_metrics) {
    for (spec in spectral_metrics) {
      model <- try(glmmTMB(
        formula = as.formula(paste(tax, "~ scale(", spec, ") + year + (1 | site/subplot_id)")),
        data = df_filtered
      ), silent = TRUE)
      
      if (!inherits(model, "try-error")) {
        r2_values <- r2_nakagawa(model)
        
        # Calculate Pearson correlation between observed and predicted values
        predicted_values <- predict(model)
        correlation_coefficient <- cor(df_filtered[[tax]], predicted_values, use = "complete.obs")
        
        # Store the results
        rs_results_2 <- rbind(rs_results, data.frame(
          year = "2016 & 2024",
          image_type = img_type,
          taxonomic_metric = tax,
          spectral_metric = spec,
          r2_marginal = r2_values$R2_marginal,
          r2_conditional = r2_values$R2_conditional,
          correlation_coefficient = correlation_coefficient
        ))
      } else {
        cat("Error fitting model for image type:", img_type, ", tax metric:", tax, 
            ", spectral metric:", spec, "\n")
      }
    }
  }
}






