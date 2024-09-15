library(tidyverse)
library(glmmTMB)
library(easystats)
library(performance)
library(lme4)

library(tidyverse)
library(glmmTMB)
library(performance)


df <- tax_and_spec_diversity_values_24 %>%
  filter(image_type == 'masked')

model <- glmmTMB(species_richness ~ scale(CV_rf) + (1|site), data = df)
summary(model)

tax_and_spec_div_values <- read_csv('data_out/tax_and_spec_diversity_values.csv')

taxonomic_metrics <- c("species_richness", 
                       "shannon_diversity", 
                       "simpson_diversity", 
                       "pielou_evenness")
spectral_metrics <- c("CV_rf",
                      "SV", 
                      "CHV")
image_types <- c("masked", 
                 "unmasked")
years <- c("2016", 
           "2024")

# Initialize the results data frame
results_df <- data.frame(
  year = integer(),
  image_type = character(),
  taxonomic_metric = character(),
  spectral_metric = character(),
  r2_marginal = numeric(),
  r2_conditional = numeric(),
  correlation_coefficient = numeric(),
  assumption_check = character(),
  stringsAsFactors = FALSE
)


plot_dir <- "model_assumptions_plots" # Directory to save the plots
dir.create(plot_dir, showWarnings = FALSE) # Create directory if it doesn't exist

for (year in years) {
  for (img_type in image_types) {
    # filter for year and image type
    df_filtered <- tax_and_spec_div_values %>%
      filter(year == !!year, image_type == img_type)
    
    for (tax in taxonomic_metrics) {
      for (spec in spectral_metrics) {
        
        # fit model
        model <- try(glmmTMB(
          formula = as.formula(paste(tax, "~ scale(", spec, ") + (1 | site)")),
          data = df_filtered
        ), silent = TRUE)
        
        # handle any errors in fitting the model
        if (!inherits(model, "try-error")) {
          # Calculate R2 values using r2_nakagawa
          r2_values <- r2_nakagawa(model)
          
          r2_marginal <- format(r2_values$R2_marginal, scientific = FALSE)
          r2_conditional <- format(r2_values$R2_conditional, scientific = FALSE)
          
          # calculate predicted values from the model
          predicted_values <- predict(model)
          
          # calculate Pearson correlation between observed and predicted values
          correlation_coefficient <- cor(df_filtered[[tax]], predicted_values, use = "complete.obs")
          
          # check model assumptions
          assumption_check <- tryCatch({
            # use check_model to check all assumptions
            check_results <- check_model(model)
            
            # visualize and save the check results
            plot_file_name <- paste0(plot_dir, "/", year, "_", img_type, "_", tax, "_vs_", spec, ".png")
            ggsave(plot_file_name, plot(check_results), width = 10, height = 8)
            
            # return a summary or a flag if assumptions are met
            if (all(summary(check_results)$OK)) "Assumptions OK" else "Issues Detected"
          }, error = function(e) {
            "Assumption Check Error"
          })
          
          # store the results
          results_df <- rbind(results_df, data.frame(
            year = year,
            image_type = img_type,
            taxonomic_metric = tax,
            spectral_metric = spec,
            r2_marginal = r2_marginal,
            r2_conditional = r2_conditional,
            correlation_coefficient = correlation_coefficient,
            assumption_check = assumption_check
          ))
        }
      }
    }
  }
}

results_df %>%
  View()

write.csv(results_df, 'data_out/modelling_results2.csv')

# alive -------------------------------------------------------------------

library(glmmTMB)
library(performance)
library(dplyr)
library(readr)

alive_tax_and_spec_diversity_values <- read_csv('data_out/alive_tax_and_spec_diversity_values.csv')

taxonomic_metrics <- c("species_richness", "shannon_diversity", "simpson_diversity", 
                       "pielou_evenness", "exp_shannon", "inv_simpson")
spectral_metrics <- c("CV", "SV", "CHV", "CHV_nopca")
image_types <- c("masked", "unmasked")
years <- c('2016', '2024')


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

# loop through years, tax and spec metrics
for (year in years) {
  for (img_type in image_types) {
    # filter data for year and image type
    df_filtered <- alive_tax_and_spec_diversity_values %>%
      filter(year == !!year, image_type == img_type)
    
    for (tax in taxonomic_metrics) {
      for (spec in spectral_metrics) {
        
        # fit model
        model <- try(glmmTMB(
          formula = as.formula(paste(tax, "~ scale(", spec, ") + (1 | site)")),
          data = df_filtered
        ), silent = TRUE)
        
        # handle any errors in fitting the model
        if (!inherits(model, "try-error")) {
          # calculate R2 values using r2_nakagawa
          r2_values <- r2_nakagawa(model)
          
          r2_marginal <- format(r2_values$R2_marginal, scientific = FALSE)
          r2_conditional <- format(r2_values$R2_conditional, scientific = FALSE)
          
          # calculate Pearson correlation between observed and predicted values
          predicted_values <- predict(model)
          correlation_coefficient <- cor(df_filtered[[tax]], predicted_values, use = "complete.obs")
          
          # store results
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


print(alive_results_df)

write.csv(alive_results_df, 'data_out/alive_modelling_results.csv')


# resampled data ----------------------------------------------------------

#okay so first we will do year on year
#then combine and do one model for all years but use 'year' as a random effect??



rs_tax_and_spec_div_values <- read_csv('data_out/resampled24_and_16_spec_and_div_values.csv')

taxonomic_metrics <- c("species_richness", "shannon_diversity", "simpson_diversity", 
                       "pielou_evenness", "exp_shannon", "inv_simpson")
spectral_metrics <- c("CV", "SV", "CHV")
image_types <- c("masked", "unmasked")
years <- c('2016', '2024')


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


for (year in years) {
  for (img_type in image_types) {
    # filter data for year and image type
    df_filtered <- rs_tax_and_spec_div_values_1624 %>%
      filter(year == !!year, image_type == img_type)
    
    
    for (tax in taxonomic_metrics) {
      for (spec in spectral_metrics) {
        
        # fit model
        model <- try(glmmTMB(
          formula = as.formula(paste(tax, "~ scale(", spec, ") + (1 | site)")),
          data = df_filtered
        ), silent = TRUE)
        
        # handle any errors in fitting the model
        if (!inherits(model, "try-error")) {
          # calculate R2 values using r2_nakagawa
          r2_values <- r2_nakagawa(model)
          
          r2_marginal <- format(r2_values$R2_marginal, scientific = FALSE)
          r2_conditional <- format(r2_values$R2_conditional, scientific = FALSE)
          
          # calculate Pearson correlation between observed and predicted values
          predicted_values <- predict(model)
          correlation_coefficient <- cor(df_filtered[[tax]], predicted_values, use = "complete.obs")
          
          # store the results
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

write.csv(rs_results, 'data_out/rs_modelling_results.csv')

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
        
        
        # Store the results
        rs_results_2 <- rbind(rs_results_2, data.frame(
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



# tall plant analysis -----------------------------------------------------


tall_tax_and_spec_div_values <- read_csv('data_out/tall_tax_and_spec_diversity_values.csv')

taxonomic_metrics <- c("species_richness", "shannon_diversity", "simpson_diversity", 
                       "pielou_evenness", "exp_shannon", "inv_simpson")
spectral_metrics <- c("CV", "CV_rf","SV", "CHV", "CHV_nopca")
image_types <- c("masked", "unmasked")
years <- c('2016', '2024')


tall_results_df <- data.frame(
  year = integer(),
  image_type = character(),
  taxonomic_metric = character(),
  spectral_metric = character(),
  r2_marginal = numeric(),
  r2_conditional = numeric(),
  correlation_coefficient = numeric(),
  stringsAsFactors = FALSE
)

# loop through years, tax and spec metrics
for (year in years) {
  for (img_type in image_types) {
    # filter data for year and image type
    df_filtered <- tall_tax_and_spec_div_values %>%
      filter(year == !!year, image_type == img_type)
  
    for (tax in taxonomic_metrics) {
      for (spec in spectral_metrics) {
        
        # fit model
        model <- try(glmmTMB(
          formula = as.formula(paste(tax, "~ scale(", spec, ") + (1 | site)")),
          data = df_filtered
        ), silent = TRUE)
        
        # handle any errors in fitting the model
        if (!inherits(model, "try-error")) {
          # calculate R2 values using r2_nakagawa
          r2_values <- r2_nakagawa(model)
          
          r2_marginal <- format(r2_values$R2_marginal, scientific = FALSE)
          r2_conditional <- format(r2_values$R2_conditional, scientific = FALSE)
          
          # calculate Pearson correlation between observed and predicted values
          predicted_values <- predict(model)
          correlation_coefficient <- cor(df_filtered[[tax]], predicted_values, use = "complete.obs")
          
          # store the results
          tall_results_df <- rbind(tall_results_df, data.frame(
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
print(tall_results_df)

write.csv(alive_results_df, 'data_out/alive_modelling_results.csv')




