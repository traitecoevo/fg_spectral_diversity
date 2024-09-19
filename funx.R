

# COMBINE TIFS FUNCTION (FOR CREATING MULTI BAND IMAGE) -----------------------------------------------------------------------

create_multiband_image <- function(mosaics_dir, desired_band_order){
# folder list | recursive = won't pick folders within folders
folders <- list.dirs(mosaics_dir, full.names = FALSE, recursive = FALSE)
# need to make this part an argument e.g. option to exclude certain folders or include certain folders
folders <- folders[folders != "point_clouds"]

## NOTE: spectral band image tif file names must be named after their band (e.g., blue, nir, etc), 
#  otherwise change 'desired_band_order' to match file names
#  should be combined in wavelength order, esp for biodivmapR processes (i.e. as above)

# loop thru each folder 
for (folder in folders) {
  # create path
  folder_path <- file.path(mosaics_dir, folder)
  
  # list of tif files | \\. represents . (dots need to be escaped w \, \ need to be escaped with  \). $means at end of file name/string
  tif_files <- list.files(folder_path, pattern = "\\.tif$", full.names = TRUE)
  
  # load as raster
  rasters <- lapply(tif_files, terra::rast)
  
  # extract band names from file names
  band_names <- tools::file_path_sans_ext(basename(tif_files))
  
  # stack rasters and assign band names
  combined_image <- terra::rast(rasters)
  names(combined_image) <- band_names
  
  # reorder the bands based on the desired band order
  combined_image <- combined_image[[match(desired_band_order, band_names)]]
  
  #create output directory folder if it doesn't exist
  output_dir <- file.path("data_out/combined_rasters", substr(basename(mosaics_dir), 1, 4))
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # create output file as .tif and as .envi
  output_filename <- file.path(output_dir, paste0(folder, "_combined_image"))
  # write .tif file
  terra::writeRaster(combined_image, filename = paste0(output_filename, '.tif'), 
                     filetype = "GTiff", gdal = c("INTERLEAVE=BAND"), overwrite = TRUE)
  plot(combined_image)
}
}

# FIND OPTIMUM THRESHOLDS FUNCTION 
# class = class column for classifications (e.g. veg, ground etc) (col name)
# value = ndvi, nir - what threshold are you seeking (col name)
# site = plot reference (col name)

find_optimum_thresholds <- function(df, class, value, site, class_value) {
  # empty df
  threshold_df <- data.frame(site = character(), threshold = numeric())
  
  # iterate over sites
  for (site_name in unique(df[[site]])) {
    
    # filter current location
    site_data <- subset(df, df[[site]] == site_name)
    
    # binary outcome variable for veg and non-veg
    site_data$binary_class <- ifelse(site_data[[class]] == class_value, 1, 0)
    
    # ROC curve
    roc_result <- pROC::roc(site_data$binary_class, site_data[[value]])
    
    # find optimum threshold
    best_threshold <- coords(roc_result, 'best')$threshold
    
    # append to the result data frame
    threshold_df <- rbind(threshold_df, data.frame(site = site_name, threshold = best_threshold))
  }
  
  return(threshold_df)
}

# CREATE_MASKED_RASTER FUNCTION 
#input can be directory with a number of files, a single file, or string of files.
#ndvi and nir thresholds can be provided as a df, if there are diff optimum values per site
# or as a single value for all sites
# this function assumes that layers are stacked in WAVELENGTH ORDER
#think about how you can make this more general for users - e.g. it requires the plot id to
# be in the file name currently - think about usability
create_masked_raster <- function(input, output_dir, 
                                 band_wavelengths, 
                                 NDVI_Thresh = 0.2, NIR_Thresh = 0.2,
                                 Red_band = 660, NIR_band = 790,
                                 NDVI_Thresh_df = NULL, NIR_Thresh_df = NULL) {
  
  if (dir.exists(input)) {
    # list all ENVI or TIF files in the directory
    files <- list.files(file.path(input), pattern = '\\.(envi|tif)$', full.names = TRUE)
  } else if (file.exists(input) || is.character(input)) {
    # single file input or string of files
    files <- input
  } else {
    stop("Invalid input provided.")
  }
  
  print(paste("Files found:", files))
  
  if (length(files) == 0) {
    stop("No files found.")
  }
  
  for (file in files) {
    # Extract the site identifier from file name
    site_id <- strsplit(basename(file), "_")[[1]][1]
    
    # Check if NDVI_Thresh_df is provided and extract the relevant threshold values
    if (!is.null(NDVI_Thresh_df)) {
      if (site_id %in% NDVI_Thresh_df$site) {
        NDVI_Thresh <- NDVI_Thresh_df$threshold[NDVI_Thresh_df$site == site_id]
      } else {
        stop(paste("No NDVI threshold values found for site", site_id))
      }
    }
    
    # Check if NIR_Thresh_df is provided and extract the relevant threshold values
    if (!is.null(NIR_Thresh_df)) {
      if (site_id %in% NIR_Thresh_df$site) {
        NIR_Thresh <- NIR_Thresh_df$threshold[NIR_Thresh_df$site == site_id]
      } else {
        stop(paste("No NIR threshold values found for site", site_id))
      }
    }
    
    # Read the raster stack
    raster_data <- stack(file)
    
    # Identify the bands for Red and NIR
    red_layer <- which(band_wavelengths == Red_band)
    nir_layer <- which(band_wavelengths == NIR_band)
    
    if (length(red_layer) == 0 | length(nir_layer) == 0) {
      stop("Red or NIR band not found in the raster.")
    }
    
    # Extract Red and NIR bands
    red <- raster_data[[red_layer]]
    nir <- raster_data[[nir_layer]]
    
    # Calculate NDVI
    ndvi <- (nir - red) / (nir + red)
    
    # Create a mask based on NDVI and NIR thresholds
    mask <- (ndvi < NDVI_Thresh) | (nir < NIR_Thresh)
    
    # Apply the mask to the raster data
    raster_data_masked <- mask(raster_data, mask, maskvalue = TRUE, updatevalue = NA)
    
    # Save the masked raster
    masked_filename <- file.path(output_dir, paste0(file_path_sans_ext(basename(file)), '_masked.tif'))
    writeRaster(raster_data_masked, filename = masked_filename, format = "GTiff", overwrite = TRUE)
    
    print(paste("Masked raster saved to:", masked_filename))
    }
  }



## EXTRACT PIXEL VALUES FUNCTION

extract_pixel_values <- function(raster_files, subplot_files, wavelength_names){
  
  all_pixel_values_list <- list()
  
  for (raster_file in raster_files) {
    
    # identify the string that represents the site name
    identifier <- str_extract(basename(raster_file), "^[^_]+")
    
    #choose the corresponding subplot file
    subplot_file <- subplot_files[grep(paste0('^', identifier), basename(subplot_files))]
    
    # read in subplot file and select geometries
    subplots <- read_sf(subplot_file) %>% 
      dplyr::select('geometry')
    
    # apply subplot ids
    subplots$subplot_id <- unlist(lapply(1:5, function(i) paste(i, 1:5, sep="_")))
    
    # read in raster file
    raster_data <- stack(raster_file)
    
    # apply names - should be saved in wavelength order as per sect 1 of this script
    names(raster_data) <- wavelength_names
    
    # create empty list
    pixel_values_list <- list()
    
    for (i in 1:nrow(subplots)){
      
      # select the i-th subplot and its id
      subplot <- subplots[i, ]
      subplot_id <- subplot$subplot_id
      
      # convert to spatial object
      subplot_sp <- as(subplot, "Spatial")
      
      # crop and mask raster using current subplot 
      cropped_raster <- crop(raster_data, subplot_sp)
      masked_raster <- mask(cropped_raster, subplot_sp)
      
      # extract pixel values
      pixel_values  <- as.data.frame(getValues(masked_raster))
      
      # add subplot id to pixel values df
      pixel_values$subplot_id <- subplot_id
      
      #add to list
      pixel_values_list[[i]] <- pixel_values
      
    }
    # combined all pixel values into one df for current raster
    all_pixel_values <- bind_rows(pixel_values_list) %>%
      na.omit()
    
    # add to overall list with all raster data pixel values 
    all_pixel_values_list[[identifier]] <- all_pixel_values
  }
  combined_values <- bind_rows(all_pixel_values_list, .id = 'identifier')
  
  return(combined_values)
}



## CV, CHV, SV metric functions appropriated from https://github.com/ALCrofts/CABO_SVH_Forest_Sites/tree/v1.0
# cv function
calculate_cv <- function(pixel_values_df, subplots, wavelengths) {
  cv <- pixel_values_df %>%
    select(c({{subplots}}, {{wavelengths}})) %>%
    group_by({{subplots}}) %>%
    summarise_all(~sd(.)/abs(mean(.))) %>%
    rowwise({{subplots}}) %>%
    summarise(CV = sum(c_across(cols = everything()), na.rm = T) / (ncol(.) - sum(is.na(c_across(everything())))))
  
  return(cv)
}

#rarefraction cv function
cv_rf <- function(pixel_values_df, # Dataframe with spectral reflectance values 
                  subplots, # What you want to calculate spectral diversity for, here it's plots.
                  wavelengths, # Cols where spectral reflectance values are
                  n = 999, # Number of random resampling events, if rarefraction = T.
                  min_points # minimum number of pixels (ie. the min # of pixels in any subplot)
){ 
  # convert to datatable (more efficient performance)
  setDT(spectral_df)
  
  # create a list to store CV values for each replication
  cv_list <- vector("list", n)
  
  # b) calculate CV for each resampling event
  for (i in seq_len(n)) {
    # sample to the minimum number of points per subplot
    sampled_df <- spectral_df[, .SD[sample(.N, min_points)], by = areas_of_interest, .SDcols = wavelengths]
    
    # calculate CV for each wavelength within each subplot
    cv_data <- sampled_df[, lapply(.SD, function(x) sd(x) / abs(mean(x, na.rm = TRUE))), by = areas_of_interest]
    
    # sum across wavelengths and normalize by the number of bands (ignoring NAs)
    cv_data[, CV := rowSums(.SD, na.rm = TRUE) / (length(wavelengths) - rowSums(is.na(.SD))), .SDcols = wavelengths]
    
    # store cv values
    cv_list[[i]] <- cv_data[, .(CV), by = areas_of_interest]
  }
  
  # c) Collapse list of CV data tables into a single data table and calculate average CV for each area of interest
  CV <- rbindlist(cv_list)[, .(CV = mean(CV, na.rm = TRUE)), by = areas_of_interest]
  
  return(CV)
}

# sv function
calculate_sv <- function(pixel_values_df, subplots, wavelengths) {
  spectral_points <- pixel_values_df %>%
    group_by({{subplots}}) %>%
    summarise(points = n())
  
  sv <- pixel_values_df %>%
    select(c({{wavelengths}}, {{subplots}})) %>%
    group_by({{subplots}}) %>%
    summarise_all(~sum((.x - mean(.x))^2)) %>%
    rowwise({{subplots}}) %>%
    summarise(SS = sum(c_across(cols = everything()))) %>%
    left_join(spectral_points) %>%
    summarise(SV = SS / (points - 1))
  
  return(sv)
}

# chv function
calculate_chv <- function(df, dim) {
  CHV_df <- df %>%
    select(1:dim)
  
  # convert to matrix
  CHV_matrix <- as.matrix(CHV_df)
  
  # calculate chv
  CHV <- geometry::convhulln(CHV_matrix, option = "FA")
  return(CHV)
}

# function to calculate chv for each subplot
calculate_chv_for_subplots <- function(df, wavelengths, dim = 3, subplots = 'subplot_id', rarefraction = TRUE, n = 999) {
  results <- tibble(subplot_id = character(), CHV = double())
  
  # Perform PCA for specified wavelengths
  PCA <- df %>%
    select(all_of(wavelengths)) %>%
    vegan::rda(scale = FALSE)
  
  # Add subplot id as column to PCA df
  pca_results <- data.frame(PCA$CA$u) %>%
    bind_cols(subplot_id = df[[subplots]])
  
  # Compute the minimum number of points across all subplots
  min_points <- pca_results %>%
    group_by(subplot_id) %>%
    summarise(points = n()) %>%
    summarise(min_points = min(points)) %>%
    pull(min_points)
  
  # Loop through each subplot 
  for (subplot in unique(df[[subplots]])) {
    # Subset data for current subplot
    subplot_sample <- pca_results %>%
      filter(subplot_id == subplot)
    
    if (rarefraction) {
      # Resample CHV n times and calculate the mean
      chv_values <- replicate(n, {
        resampled <- subplot_sample %>%
          select(-subplot_id) %>%
          sample_n(min_points, replace = FALSE)
        
        chv_out <- calculate_chv(resampled, dim = dim)
        return(chv_out$vol)
      })
      
      mean_chv <- mean(chv_values)
    } else {
      # calculate CHV without resampling
      chv_out <- calculate_chv(subplot_sample, dim = dim)
      mean_chv <- chv_out$vol
    }
    
    # store results
    results <- results %>%
      add_row(subplot_id = subplot, CHV = mean_chv)
  }
  
  return(results)
}


## FUNCTION FOR CALCULATING ALL METRICS

calculate_metrics <- function(pixel_values_df, masked = TRUE, wavelengths) {
  results <- list()
  
  # loop through each site (represented as 'identifier' from file name)
  for (identifier in unique(pixel_values_df$identifier)) {
    
    # filter pixel values for the current identifier
    site_pixel_values <- pixel_values_df %>% filter(identifier == !!identifier)
    
    # calculate metrics (CV, SV, CHV)
    cv <- calculate_cv(site_pixel_values, subplot_id, wavelengths)
    sv <- calculate_sv(site_pixel_values, subplot_id, wavelengths)
    chv <- calculate_chv_for_subplots(site_pixel_values, wavelengths)
      
    # store results
    results[[identifier]] <- list(CV = cv, SV = sv, CHV = chv, CHV_nopca = chv_nopca)
  }
  
  # combine metrics into data frames
  combined_cv <- bind_rows(lapply(results, function(x) x$CV), .id = 'identifier')
  combined_sv <- bind_rows(lapply(results, function(x) x$SV), .id = 'identifier')
  combined_chv <- bind_rows(lapply(results, function(x) x$CHV), .id = 'identifier')
  
  # create a data frame for combined metrics
  combined_metrics <- combined_cv %>%
    left_join(combined_sv, by = c("identifier", "subplot_id")) %>%
    left_join(combined_chv, by = c("identifier", "subplot_id")) 
  
  # add image_type column based on masked argument
  if (masked) {
    combined_metrics <- combined_metrics %>%
      mutate(image_type = 'masked')
  } else {
    combined_metrics <- combined_metrics %>%
      mutate(image_type = 'unmasked')
  }
  
  return(combined_metrics)
}


# CALCULATE_FIELD_DIVERSITY -----------------------------------------------

calculate_field_diversity <- function(survey_data){
  # get unique site names
  ausplot_sites <- unique(survey_data$site_location_name) 
  ausplot_sites <- ausplot_sites[ausplot_sites != ""]
  
  # list to store results for all lists
  all_site_results <- list()
  
  # list to store community matrices - useful to check that community matrices are correct :)
  community_matrices <- list()
  
  # list to store presence absense matrices
  presence_absence_matrices <- list()

  # Loop through each unique site
  for (site in ausplot_sites) {
    # Filter data for the current site
    site_survey_data <- survey_data %>%
      filter(site_location_name == site)
    
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
    
    site_survey_data <- site_survey_data %>%
      mutate(
        point_number2 = case_when(
          transect_direction == "E-W" ~ 100 - point_number,
          transect_direction == "N-S" ~ 100 - point_number,
          TRUE ~ point_number
        ),
        transect_direction2 = case_when(
          transect_direction %in% c("W-E", "E-W") ~ "W-E",
          transect_direction %in% c("N-S", "S-N") ~ "S-N"
        )
      )
    
    # Assign plotXY coordinates based on 'transect_direction2' and 'transect_number'
    site_survey_data <- site_survey_data %>%
      mutate(
        X_plot = case_when(
          transect_direction2 == "W-E" ~ point_number2,
          transect_direction2 == "S-N" & transect_number == 1 ~ 10,
          transect_direction2 == "S-N" & transect_number == 2 ~ 30,
          transect_direction2 == "S-N" & transect_number == 3 ~ 50,
          transect_direction2 == "S-N" & transect_number == 4 ~ 70,
          transect_direction2 == "S-N" & transect_number == 5 ~ 90
        ),
        Y_plot = case_when(
          transect_direction2 == "S-N" ~ point_number2,
          transect_direction2 == "W-E" & transect_number == 1 ~ 10,
          transect_direction2 == "W-E" & transect_number == 2 ~ 30,
          transect_direction2 == "W-E" & transect_number == 3 ~ 50,
          transect_direction2 == "W-E" & transect_number == 4 ~ 70,
          transect_direction2 == "W-E" & transect_number == 5 ~ 90
        )
      )
    
    # subplot rows and columns - +1 ensures 0 point values fall into correct subplot, 
    # pmin ensures 100 point values falls in correct subplot given +1
    site_survey_data$subplot_row <- pmin(ceiling((site_survey_data$Y_plot + 1) / 20), 5)
    site_survey_data$subplot_col <- pmin(ceiling((site_survey_data$X_plot + 1) / 20), 5)
    
    # single ID for subplot row and column
    site_survey_data$subplot_id <- paste(site_survey_data$subplot_row, site_survey_data$subplot_col, sep = "_")
    
    subplot_diversity <- site_survey_data %>%
      drop_na(standardised_name) %>%
      filter(!standardised_name %in% c('Dead grass', 'Dead shrub')) %>% 
      group_by(subplot_id) %>%
      summarise(species_richness = n_distinct(standardised_name))
    
    community_matrix <- site_survey_data %>%
      drop_na(standardised_name) %>%
      filter(!standardised_name %in% c('Dead grass', 'Dead shrub')) %>% 
      count(subplot_id, standardised_name) %>%
      spread(standardised_name, n, fill = 0)
    
    
    # remove unwanted column if it exists -- WHY DOES THIS COLUMN EXIST~!!!!>???>
 #   if ("V1" %in% colnames(community_matrix)) {
  #    community_matrix <- community_matrix %>% 
   #     dplyr::select(-V1)
  #  }
    
    # store the community matrix  in the list
    community_matrices[[site]] <- community_matrix
    
    
    # calculate presence-absence values
    presence_absence_matrix <- site_survey_data %>%
      group_by(subplot_id) %>%
      summarise(presence_proportion = sum(!is.na(standardised_name)) / n()) 
    
    # store the presence-absence matrix
    presence_absence_matrices[[site]] <- presence_absence_matrix
    
    # calculate diversity indices
    shannon_diversity <- diversity(community_matrix[, -1], index = "shannon")
    simpson_diversity <- diversity(community_matrix[, -1], index = "simpson")
    
    subplot_diversity <- subplot_diversity %>%
      mutate(shannon_diversity = shannon_diversity,
             simpson_diversity = simpson_diversity,
             pielou_evenness = shannon_diversity / log(species_richness),
             exp_shannon = exp(shannon_diversity),
             inv_simpson = 1 / (simpson_diversity),
             site = site)
    
    # store  result for  current site
    all_site_results[[site]] <- subplot_diversity
  }
  
  # combine into one df
  final_results <- bind_rows(all_site_results, .id = "site")
  
  return(list(
    final_results = final_results, 
    community_matrices = community_matrices, 
    presence_absence_matrices = presence_absence_matrices  
  ))
}

