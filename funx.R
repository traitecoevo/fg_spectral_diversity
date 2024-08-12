

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
  rasters <- lapply(tif_files, raster)
  
  # extract band names from file names
  band_names <- tools::file_path_sans_ext(basename(tif_files))
  
  # stack rasters and assign band names
  combined_image <- stack(rasters)
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
 # writeRaster(combined_image, filename = paste0(output_filename,'.tif'), format = "GTiff", options="INTERLEAVE=BAND", overwrite = TRUE)
  # write .envi file
  writeRaster(combined_image, filename = paste0(output_filename,'.envi'), format = "ENVI", options="INTERLEAVE=BAND", overwrite = TRUE)
  
  plot(combined_image)
}
}

# CREATE_MASKED_RASTER FUNCTION 
create_masked_raster <- function(envi_dir, mask_dir, masked_raster_dir, year,
                                 band_names, band_wavelengths,
                                 window_size = 100, NDVI_Thresh, NIR_Thresh, Blue_Thresh = 0,
                                 Blue, Red, NIR,
                                 Continuum_Removal = TRUE, TypePCA = 'SPCA', FilterPCA = TRUE, ExcludedWL = NA){
  
  envi_files <- list.files(file.path(envi_dir, paste0(year)), pattern = 'combined_image.envi$', full.names = TRUE)
  
  # if there are no files saved, print so the user knows
  print(paste("ENVI files found:", envi_files))
  
  if (length(envi_files) == 0) {
    stop("No ENVI files found.")
  }
  
  for (envi_file in envi_files){
    
    #open the .hdr file so you can add required info to it
    HDR <- read_ENVI_header(HDRpath = paste0(file_path_sans_ext(envi_file),'.hdr'))
    
    HDR$wavelength <- band_wavelengths
    HDR$'band names' <- band_names
    HDR$'wavelength units' <- 'Nanometers'
    
    write_ENVI_header(HDR = HDR, HDRpath = paste0(file_path_sans_ext(envi_file),'.hdr'))
    
    Input_Image_File <- envi_file
    
    Input_Mask_File <- perform_radiometric_filtering(Image_Path = Input_Image_File, Mask_Path = FALSE,
                                                     Output_Dir = mask_dir, TypePCA = TypePCA,
                                                     NDVI_Thresh = NDVI_Thresh, Blue_Thresh = Blue_Thresh,
                                                     NIR_Thresh = NIR_Thresh,
                                                     Blue = Blue,
                                                     Red = Red,
                                                     NIR = NIR)
    
    print(paste("Mask file created:", Input_Mask_File))
    
    mask <- read_stars(Input_Mask_File)
    mask_raster <- as(mask, "Raster")
    
    # print range of mask values 
    print(paste("Mask value range:", range(values(mask_raster), na.rm = TRUE)))
    
    mask_raster[mask_raster == 0] <- NA
    
    raster_data <- stack(envi_file)
    raster_data_masked <- mask(raster_data, mask_raster)
    
    masked_year_dir <- file.path(masked_raster_dir, year)
    if (!dir.exists(masked_year_dir)) {
      dir.create(masked_year_dir, recursive = TRUE)
    }
    
    masked_filename <- file.path(masked_year_dir, paste0(file_path_sans_ext(basename(envi_file)),'_masked'))
    
    print(paste("Saving masked raster to:", masked_filename))
    
    writeRaster(raster_data_masked, filename = masked_filename, format = "GTiff", options="INTERLEAVE=BAND", overwrite = TRUE)
  }
}


## EXTRACT PIXEL VALUES FUNCTION

extract_pixels_values <- function(raster_files, subplot_files, wavelength_names){
  
  all_pixel_values_list <- list()
  
  for (raster_file in raster_files) {
    
    # identify the string that represents the site name
    identifier <- str_extract(basename(raster_file), "^[^_]+")
    
    #choose the corresponding subplot file
    subplot_file <- subplot_files[grep(paste0('^', identifier), basename(subplot_files))]
    
    # read in subplot file and select geometries
    subplots <- read_sf(subplot_file) %>% 
      select('geometry')
    
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
  return(all_pixel_values_list)
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
calculate_chv <- function(pixel_values_df, subplots, wavelengths) {
  subplot_of_interest <- deparse(substitute(subplots))
  
  PCA <- pixel_values_df %>%
    select(c({{wavelengths}})) %>%
    rda(scale = F)
  
  CHV <- data.frame(PCA$CA$u) %>%
    select(c(1:3)) %>%
    cbind({{pixel_values_df}}) %>%
    select(c('PC1', 'PC2', 'PC3', {{subplots}})) %>%
    group_split({{subplots}}) %>%
    set_names(map(., ~unique(.[[subplot_of_interest]]))) %>%
    map(~convhulln(.x[-4], option = 'FA')) %>%
    map_dbl('vol') %>%
    as.data.frame() %>%
    rownames_to_column('subplot_id') %>%
    rename(CHV = '.')
  
  rm(subplot_of_interest) 
  
  return(CHV)
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
    chv <- calculate_chv(site_pixel_values, subplot_id, wavelengths)
    
    # store results
    results[[identifier]] <- list(CV = cv, SV = sv, CHV = chv)
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
  
  # list to store community matrices- this is a temporary step, to check that community matrices are correct :)
  community_matrices <- list()
  
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
             pielou_evenness = shannon_diversity / log(species_richness),
             site = site)
    
    # store  result for  current site
    all_site_results[[site]] <- subplot_diversity
  }
  
  # combine into one df
  final_results <- bind_rows(all_site_results, .id = "site")
  
  return(final_results)
}
