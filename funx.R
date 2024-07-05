

# COMBINE TIFS FUNCTION (FOR CREATING MULTI BAND IMAGE) -----------------------------------------------------------------------

create_multiband_image <- function(folder_with_tifs){
# folder list | recursive = won't pick folders within folders
folders <- list.dirs(folder_with_tifs, full.names = FALSE, recursive = FALSE)
folders <- folders[folders != "point_clouds"]

# define order of bands
desired_band_order <- c("blue", "green", "red", "red_edge", "nir")

## NOTE: spectral band image tif file names must be named after their band (e.g., blue, nir, etc), 
#  otherwise change 'desired_band_order' to match file names
#  should be combined in wavelength order, esp for biodivmapR processes (i.e. as above)

# loop thru each folder 
for (folder in folders) {
  # create path
  folder_path <- file.path(folder_with_tifs, folder)
  
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
  
  # create output file
  output_filename <- file.path("data_out/combined_rasters", paste0(folder, "_combined_image.tif"))
  writeRaster(combined_image, filename = output_filename, format = "GTiff", options="INTERLEAVE=BAND", overwrite = TRUE)
  
  plot(combined_image)
}
}

# CREATE_MASKED_RASTER FUNCTION - needs work, want to add arguments e.g. wavelengths, thresholds, etc

#THIS FUNCTION ASSUMES THAT WAVELENGTHS DETECTED BY DRONE ARE:450,560,650,730,840.
create_masked_raster <- function(raster_files, envi_dir, mask_dir, masked_raster_dir, window_size = 100){
  for (raster_file in raster_files){
    
    #read in raster file
    raster_data <- stack(raster_file)
    
    envi_filename <- paste0(envi_dir, file_path_sans_ext(basename(raster_file)))
    
    # write raster as an envi file. overwrite = T so it doesn't give annoying errors erry time
    terra::writeRaster(x = raster_data, filename = envi_filename, filetype = 'ENVI', overwrite = T)
    
    #open the .hdr file so you can add required info to it
    HDR <- read_ENVI_header(HDRpath = paste0(envi_filename,'.hdr'))
    
    HDR$wavelength <- c(450,560,650,730,840)
    
    HDR$'band names' <- c("blue", "green", "red", "red_edge", "nir")
    
    HDR$'wavelength units' = 'Nanometers'
    
    write_ENVI_header(HDR = HDR, HDRpath = paste0(envi_filename,'.hdr'))
    
    name_raster <- file_path_sans_ext(basename(raster_file))
    
    tif_file <- file.path(envi_dir,name_raster)
    
    Input_Image_File <- tif_file
    
    NDVI_Thresh <- 0.05
    Blue_Thresh <- 0.15
    NIR_Thresh <- 0.02
    
    # continuum removal is a normalisation procedure which reduces multiplicative effects
    Continuum_Removal <- TRUE
    
    TypePCA <- 'SPCA'
    
    # PCA FILTERING:        Set to TRUE if you want second filtering based on PCA outliers to be processed.
    # Slower process
    # Automatically set to FALSE if TypePCA     = 'MNF'
    FilterPCA <- TRUE
    
    Excluded_WL <- NA
    
    Input_Mask_File <- perform_radiometric_filtering(Image_Path = Input_Image_File, Mask_Path = FALSE,
                                                     Output_Dir = mask_dir, TypePCA = TypePCA,
                                                     NDVI_Thresh = NDVI_Thresh, Blue_Thresh = Blue_Thresh,
                                                     NIR_Thresh = NIR_Thresh,
                                                     Blue = 450,
                                                     Red = 650,
                                                     NIR = 840)
    mask <- read_stars(Input_Mask_File)
    
    mask_raster <- as(mask, "Raster")
    
    mask_raster[mask_raster == 0] <- NA
    
    raster_data_masked <- mask(raster_data, mask_raster)
    
    masked_filename <- file.path(masked_raster_dir, paste0(name_raster,'_masked'))
    
    writeRaster(raster_data_masked, filename = masked_filename, format = "GTiff", options="INTERLEAVE=BAND", overwrite = TRUE)
  }
}

## EXTRACT PIXEL VALUES FUNCTION

extract_pixels_values <- function(raster_files, subplot_files){
  
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
    names(raster_data) <- c('blue', 'green', 'red', 'red_edge', 'nir')
    
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

calculate_metrics <- function(pixel_values_df, masked = TRUE) {
  results <- list()
  
  # loop through each site (represented as 'identifier' from file name)
  for (identifier in unique(pixel_values_df$identifier)) {
    
    # filter pixel values for the current identifier
    site_pixel_values <- pixel_values_df %>% filter(identifier == !!identifier)
    
    # calculate metrics (CV, SV, CHV)
    cv <- calculate_cv(site_pixel_values, subplot_id, c('blue', 'green', 'red', 'red_edge', 'nir'))
    sv <- calculate_sv(site_pixel_values, subplot_id, c('blue', 'green', 'red', 'red_edge', 'nir'))
    chv <- calculate_chv(site_pixel_values, subplot_id, c('blue', 'green', 'red', 'red_edge', 'nir'))
    
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

