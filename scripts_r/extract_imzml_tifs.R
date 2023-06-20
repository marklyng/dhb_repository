#### Extract non-compressed tiff from imzML 

#### Packages ####
library(here)
library(tidyverse)
library(Cardinal)
library(ijtiff)

#### Functions ####
# Takes a Cardinal::MSImagingExperiment, dataframe of mz-values, and an export directory and returns 32-bit tiff-images
# for import in e.g. imageJ
extract_peak_imgs <- function(msi_data, features, out_dir, overwrite = T, bits_per_sample = 32) {
  # Map over unique mz-values
  map(unique(pull(features, mz)), function(mz) {
    mz_formatted <- format(round(mz, 4), nsmall = 4) # Format mz-value to four decimal points
    dir.create(str_c(out_dir, "/", mz_formatted)) # Create a directory for each mz-value
    
    # Convert intensity array to tif  
    tif <- slice(msi_data, mz = mz) %>%
      ijtiff_img()
    
    # Format image name
    image_name <- str_c(msi_data@metadata$name, "_", mz_formatted)
    
    # Write image
    write_tif(tif, str_c(out_dir, "/", mz_formatted, "/", image_name, ".tif"),
              overwrite = overwrite,
              msg = F, # still outputs messages at the end?
              bits_per_sample = bits_per_sample
    )
  })
}

#### Data import ####
# Directory containing .imzML and .ibd files. Prefixes MUST be identical (e.g. XX.imzML and XX.ibd, NOT XX.imzML and XY.ibd)
in_dir <- here("data/msi/imzml") 
in_files <- list.files(in_dir, pattern = "imzML") # Get input file names

# Read imzML
sample_list <- map(list.files(in_dir, pattern = "imzML", full.names = T), readMSIData) %>%
  set_names(str_remove(in_files, ".imzML"))

# List of mz-values to extract. Must contain a column named "mz"
# Here I've downloaded the metaSPACE annotation list as a csv
features <- read_csv(here("md/metaspace_annotations.csv"), skip = 2) 

# Where to put the images
out_dir <- here("images/msi_280322/msi_processed")

#### Data wrangle ####

# The feature list can contain any list of mz-values, as long as they exist in the data. If they do not, no image will be exported for the non-existing mz
# E.g. you can filter the metaSPACE annotation .csv for specific molecule names
features_namefiltered <- features %>% 
  filter(str_detect(moleculeNames, pattern = "Bacillibactin|Corynebactin|Viscosin|Pyoverdin"))

# Or you can filter on specific mz-values
features_mzfiltered <- features %>% 
  filter(mz > 1200 & mz < 1500)



# Normalize (can be "tic" or "rms) - takes A WHILE
sample_list_norm <- map(sample_list, function(x) {
  x %>% 
    Cardinal::normalize(method = "rms") %>% 
    process()
})


#### Extract tifs ####
map(sample_list, function(sample) {
  extract_peak_imgs(sample, features, out_dir)
})








#### imzML metadata ####
# Get metadata from the imzml file (e.g. pixel resolution)
# imzml_md <- read_delim(here(in_dir, str_c(sample_name, ".imzML")),
#                        col_names = F)
# 
# pixel_size_x <- imzml_md %>% 
#   filter(str_detect(X3, "1000046")) %>% 
#   mutate(X3 = str_extract(X3, pattern = "(?<=value=\")[:digit:]*")) %>% 
#   pull()
# 
# pixel_size_y <- imzml_md %>% 
#   filter(str_detect(X3, "1000047")) %>% 
#   mutate(X3 = str_extract(X3, pattern = "(?<=value=\")[:digit:]*")) %>% 
#   pull()
# 
# pixel_size_unit <- imzml_md %>% 
#   filter(str_detect(X3, "1000045")) %>% 
#   mutate(X3 = str_extract(X3, pattern = "(?<=unitName=\")[:alpha:]*")) %>% 
#   pull()
