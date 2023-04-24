
# reading libraries -------------------------------------------------------

library(dplyr)
library(magrittr)
library(reshape2)
library(reshape)
library(purrr)
library(raster)
library(terra)

# reading data ------------------------------------------------------------

name_files <- paste("ras_", c("fdi", "fev", "mpd", "pd", "pdfunc", "riq", "vpd"), ".rds", sep = "")
list_files <- vector(mode = "list", length = length(name_files))
list_m <- vector(mode = "list", length = length(name_files))
list_data <- lapply(name_files, function(x) readRDS(file = here::here("data", "processed", x)))

list_data_names <- vector(mode = "list", length = length(list_data))
for(i in 1:length(list_data)){
  names <- rownames(list_data[[i]])
  list_data_names[[i]] <- data.frame(list_data[[i]], ID = names)
}


# joining all lists -------------------------------------------------------

m <-   
list_data_names %>% 
  reduce(inner_join, by = "ID")

# croping to remove duplicated coordinates
m_crop <- m[, colnames(m)[c(4, 1, 2, 3, 7, 10, 13, 16, 19, 22)]]

# renaming 
colnames(m_crop) <- c("ID", "x", "y", "ras_FDi", "ras_FEv", "ras_mpd", "ras_pd", "ras_pdfunc", "RIQUEZA", "ras_vpd")
m_crop_coords <- m_crop[, c("x", "y", "ID", colnames(m_crop)[4:ncol(m_crop)])]


# getting latitudinal bands -----------------------------------------------

dfs_extent <- raster::extent(m_crop[, c(2, 3)]) 
lat_limits <- dfs_extent[c(3, 4)]
long_limits <- dfs_extent[c(1,2)]
band <- seq(lat_limits[1], lat_limits[2], by = 0.5)

# cropping latitudinal bands
comm_bands <- lapply(band, function(i) subset(m_crop, y >= i & y <=  i + 1))

