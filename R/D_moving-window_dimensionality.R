
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

# obtaining communities in moving window ----------------------------------


rast_m <- raster::rasterFromXYZ(xyz = m_crop[, c("x", "y", "ID", 
                                                 colnames(m_crop)[4:ncol(m_crop)])], crs = "+proj=longlat +datum=WGS84 +no_defs s")

mw_1 <- extract(rast_m, m_crop_coords[, c("x", "y")], buffer= 1000000, na.rm = TRUE) # buffer with specified size in meters

lapply(mw_1, function(x){
  x[!is.na(x), ]
})
x[!is.na(mw_1[[1]]), ]
mw_1[[1]][apply(mw_1[[1]], 1, Compose(is.finite, all)),]


mw_1_clean <- 
  lapply(mw_1, function(x){
    x[apply(x, 1, function(y) any(!is.na(y))), ]
  })

unlist(lapply(mw_1_clean, function(x) dim(x)[1])) # vector with the number of community in each buffer



# calculating dimensionality from moving window units ---------------------

devtools::install_github("GabrielNakamura/Dimensionality_package")
library(Dimensionality)

# EE values
EE_mw1 <-
  lapply(mw_1_clean, function(x){
    Dimensionality::EvennessEigen(matrix.M = x[, 2:ncol(x)], 
                                  scale = FALSE, 
                                  evenness = "Camargo")
  }
  )
IV_mw1 <-
  lapply(mw_1_clean, function(x){
    Dimensionality::EvennessEigen(matrix.M = x[, 2:ncol(x)], 
                                  scale = FALSE, 
                                  evenness = "Camargo")
  }
  )

IVs_all_mw1 <- 
  lapply(mw_1_clean, function(x){
    colSums(Dimensionality::ImportanceVal(matrix.M = x[, -1], scale = TRUE)$IV.obs_stopRule)
  })
IVs_all_mw1_sum <- do.call(rbind, IVs_all_mw1)

colnames(IVs_all_mw1_sum) <- paste("IV", colnames(IVs_all_mw1_sum), sep = "_")



EE_mw1 <- unlist(EE_mw1)
m_crop_org <- m_crop[, c("x", "y", "ID", colnames(m_crop)[4:ncol(m_crop)])]
m_crop_org <- data.frame(m_crop_org, EE_mw1, IVs_all_mw1_sum)
plot(rast(m_crop_org, type = "xyz"))


# getting latitudinal bands -----------------------------------------------

dfs_extent <- raster::extent(m_crop[, c(2, 3)]) 
lat_limits <- dfs_extent[c(3, 4)]
long_limits <- dfs_extent[c(1,2)]
band <- seq(lat_limits[1], lat_limits[2], by = 0.5)

# cropping latitudinal bands
comm_bands <- lapply(band, function(i) subset(m_crop, y >= i & y <=  i + 1))

