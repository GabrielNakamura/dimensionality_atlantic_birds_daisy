
name_files <- paste("ras_", c("fdi", "fev", "mpd", "pd", "pdfunc", "riq", "vpd"), ".rds", sep = "")
list_files <- vector(mode = "list", length = length(name_files))
list_m <- vector(mode = "list", length = length(name_files))
list_data <- lapply(name_files, function(x) readRDS(file = here::here("data", "processed", x)))

lapply(list_data, function(x) dim(x))
head(list_data[[5]])
head(list_data[[2]])

# minimmum extent
do.call(rbind, lapply(list_data, function(x){
  dfs_extent <- raster::extent(x) 
  lat_limits <- dfs_extent[c(3, 4)]
  # long_limits <- dfs_extent[c(1,2)]
}))

for(i in 1:length(name_files)){
  # i = 3
  df <- readRDS(file = here::here("data", "processed", name_files[i]))
  dfs <- raster::rasterFromXYZ(df)
  dfs_extent <- raster::extent(dfs) 
  lat_limits <- dfs_extent[c(3, 4)]
  lat_limits[2] <- -5.25
  long_limits <- dfs_extent[c(1,2)]
  band <- seq(lat_limits[1], lat_limits[2], by = 0.5)
  band_data <- 
  lapply(1:(length(band) - 1), function(i){
   x <-  raster::crop(dfs, raster::extent(c(long_limits, band[i], band[i + 1])))
   x <- raster::as.data.frame(x)
   #x <-  na.omit(raster::as.data.frame(x))
   return(x)
  })
  list_m[[i]] <- band_data
}

do.call(rbind, lapply(list_m, function(x){
  x[1]
})
)
unlist(lapply(list_m, function(x){
  x[2]
})
)

for(i in 1:length(list_m[[1]])){
  # i = 3
  do.call(cbind, lapply(list_m, function(x) do.call(data.frame, x[[3]])))
}
list_m2 <- vector(mode = "list", length = 7)
for(i in 1:length(list_m)){
  #i = 1
  for(j in 1:7)
  lapply(list_m, `[[`, i)[[j]]
}
lapply(list_m, `[[`, 1)
length(list_m[[1]][1])
lapply(list_m, function(x){
  x[[1]][[1]]
})
for(i in 1:length(list_m)){
  lapply(list_m, function(x){
    x[]
  })
}

lapply(list_m, function(x){
  x[1]
})

for(i in 1:7){
  i = 1
  list_m[[i]][[45]]
}
