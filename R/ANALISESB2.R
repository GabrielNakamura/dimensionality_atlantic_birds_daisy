

step <- 2

analyses <- function(presab, func, morf, filo, rast,lat_limits, long_limits, step, all_res = FALSE){
  

# adicionando o shape da Mata Atlântica #
# transformando em SpatialPolygonsDataFrame
shape <- raster::shapefile ("./Datasets/shape_mata_atlantica_IBGE_5milhoes_policonica_sirgas2000.shp")
proj4string(taxon.diver$Richness_Raster)
crs(shape)<- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
matriz <- lets.pamcrop(taxon.diver, shape, remove.sp = TRUE)

# adicionando o shape da Mata Atlântica #
# transformando em SpatialPolygonsDataFrame
shape <- raster::shapefile ("./Datasets/shape_mata_atlantica_IBGE_5milhoes_policonica_sirgas2000.shp")
proj4string(taxon.diver$Richness_Raster)
crs(shape)<- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
matriz <- lets.pamcrop(taxon.diver, shape, remove.sp = TRUE)


#CALCULANDO A DIVERSIDADE FUNCIONAL

#adicionando _ as esp?cies e nomeando as linhas com os nomes das esp?cies para usar com a matriz de presen?a/ aus?ncia

rownames(func) <- func[ , 1]

 func <- func [, -1]

dist.gower.func <- as.matrix(gowdis(func))
rownames(dist.gower.func)<- rownames(func)
colnames(dist.gower.func) <- rownames(func)

#calculo do MPD (func)
func.mpd <- mpd(presab[, -c(1:2)], dist.gower.func)
#calculo do MNTD (func)  
func.mntd <- mntd(presab[, -c(1:2)] ,dist.gower.func )
#Organizando as coordenadas 
coord.func.mpd.mntd <- cbind(presab[,1:2],func.mpd, func.mntd)
#tirando as linhas que n?o tem informa??es (ficamos com 262 esp?cies)
coord.func.mpd.mntd <- na.omit(coord.func.mpd.mntd)
coord.func.mpd.mntd <- as.data.frame(coord.func.mpd.mntd)
#Organizando as informa??es para o raster de MPD da diversidade funcional 
x <- coord.func.mpd.mntd$`Longitude(x)`
y <- coord.func.mpd.mntd$`Latitude(y)`
MPD.func <- coord.func.mpd.mntd$func.mpd
xy <- cbind(x, y)

temp <- taxon.diver$Richness_Raster
temp[!is.na(temp)] <- 0
cels <- cellFromXY(temp,xy)
temp[cels] <- MPD.func


plot(temp, xlab = "Longitude", ylab = "Latitude")
raster::writeRaster(temp, filename = "./FUNCMPD.tiff", format = "GTiff", overwrite=TRUE)

#Organizando as informa??es para o raster de MNTD da diversidade funcional 
x <- coord.func.mpd.mntd$`Longitude(x)`
y <- coord.func.mpd.mntd$`Latitude(y)`
MNTD.func <- coord.func.mpd.mntd$func.mntd
xy <- cbind(x, y)

temp <- taxon.diver$Richness_Raster
temp[!is.na(temp)] <- 0
cels <- cellFromXY(temp,xy)
temp[cels] <- MNTD.func


plot(temp, xlab = "Longitude", ylab = "Latitude")
raster::writeRaster(temp, filename = "./FUNCMNTD.tiff", format = "GTiff", overwrite=TRUE)

#CALCULO DA DIVERSIDADE FILOGEN?TICA

#matrix de dist?ncia
filo.matriz <- lapply(filo, cophenetic.phylo)

#calculando MPD
filo.mpd <- sapply(filo.matriz, function(i) mpd(presab[, -c(1:2)] , i))
average_mpd <- rowSums(filo.mpd)
sd_mpd <- apply(filo.mpd, 1, sd)
#calculando MNTD
filo.mntd <- sapply(filo.matriz, function(i) mntd(presab[, -c(1:2)], i))
average_mntd <- rowSums(filo.mntd)
sd_mntd <- apply(filo.mntd, 1, sd)
# juntando latitude, longitude, riqueza e frequ?ncia funcional
coord.filo.mpd.mntd <- cbind(taxon.diver.data.frame[,1:2], average_mpd, average_mntd, sd_mpd, sd_mntd)
#tirando as linhas que n?o tem informa??es (ficamos com 262 esp?cies)
coord.filo.mpd.mntd <- na.omit(coord.filo.mpd.mntd)
coord.filo.mpd.mntd <- as.data.frame(coord.filo.mpd.mntd)
#Organizando as informa??es para o raster de MPD da diversidade filogenetica 
x <- coord.filo.mpd.mntd$`Longitude(x)`
y <- coord.filo.mpd.mntd$`Latitude(y)`
MPD.filo <- coord.filo.mpd.mntd$average_mpd
xy <- cbind(x, y)

temp <- taxon.diver$Richness_Raster
temp[!is.na(temp)] <- 0
cels <- cellFromXY(temp,xy)
temp[cels] <- MPD.filo


plot(temp, xlab = "Longitude", ylab = "Latitude")
raster::writeRaster(temp, filename = "./FILOMPD.tiff", format = "GTiff", overwrite=TRUE)

#Organizando as informa??es para o raster de MNTD da diversidade filogenetica  
x <- coord.filo.mpd.mntd$`Longitude(x)`
y <- coord.filo.mpd.mntd$`Latitude(y)`
MNTD.filo <- coord.filo.mpd.mntd$average_mntd
xy <- cbind(x, y)

temp <- taxon.diver$Richness_Raster
temp[!is.na(temp)] <- 0
cels <- cellFromXY(temp,xy)
temp[cels] <-MNTD.filo


plot(temp, xlab = "Longitude", ylab = "Latitude")
raster::writeRaster(temp, filename = "./FILOMNTD.tiff", format = "GTiff", overwrite=TRUE)


# DIVERSIDADE MORFOLOGICA
#calculo original no meu script
rownames(morf)<- morf[,1]
morf <- morf [, -1]


matriz.dist.morf <- vegdist(morf,  method="mahalanobis")

matriz.mah <- as.matrix(matriz.dist.morf )

#calculo de MPD
morf.mpd <- mpd(presab[, -c(1:2)] , matriz.mah)
#calculo do MNTD  
morf.mntd <- mntd(presab[, -c(1:2)] , matriz.mah)
#juntando as informa??es dos indices com as coordenadas geogr?ficas
coord.morf.mpd.mntd <- cbind(presab[,1:2],morf.mpd, morf.mntd)

#tirando as linhas que n?o tem informa??es (ficamos com 262 esp?cies)
coord.morf.mpd.mntd <- na.omit(coord.morf.mpd.mntd)
coord.morf.mpd.mntd <- as.data.frame(coord.morf.mpd.mntd)
#Organizando as informa??es para o raster de MPD da diversidade morfologica 
x <- coord.morf.mpd.mntd$`Longitude(x)`
y <- coord.morf.mpd.mntd$`Latitude(y)`
MPD.morf <- coord.morf.mpd.mntd$morf.mpd
xy <- cbind(x, y)

temp <- taxon.diver$Richness_Raster
temp[!is.na(temp)] <- 0
cels <- cellFromXY(temp,xy)
temp[cels] <- MPD.morf


plot(temp, xlab = "Longitude", ylab = "Latitude")
raster::writeRaster(temp, filename = "./MORFMPD.tiff", format = "GTiff", overwrite=TRUE)

#Organizando as informa??es para o raster de MNTD da diversidade morfologica  
x <- coord.morf.mpd.mntd$`Longitude(x)`
y <- coord.morf.mpd.mntd$`Latitude(y)`
MNTD.morf <- coord.morf.mpd.mntd$morf.mntd
xy <- cbind(x, y)

temp <- taxon.diver$Richness_Raster
temp[!is.na(temp)] <- 0
cels <- cellFromXY(temp,xy)
temp[cels] <- MNTD.morf


plot(temp, xlab = "Longitude", ylab = "Latitude")
raster::writeRaster(temp, filename = "./MORFMNTD.tiff", format = "GTiff", overwrite=TRUE)

#Montando os raster's
rast.mpd.fun <- rasterize(coord.func.mpd.mntd[,c(1,2)], rast, coord.func.mpd.mntd[, 3])
rast.mntd.fun <- rasterize(coord.func.mpd.mntd[,c(1,2)], rast, coord.func.mpd.mntd[, 4])
rast.mpd.filo <- rasterize(coord.filo.mpd.mntd[,c(1,2)], rast, coord.filo.mpd.mntd[, 3])
rast.mntd.filo <- rasterize(coord.filo.mpd.mntd[,c(1,2)], rast, coord.filo.mpd.mntd[, 4])
rast.mpd.morf <- rasterize(coord.morf.mpd.mntd[,c(1,2)], rast, coord.morf.mpd.mntd[, 3])
rast.mntd.morf <- rasterize(coord.morf.mpd.mntd[,c(1,2)], rast, coord.morf.mpd.mntd[, 4])
diver <- stack(list(rast, rast.mntd.fun, rast.mpd.fun, rast.mntd.filo, rast.mpd.filo, rast.mpd.morf, rast.mntd.morf))

band <- seq(lat_limits[1], lat_limits[2], by = step)
list_band <-
  lapply(1:(length(band) - 1), function(i) {
    x <- crop(diver, extent(c(long_limits, band[i], band[i + 1])))
    x <-  na.omit(as.data.frame(x))
    return(x)
    
  })
pca.result <- lapply(list_band, function(i) prcomp( scale(i))$sdev^2)
ca_index <- sapply(pca.result, camargo_index)
latitude <- sapply(1:(length(band)-1), function(i) (band[i] + band[i + 1])/2)
area <- sapply(list_band, nrow)
reg <- lm(ca_index ~latitude + area)$coef
if(all_res == FALSE){return(reg)}else{
return(list(fun= coord.func.mpd.mntd, filo = coord.filo.mpd.mntd, morf= coord.morf.mpd.mntd, camargo_lat = cbind(ca_index, latitude), reg = reg))
}
}

camargo_index <- function(i){
  
  pro <- i/sum(i)
  camargo <- matrix(NA, nrow =  length(pro), ncol = length(pro))
  for(p in 1:(length(pro)-1)){
    for(q in (p+1):length(pro)){
    camargo[p, q] <- (pro[p] - pro[q])/length(pro)
  }
}

  return(sum(camargo, na.rm = TRUE))
}         




