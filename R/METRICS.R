########## Calculo das métricas ########


# Phylogenetic Diversity --------------------------------------------------

tc_ma = as.matrix(tc_frame[, -c(1:2)])
#adicionando nomes nas linhas 
rownames(tc_ma) = as.numeric(1:324)

obj = comparative.comm(
  phy= phylo,
  comm= tc_ma,
  traits = NULL,
  env = NULL,
  warn = TRUE,
  force.root = -1
)


#VPD
vpd = .vpd(obj)
# Transformando e dataframe
vpd_frame = as.data.frame(vpd)
vpd_f = rownames_to_column(vpd_frame)
write.csv(vpd_f, "vpd_f.csv")
vpd_frame = read.csv("vpd_f.csv")
write.csv(vpd_frame, "vpd_frame.csv")
# Junção das coordendas com os valores de MPD
lonlat = tc_matrix$Presence_and_Absence_Matrix[, 1:2]
coorvpd = cbind(lonlat, vpd_frame)
coorvpd = na.omit(coorvpd)
# coordenadas
x = coorvpd$`Longitude(x)`
y = coorvpd$`Latitude(y)`
xy = cbind(x,y)
#preenchendo células 
temp <- tc_matrix$Richness_Raster
temp[!is.na(temp)] <- 0
cels <- cellFromXY(temp,xy)
temp[cels] <- coorvpd$vpd
# PLOT MAP
ras_vpd =  plot(temp, xlab = "Longitude", ylab = "Latitude")
# Salvando raster
ras_vpd =raster::writeRaster(temp, filename = "ras_vpd", 
                             format = "GTiff", overwrite=TRUE)

#MPD
mpd = .mpd(obj)
# Transformando e dataframe
mpd_frame = as.data.frame(mpd)
# Junção das coordendas com os valores de MPD
lonlat = tc_matrix$Presence_and_Absence_Matrix[, 1:2]
coormpd = cbind(lonlat, mpd_frame)
coormpd = na.omit(coormpd)
# coordenadas
x = coormpd$`Longitude(x)`
y = coormpd$`Latitude(y)`
xy = cbind(x,y)
#preenchendo células 
temp <- tc_matrix$Richness_Raster
temp[!is.na(temp)] <- 0
cels <- cellFromXY(temp,xy)
temp[cels] <- coormpd$mpd

# PLOT MAP
ras_mpd =  plot(temp, xlab = "Longitude", ylab = "Latitude")
# Salvando raster
ras_mpd =raster::writeRaster(temp, filename = "ras_mpd", 
                    format = "GTiff", overwrite=TRUE)


# PD
matriz_pa = as.matrix(tc_frame[, -c(1,2)])
sespd = ses.pd(samp=matriz_pa, 
               tree= phylo, null.model = 'taxa.labels',
               runs=1000, iterations = 1000)
pd = as.data.frame(sespd$pd.obs.z)
# Junção das coordendas com os valores de PD
lonlat = tc_matrix$Presence_and_Absence_Matrix[, 1:2]
coorpd = cbind(lonlat, pd)
coorpd = na.omit(coorpd)
# coordenadas
x = coorpd$`Longitude(x)`
y = coorpd$`Latitude(y)`
xy = cbind(x,y)
#preenchendo células 
temp <- tc_matrix$Richness_Raster
temp[!is.na(temp)] <- 0
cels <- cellFromXY(temp,xy)
temp[cels] <- coorpd$`sespd$pd.obs.z`

# PLOT MAP
ras_pd =  plot(temp, xlab = "Longitude", ylab = "Latitude")
# Salvando raster
ras_pd =raster::writeRaster(temp, filename = "ras_pd", 
                             format = "GTiff", overwrite=TRUE)

# FUNCTIONAL ----------------------------------------------------

#Selecionando atributos funcionais
fc = tfm_coor %>% distinct(Species, Diet.5Cat, Body_length.mm.mean, Wing_length.mm.mean,
                           Tail_length.mm.mean, Tarsus_length.mm.mean, Bill_length.mm.mean,
                           Bill_depth.mm.mean, Bill_width.mm.mean)

distgower = gowdis(fc)
distgower = as.dist(distgower)
#Dendograma

dendo = hclust(distgower, method= 'average')
#transforma o dendrograma em obj phylo
phylodendo = as.phylo(dendo)
# modificando o nome das espécies
tippn = as.data.frame(phylodendo$tip.label)
tippn = cbind(tippn, as.data.frame(phylo$tip.label))
phylodendo = sub.taxa.label(phylodendo, tippn)

#aplica o ses.pd
pdfunc = ses.pd(samp= matriz_pa, 
                tree= phylodendo,
                null.model = 'taxa.labels',
                runs=1000, iterations = 1000)
pdfunc = as.data.frame(pdfunc$pd.obs.z)
# Junção das coordendas com os valores de PD
lonlat = tc_matrix$Presence_and_Absence_Matrix[, 1:2]
coorpdfunc = cbind(lonlat, pdfunc)
coorpdfunc = na.omit(coorpdfunc)
# coordenadas
x = coorpdfunc$`Longitude(x)`
y = coorpdfunc$`Latitude(y)`
xy = cbind(x,y)
#preenchendo células 
temp <- tc_matrix$Richness_Raster
temp[!is.na(temp)] <- 0
cels <- cellFromXY(temp,xy)
temp[cels] <- coorpdfunc$`pdfunc$pd.obs.z`

# PLOT MAP
ras_pdfunc =  plot(temp, xlab = "Longitude", ylab = "Latitude")
# Salvando raster
ras_pdfunc =raster::writeRaster(temp, filename = "ras_pdfunc", 
                            format = "GTiff", overwrite=TRUE)

#FDi
mtfun = as.matrix(distgower)
dimnames(mtfun) = dimnames(matriz_pa)
rownames(mtfun) = colnames(matriz_pa)
medfunc <- dbFD(mtfun, matriz_pa)

FDi= as.data.frame(medfunc$FDiv)
# Junção das coordendas com os valores de MPD
lonlat = tc_matrix$Presence_and_Absence_Matrix[, 1:2]
coorFDi = cbind(lonlat, FDi)
coorFDi = na.omit(coorFDi)
# coordenadas
x = coorFDi$`Longitude(x)`
y = coorFDi$`Latitude(y)`
xy = cbind(x,y)
#preenchendo células 
temp <- tc_matrix$Richness_Raster
temp[!is.na(temp)] <- 0
cels <- cellFromXY(temp,xy)
temp[cels] <- coorFDi$`medfunc$FDiv`
# PLOT MAP
ras_FDi =  plot(temp, xlab = "Longitude", ylab = "Latitude")
# Salvando raster
ras_FDi =raster::writeRaster(temp, filename = "ras_FDi", 
                             format = "GTiff", overwrite=TRUE)

#FEv
FEv = as.data.frame(medfunc$FEve)
# Junção das coordendas com os valores de MPD
lonlat = tc_matrix$Presence_and_Absence_Matrix[, 1:2]
coorFEv = cbind(lonlat, FEv)
coorFEv = na.omit(coorFEv)
# coordenadas
x = coorFEv$`Longitude(x)`
y = coorFEv$`Latitude(y)`
xy = cbind(x,y)
#preenchendo células 
temp <- tc_matrix$Richness_Raster
temp[!is.na(temp)] <- 0
cels <- cellFromXY(temp,xy)
temp[cels] <- coorFEv$`medfunc$FEve`
# PLOT MAP
ras_FEv =  plot(temp, xlab = "Longitude", ylab = "Latitude")
# Salvando raster
ras_FEv =raster::writeRaster(temp, filename = "ras_FEv", 
                             format = "GTiff", overwrite=TRUE)

