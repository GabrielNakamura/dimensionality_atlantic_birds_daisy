########### MAPAS DIMENSIONALIDADE ###########



# Shafile Forest Atlantic -------------------------------------------------
shp_mt = shapefile("shape_mata_atlantica_IBGE_5milhoes_policonica_sirgas2000.shp")
#transformando a projeção do shapefile da MT igual ao raster 
shp_mt = spTransform(shp_mt, CRS("+proj=longlat +datum=WGS84 +no_defs 
                         s") )
shp_mt = st_as_sf(shp_mt)


# RICHNESS ----------------------------------------------------------------
ras_riq = as.data.frame(riq, xy = TRUE)
ras_riq$RIQUEZA[ras_riq$RIQUEZA == 0] = NA
ras_riq = na.omit(ras_riq)
ggplot()+
  geom_raster(aes(x = x, y=y, fill = RIQUEZA),
              data = ras_riq)  +
  geom_sf(fill = 'transparent',data = shp_mt) +
  theme_classic() +
  scale_fill_viridis_b(name = 'RIQUEZA', direction = -1) +
  theme_void() 


# PHYLOGENTIC -------------------------------------------------------------


#raster MPD transformado em dataframe
rasmpd_frame = as.data.frame(ras_mpd, xy= TRUE)
# transformando 0 em NA
rasmpd_frame$ras_mpd[rasmpd_frame$ras_mpd == 0] = NA
rasmpd_frame = na.omit(rasmpd_frame)

#MAPA MPD
ggplot()+
  geom_raster(aes(x = x, y=y, fill = ras_mpd),
              data = rasmpd_frame)  +
  geom_sf(fill = 'transparent',data = shp_mt) +
  theme_classic() +
  scale_fill_viridis_b(name = 'MPD', direction = -1) +
  theme_void() 
  
#MAPA VPD
rasvpd_frame = as.data.frame(ras_vpd, xy = TRUE)
rasvpd_frame$ras_vpd[rasvpd_frame$ras_vpd == 0] = NA
rasvpd_frame = na.omit(rasvpd_frame)

ggplot()+
  geom_raster(aes(x = x, y=y, fill = ras_vpd),
              data = rasvpd_frame)  +
  geom_sf(fill = 'transparent',data = shp_mt) +
  theme_classic() +
  scale_fill_viridis_b(name = 'VPD', direction = -1) +
  theme_void() 

# MAPA PD
raspd_frame = as.data.frame(ras_pd, xy = TRUE)
raspd_frame$ras_pd[raspd_frame$ras_pd == 0] = NA
raspd_frame = na.omit(raspd_frame)

ggplot()+
  geom_raster(aes(x = x, y=y, fill = ras_pd),
              data = raspd_frame)  +
  geom_sf(fill = 'transparent',data = shp_mt) +
  theme_classic() +
  scale_fill_viridis_b(name = 'PD', direction = -1) +
  theme_void() 

# MAPA PD FUNC ------------------------------------------------------------

# MAPA PD
ras_pdfunc_frame = as.data.frame(ras_pdfunc, xy = TRUE)
ras_pdfunc_frame$ras_pdfunc[ras_pdfunc_frame$ras_pdfunc == 0] = NA
ras_pdfunc_frame = na.omit(ras_pdfunc_frame)

ggplot()+
  geom_raster(aes(x = x, y=y, fill = ras_pdfunc),
              data = ras_pdfunc_frame)  +
  geom_sf(fill = 'transparent',data = shp_mt) +
  theme_classic() +
  scale_fill_viridis_b(name = 'PD', direction = -1) +
  theme_void() 
