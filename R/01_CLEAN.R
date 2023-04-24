ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("taxize", "FD", "reshape", "iterators",
              "dplyr", "RStoolbox", "devtools", "maptools", 
              "ape", "vegan", "picante", "rgdal", "raster", 
              "hillR", "factoextra", "FactoMineR", "gower", 
              "cluster", "StatMatch" , "bdc" , "sf",
              "data.table","ggplot2","geobr", "tmaptools", "terra",        
              "rgbif", "tibble", "spData", "pez", "letsR", "rnaturalearth", 
              "ggthemes", "phytools", "phangorn", "stats")
ipak(packages)
rm(ipak,packages)

# Directory ---------------------------------------------------------------

setwd("G:/Meu Drive/master's article")
taxon <- here::here("data", "ATLANTIC_BIRDS_qualitative.csv")

# Taxonomic Clean (150, 423)---------------------------------------------------------
taxon = read.csv("ATLANTIC_BIRDS_qualitative.csv", fileEncoding = "latin1", sep = ";")
taxon = taxon %>% dplyr::select(Record_id, Latitude_y, Longitude_x,
                         Species, Country, State, OlsonG200r)
dim(tc)
#coordinates (149974)
tc = taxon %>% filter(complete.cases(Longitude_x, Latitude_y))
# duplicats (91725)
tc = tc %>% distinct(Species, Latitude_y, Longitude_x, .keep_all = TRUE)
#deleting characters (91706)
tc = tc %>% filter(!grepl("sp.$", Species))
#Brasil (90492)
tc = tc %>% filter(!grepl("0", Country))
#Atlantic Forest (66372)
tc= tc %>% filter(OlsonG200r=="Atlantic Forests" | OlsonG200r=="Atlantic Dry Forests")

conf_tc <- gnr_resolve(as.character(unique(tc$Species)),data_source_ids = 11, canonical = TRUE)

# For ---------------------------------------------------------------------


for (i in 1:length(tc$Species)) {
  if (tc$Species[i] == "Glyphorhynchus spirurus"){
    tc$Species[i] = "Glyphorynchus spirurus"}
  if (tc$Species[i] == "Phylloscartes cecilliae"){
    tc$Species[i] = "Phylloscartes ceciliae"}
  if (tc$Species[i] == "Myrmoderus ruficaudus"){
    tc$Species[i] = "Myrmoderus ruficauda"}
  if (tc$Species[i] == "Amadonastur lacernulatus"){
    tc$Species[i] = "Leucopternis lacernulatus"}
  if (tc$Species[i] == "Lanio cristatus"){
    tc$Species[i] = "Loriotus cristatus"}
  if (tc$Species[i] == "Tyranniscus burmeisteri"){
    tc$Species[i] = "Phyllomias burmeisteri"}
  if (tc$Species[i] == "Urubitinga urubitinga"){
    tc$Species[i] = "Buteogallus urubitinga"}
  if (tc$Species[i] == "Aramides cajaneus"){
    tc$Species[i] = "Aramides cajanea"}
  if (tc$Species[i] == "Hydropsalis forcipata"){
    tc$Species[i] = "Macropsalis creagra"}
  if (tc$Species[i] == "Hydropsalis parvula"){
    tc$Species[i] = "Caprimulgus parvulus"}
  if (tc$Species[i] == "Buteogallus coronatus"){
    tc$Species[i] = "Harpyhaliaetus coronatus"}
  if (tc$Species[i] == "Heliodoxa rubriacauda"){
    tc$Species[i] = "Clytolaema rubricauda"}
}

############# FUNCTION DIVERSITY #############
# Data base ---------------------------------------------------------------
func = read.csv("BirdFuncDat.csv", sep = ";")
# Selecting species and guilda  -------------------------------------------------------------
func = func[, c(8, 20)]

# Selecting morphologyc datas ---------------------------------------------

morf = read.csv("ATLANTIC_BIRD_TRAITS_Spp_Info.csv")

# merge between morf and func ---------------------------------------------

funmor = merge(func, morf, by.x = "Scientific", by.y = "Binomial")
# merge between taxonomic and function diversity --------------------------
tc_unique = as.data.frame(unique(tc$Species))
colnames(tc_unique)[1] = "Species"
tfm = merge(tc_unique, funmor, by.x = "Species", 
           by.y = "Scientific" )

tfm_unique = as.data.frame(unique(tfm$Species))
colnames(tfm_unique)[1] = "Species"

# Joing information -------------------------------------------------------
tfm_coor = merge(tfm,tc,  by.x = "Species", by.y = "Species")

tfm_coor = as.data.frame(tfm_coor)
#Selecting informations 
tfm_coor = tfm_coor[, c(1, 2, 15, 22,  29, 36, 43, 50, 57, 64, 65)]
# excluindo células sem informação (310 espécies)
tfm_coor = na.omit(tfm_coor)
write.csv(unique(tfm_coor$Species), "tfmuni.csv")
# PHYLOGENETIC DIVERSITY --------------------------------------------------
phylo = read.tree("AllBirdsEricson1.tre")
# máxima credibilidade
phylo = maxCladeCred(phylo, tree = TRUE, part = NULL, rooted = TRUE)
phylo$tip.label = gsub("_", " ", phylo$tip.label)
 # cortando a árvore
phylo = keep.tip(phylo, tips)




# Matrix presence/ ausence ------------------------------------------------

# Objetos para o raster de riqueza ----------------------------------------
x = tfm_coor$Longitude_x
y = tfm_coor$Latitude_y
species = tfm_coor$Species
xy = cbind(x, y)
#raster riqueza, matriz de p/a, 
tc_matrix =lets.presab.points(xy, species, xmn = -60, xmx =  -33, ymn = -35, ymx = 0, resol = 1/2, remove.cells = TRUE  )
#transformando minha matriz em data frame
tc_frame = as.data.frame(tc_matrix[["Presence_and_Absence_Matrix"]])
#soma das presenças das espécies para cada célula
tc_sum = rowSums(tc_frame[,-c(1,2)])#muitas células com apenas uma espécie
tc_sum = as.data.frame(tc_sum)
tc_sum = cbind(tc_frame[, c(1,2)], tc_sum)
plot(tc_matrix, xlab = "Longitude", ylab = "Latitude")


riq = raster::writeRaster(tc_matrix$Richness_Raster, filename = "./RIQUEZA.tiff", format = "GTiff", overwrite = TRUE)





