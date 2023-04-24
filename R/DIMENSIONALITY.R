########## DIMENSIONALITY ############

# Matrix with metrics -----------------------------------------------------

metrics = cbind(tc_sum, vpd_frame, mpd_frame, pd, pdfunc, FDi, FEv)
metrics = na.omit(metrics)
metrics = metrics [, -c(1,2)]

# Install Packages --------------------------------------------------------


devtools::install_github("GabrielNakamura/Dimensionality_package", force = TRUE)

## Dimensionality analysis (EE and IV)


# Read libraries ----------------------------------------------------------

library(Dimensionality)


# dimensionality analysis  ------------------------------------------------

EE_metric = Dimensionality::EvennessEigen(matrix.M = as.matrix(metrics),
                                           scale = TRUE,
                                           method = "standardize",
                                           evenness = "Camargo")
IVs_diversity = Dimensionality::ImportanceVal(matrix.M = as.matrix(metrics),
                                               scale= TRUE,
                                               method = "max",
                                               stopRule = TRUE)

IVs_diversity$IV.obs_stopRule


# by realm ----------------------------------------------------------------

reino <- metric_matrix$reino
dimensionlity_all_realms <-
  matrix(data = unlist(lapply(unique(reino),
                              function(x){
                                EE_metric_realm <- Dimensionality::EvennessEigen(matrix.M = as.matrix(div_metrics[which(metric_matrix$reino == x), ]),
                                                                                 scale = TRUE,
                                                                                 method = "standardize",
                                                                                 evenness = "Camargo")
                                IVs_diversity_realm <- Dimensionality::ImportanceVal(matrix.M = as.matrix(div_metrics[which(metric_matrix$reino == x), ]),
                                                                                     scale= TRUE,
                                                                                     method = "max",
                                                                                     stopRule = TRUE)
                                if(is.matrix(IVs_diversity_realm$IV.obs_stopRule) == TRUE){
                                  IVs_total_realm <- colSums(IVs_diversity_realm$IV.obs_stopRule)
                                } else {
                                  IVs_total_realm<- IVs_diversity_realm$IV.obs_stopRule
                                }
                                dimensionality_obs <- c(EE_metric_realm, IVs_total_realm)
                                return(dimensionality_obs)
                              }
  )
  ),
  nrow = length(unique(reino)),
  ncol = ncol(div_metrics) + 1,
  byrow = TRUE,
  dimnames = list(unique(reino), c("EE", colnames(div_metrics)
  )
  )
  )

EEs<- dimensionlity_all_realms[,1]

write.table(dimensionlity_all_realms, "dimensionality_all_realms.txt")
write.table(EEs, "EEs.txt")

###Para calcular equitabilidade dos IVs

#carregando a função
camargo.eveness <- function(n_spec, include_zeros = T){
  if(is.vector(n_spec)==FALSE){
    stop("\n n_spec must be a vector of abundance of species \n")
  }
  if (include_zeros){
    n <- n_spec
  }  else{
    n <- n_spec[n_spec > 0]
  }
  
  S <- length(n)
  camar<-matrix(nrow=length(n), ncol=length(n))
  for (i in 1:S)
  {
    for (j in 1:S)
    {
      p_i <- n[i]/sum(n)
      p_j <- n[j]/sum(n)
      camar[i,j] <- ((abs(p_i - p_j))/S)
    }
  }
  sum.camar<- abs(sum(as.dist(camar, diag= FALSE, upper= FALSE)))
  return(1-sum.camar)
}


#calculando a equitabilidade

valores<- numeric(length = 28)

for (i in 1:28){
  
  vector<- dimensionlity_all_realms[i, 2:8]
  valores[i]<- camargo.eveness(n_spec = vector, include_zeros = T)
  
}

write.table(valores, "equitabilidadeIVs.txt")


