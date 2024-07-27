

## JLW 2020 - Fit gRodon

# Load Packages ----------------------------------------------------------------

library(dplyr)
library(data.table)
library(MASS)
library(gRodon)
library(parallel)

# Helper Functions -------------------------------------------------------------

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

rgrep <- function(big,small_vec){
  small_vec[lapply(small_vec,grepl,x=big) %>% unlist()]
}

# Calculate Codon Usage Statistics ---------------------------------------------

cu_AOA_NOB <- gRodon:::getStatisticsBatch("C:/Users/jlwei/Documents/AOA_NOB/",
                                  mc.cores = 1)
setwd("C:/Users/jlwei/Documents/gRodon2/inst/extdata/")
load("CodonStatistics_GC.rda")
cu <- rbind(cu,cu_AOA_NOB)

# Load Growth Dataset ----------------------------------------------------------

load("GrowthRates_Madin.rda")
load("Accession2Species_Madin.rda")
cu <- cu %>% mutate_all(unlist)
names(d)[1] <- "Species"
d <- d %>% as.data.frame(stringsAsFactors=F)
d <- rbind(d,data.frame(Species=c("Nitrososphaera viennensis",
                                  "Candidatus Nitrosocosmicus exaquare",
                                  "Nitrosopumilus maritimus SCM1",
                                  "Nitrolancea hollandica",
                                  "Candidatus Nitrotoga arctica",
                                  "Nitrococcus mobilis",
                                  "Nitrospina gracilis",
                                  "Nitrospinae watsonii"),
                        d=c(29.2,
                            51.7,
                            21.32,
                            36.16,
                            54,
                            9.79,
                            23.77,
                            23.1),
                        OptTemp=c(45,
                                  33,
                                  28,
                                  40,
                                  16,
                                  28,
                                  28,
                                  28),
                        GrowthTemp=NA,
                        Extremophile=FALSE))
spp_acc <- rbind(spp_acc,data.frame(V1=c("GCF_918378365.1",
                                         #"GCF_000698785.1",
                                         #"GCF_001870125.1",
                                         #"GCF_000018465.1",
                                         #"GCF_000297255.1",
                                         #"GCF_000153205.1",
                                         #"GCF_000341545.2",
                                         "GCF_946900835.1"),
                                    V2=c("Candidatus Nitrotoga arctica",
                                         #"Nitrososphaera viennensis",
                                         #"Candidatus Nitrosocosmicus exaquare",
                                         #"Nitrosopumilus maritimus SCM1",
                                         #"Nitrolancea hollandica",
                                         #"Nitrococcus mobilis",
                                         #"Nitrospina gracilis",
                                         "Nitrospinae watsonii")))

# Merge datasets
rownames(spp_acc) <- spp_acc$V1 %>% gsub(pattern="[.].*",replace="")
cu$Accession <- cu$File %>% gsub(pattern="[.].*",replace="")
cu$Spp <- spp_acc[cu$Accession,"V2"]
cu$Species <- lapply(cu$Spp,rgrep,small_vec=d$Species) %>%
  lapply("[",1) %>% unlist()
cu$Species[cu$Spp %in% d$Species] <- cu$Spp[cu$Spp %in% d$Species]
cu <- merge.easy(cu,d,key="Species") %>% subset(!is.na(Species))

# Average CUB estimates over species
stat_data <- cu %>%
  subset(Extremophile == FALSE) %>%
  group_by(Species) %>%
  summarise_all(mean,na.rm=T) %>%
  subset(!is.na(Species))

# Average CUB estimates over species, including extremophiles
stat_data_extremo <- cu %>%
  group_by(Species) %>%
  summarise_all(mean,na.rm=T) %>%
  subset(!is.na(Species))
stat_data_extremo$OGT <- stat_data_extremo$OptTemp
stat_data_extremo$OGT[is.na(stat_data_extremo$OGT)] <-
  stat_data_extremo$GrowthTemp[is.na(stat_data_extremo$OGT)]

# Fit Models -------------------------------------------------------------------

model_list <- gRodon:::fitModels(stat_data=stat_data, stat_data_extremo=stat_data_extremo)

gRodon_model_base_AOANOB <- model_list[[1]]
gRodon_model_temp_AOANOB <- model_list[[2]]
gRodon_model_partial_AOANOB <- model_list[[3]]
gRodon_model_partial_temp_AOANOB <- model_list[[4]]
gRodon_model_meta_AOANOB <- model_list[[5]]
gRodon_model_meta_temp_AOANOB <- model_list[[6]]
lambda_milc_AOANOB <- model_list[[7]]

setwd("C:/Users/jlwei/Documents/gRodon2/R/")
load("sysdata.rda")
save(gRodon_model_base,
     gRodon_model_base_euk,
     gRodon_model_base_euk_i,
     gRodon_model_base_madin,
     gRodon_model_meta,
     gRodon_model_meta_madin,
     gRodon_model_meta_madin_i,
     gRodon_model_meta_temp,
     gRodon_model_meta_temp_madin,
     gRodon_model_meta_temp_madin_i,
     gRodon_model_newmeta,
     gRodon_model_newmeta_i,
     gRodon_model_newmeta_nogc,
     gRodon_model_newmeta_nogc_i,
     gRodon_model_newmeta_nogc_temp,
     gRodon_model_newmeta_nogc_temp_i,
     gRodon_model_newmeta_temp,
     gRodon_model_newmeta_temp_i,
     gRodon_model_partial,
     gRodon_model_partial_madin,
     gRodon_model_partial_temp,
     gRodon_model_partial_temp_madin,
     gRodon_model_temp,
     gRodon_model_temp_euk,
     gRodon_model_temp_euk_i,
     gRodon_model_temp_madin,
     lambda_milc,
     lambda_milc_euk,
     lambda_milc_euk_i,
     lambda_milc_madin,
     lambda_milc_madin_i,
     lambda_newmeta,
     lambda_newmeta_i,
     gRodon_model_base_AOANOB,
     gRodon_model_temp_AOANOB,
     gRodon_model_partial_AOANOB,
     gRodon_model_partial_temp_AOANOB,
     gRodon_model_meta_AOANOB,
     gRodon_model_meta_temp_AOANOB,
     lambda_milc_AOANOB,
     file="sysdata.rda")
