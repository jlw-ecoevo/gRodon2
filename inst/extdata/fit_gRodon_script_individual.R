## JLW 2020 - Fit gRodon

# Load Packages ----------------------------------------------------------------

library(dplyr)
library(data.table)
library(MASS)
# library(gRodon)
library(parallel)
library(Biostrings)
library(coRdon)
library(stringi)
library(matrixStats)

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

cu <- gRodon:::getStatisticsBatch("~/gRodon/inst/extdata/madin_genomes/",
                                  mc.cores = 6,
                                  bg = "individual")
setwd("~/gRodon/inst/extdata/")
save(cu, file = "CodonStatistics_individual.rda")

# Load Growth Dataset ----------------------------------------------------------

setwd("~/gRodon2/inst/extdata/")
load("GrowthRates_Madin.rda")
load("CodonStatistics_individual.rda")
load("Accession2Species_Madin.rda")
cu <- cu %>% mutate_all(unlist)
names(d)[1] <- "Species"
d <- d %>% as.data.frame(stringsAsFactors=F)

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

# model_list <- gRodon:::fitModels(stat_data, stat_data_extremo)
#
# gRodon_model_base_madin_i <- model_list[[1]]
# gRodon_model_temp_madin_i <- model_list[[2]]
# gRodon_model_partial_madin_i <- model_list[[3]]
# gRodon_model_partial_temp_madin_i <- model_list[[4]]
# gRodon_model_meta_madin_i <- model_list[[5]]
# gRodon_model_meta_temp_madin_i <- model_list[[6]]
# lambda_milc_madin_i <- model_list[[7]]

model_list <- gRodon:::fitGCModels(stat_data, stat_data_extremo)

gRodon_model_newmeta_i <- model_list[[1]]
gRodon_model_newmeta_temp_i <- model_list[[2]]
gRodon_model_newmeta_nogc_i <- model_list[[3]]
gRodon_model_newmeta_nogc_temp_i <- model_list[[4]]
lambda_newmeta_i <- model_list[[5]]
gRodon_model_meta_madin_i <- model_list[[6]]
gRodon_model_meta_temp_madin_i <- model_list[[7]]
lambda_milc_madin_i <- model_list[[8]]

setwd("~/gRodon2/R/")
load("sysdata.rda")
save(gRodon_model_base_madin,
     gRodon_model_temp_madin,
     gRodon_model_partial_madin,
     gRodon_model_partial_temp_madin,
     gRodon_model_meta_madin,
     gRodon_model_meta_temp_madin,
     lambda_milc_madin,
     gRodon_model_base,
     gRodon_model_temp,
     gRodon_model_partial,
     gRodon_model_partial_temp,
     gRodon_model_meta,
     gRodon_model_meta_temp,
     lambda_milc,
     gRodon_model_newmeta,
     gRodon_model_newmeta_temp,
     gRodon_model_newmeta_nogc,
     gRodon_model_newmeta_nogc_temp,
     lambda_newmeta,
     gRodon_model_meta_madin_i,
     gRodon_model_meta_temp_madin_i,
     lambda_milc_madin_i,
     gRodon_model_newmeta_i,
     gRodon_model_newmeta_temp_i,
     gRodon_model_newmeta_nogc_i,
     gRodon_model_newmeta_nogc_temp_i,
     lambda_newmeta_i,
     file="sysdata.rda")
