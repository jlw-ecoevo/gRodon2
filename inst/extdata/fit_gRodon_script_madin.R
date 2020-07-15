

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

cu <- gRodon:::getStatisticsBatch("~/gRodon/inst/extdata/madin_genomes/",
                                  mc.cores = 6)
setwd("~/gRodon/inst/extdata/")
save(cu, file = "CodonStatistics_Madin.rda")

# Load Growth Dataset ----------------------------------------------------------

setwd("~/gRodon/inst/extdata/")
load("GrowthRates_Madin.rda")
load("CodonStatistics_Madin.rda")
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

model_list <- gRodon:::fitModels(stat_data, stat_data_extremo)

gRodon_model_base_madin <- model_list[[1]]
gRodon_model_temp_madin <- model_list[[2]]
gRodon_model_partial_madin <- model_list[[3]]
gRodon_model_partial_temp_madin <- model_list[[4]]
gRodon_model_meta_madin <- model_list[[5]]
gRodon_model_meta_temp_madin <- model_list[[6]]
lambda_milc_madin <- model_list[[7]]

setwd("~/gRodon/R/")
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
     file="sysdata.rda")
