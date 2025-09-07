

## JLW 2025 - Fit gRodon

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

# Load Data --------------------------------------------------------------------

load("CodonStatistics_madin_trim.rda")
load("GrowthRates_Madin.rda")
load("Accession2Species_Madin.rda")

cu250 <- genome_df %>%
  subset(!grepl("Error",File)) %>%
  mutate(MILC=MILC.250 %>% unlist() %>% as.numeric(),
         ENCprime=ENCprime.250 %>% unlist() %>% as.numeric(),
         B=B.250 %>% unlist() %>% as.numeric(),
         SCUO=SCUO.250 %>% unlist() %>% as.numeric(),
         MCB=MCB.250 %>% unlist() %>% as.numeric(),
         nHE=nHE  %>% unlist() %>% as.numeric(),
         File=File  %>% unlist() %>% as.character()) %>%
  subset(select=c(File,MILC,ENCprime,B,SCUO,MCB,nHE)) %>%
  subset(nHE>=10)

cu250 <- genome_df %>%
  subset(!grepl("Error",File)) %>%
  mutate(MILC=MILC.250 %>% unlist() %>% as.numeric(),
         ENCprime=ENCprime.250 %>% unlist() %>% as.numeric(),
         B=B.250 %>% unlist() %>% as.numeric(),
         SCUO=SCUO.250 %>% unlist() %>% as.numeric(),
         MCB=MCB.250 %>% unlist() %>% as.numeric(),
         nHE=nHE  %>% unlist() %>% as.numeric(),
         File=File  %>% unlist() %>% as.character()) %>%
  subset(select=c(File,MILC,ENCprime,B,SCUO,MCB,nHE)) %>%
  subset(nHE>=10)

names(d)[1] <- "Species"
d <- d %>% as.data.frame(stringsAsFactors=F)

# Merge datasets
rownames(spp_acc) <- spp_acc$V1 %>% gsub(pattern="[.].*",replace="")

cu250$Accession <- cu250$File %>% gsub(pattern="[.].*",replace="")
cu250$Spp <- spp_acc[cu250$Accession,"V2"]
cu250$Species <- lapply(cu250$Spp,rgrep,small_vec=d$Species) %>%
  lapply("[",1) %>% unlist()
cu250$Species[cu250$Spp %in% d$Species] <- cu250$Spp[cu250$Spp %in% d$Species]
cu250 <- merge.easy(cu250,d,key="Species") %>% subset(!is.na(Species))

cu250$Accession <- cu250$File %>% gsub(pattern="[.].*",replace="")
cu250$Spp <- spp_acc[cu250$Accession,"V2"]
cu250$Species <- lapply(cu250$Spp,rgrep,small_vec=d$Species) %>%
  lapply("[",1) %>% unlist()
cu250$Species[cu250$Spp %in% d$Species] <- cu250$Spp[cu250$Spp %in% d$Species]
cu250 <- merge.easy(cu250,d,key="Species") %>% subset(!is.na(Species))

# Average CUB estimates over species
stat_data250 <- cu250 %>%
  subset(Extremophile == FALSE) %>%
  group_by(Species) %>%
  summarise_all(mean,na.rm=T) %>%
  subset(!is.na(Species))

# Average CUB estimates over species, including extremophiles
stat_data_extremo250 <- cu250 %>%
  group_by(Species) %>%
  summarise_all(mean,na.rm=T) %>%
  subset(!is.na(Species))
stat_data_extremo250$OGT <- stat_data_extremo250$OptTemp
stat_data_extremo250$OGT[is.na(stat_data_extremo250$OGT)] <-
  stat_data_extremo250$GrowthTemp[is.na(stat_data_extremo250$OGT)]

# Average CUB estimates over species
stat_data250 <- cu250 %>%
  subset(Extremophile == FALSE) %>%
  group_by(Species) %>%
  summarise_all(mean,na.rm=T) %>%
  subset(!is.na(Species))

# Average CUB estimates over species, including extremophiles
stat_data_extremo250 <- cu250 %>%
  group_by(Species) %>%
  summarise_all(mean,na.rm=T) %>%
  subset(!is.na(Species))
stat_data_extremo250$OGT <- stat_data_extremo250$OptTemp
stat_data_extremo250$OGT[is.na(stat_data_extremo250$OGT)] <-
  stat_data_extremo250$GrowthTemp[is.na(stat_data_extremo250$OGT)]


# Fit Models -------------------------------------------------------------------


model_list150 <- gRodon:::fitTrimModels(stat_data=stat_data150, stat_data_extremo=stat_data_extremo150)
model_list250 <- gRodon:::fitTrimModels(stat_data=stat_data250, stat_data_extremo=stat_data_extremo250)

gRodon_model_base_t150 <- model_list150[[1]]
gRodon_model_temp_t150 <- model_list150[[2]]

gRodon_model_base_t250 <- model_list250[[1]]
gRodon_model_temp_t250 <- model_list250[[2]]

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
     gRodon_model_base_t150,
     gRodon_model_temp_t150,
     gRodon_model_base_t250,
     gRodon_model_temp_t250,
     file="sysdata.rda")
