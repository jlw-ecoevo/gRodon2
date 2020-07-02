

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

cu <- gRodon:::getStatisticsBatch("~/gRodon/inst/extdata/vs_genomes", 
                         mc.cores = 4)
setwd("~/gRodon/inst/extdata/")
save("CodonStatistics.rda")

# Load Growth Dataset ----------------------------------------------------------

setwd("~/gRodon/inst/extdata/")
load("GrowthRates.rda")
load("CodonStatistics.rda")
load("Accession2Species.rda")

# Merge datasets
rownames(spp_acc) <- spp_acc$V1 %>% gsub(pattern="[.].*",replace="")
cs_df$Accession <- cs_df$File %>% gsub(pattern="[.].*",replace="")
cs_df$Species <- spp_acc[cs_df$Accession,"V2"]
cs_df$Species <- lapply(cs_df$Species,rgrep,small_vec=d$Species) %>% 
  lapply("[",1) %>% unlist()
cs_df <- merge.easy(cs_df,d,key="Species")

# Average CUB estimates over species
stat_data <- cs_df %>% group_by(Species) %>% 
  summarise_all(mean,na.rm=T) %>% 
  subset(!is.na(Species)) 

# Fit Models -------------------------------------------------------------------

model_list <- fitModels(stat_data)

gRodon_model_base <- model_list[[1]]
gRodon_model_temp <- model_list[[2]]
gRodon_model_partial <- model_list[[3]]
gRodon_model_partial_temp <- model_list[[4]]
gRodon_model_meta <- model_list[[5]]
gRodon_model_meta_temp <- model_list[[6]]
lambda_milc <- model_list[[7]]

setwd("~/gRodon/R/")
save(gRodon_model_base,
     gRodon_model_temp,
     gRodon_model_partial,
     gRodon_model_partial_temp,
     gRodon_model_meta,
     gRodon_model_meta_temp,
     lambda_milc,
     file="sysdata.rda")