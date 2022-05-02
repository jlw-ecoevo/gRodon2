

## JLW 2020 - Fit gRodon

# Load Packages ----------------------------------------------------------------

library(dplyr)
library(data.table)
library(MASS)
library(gRodon)
library(parallel)
library(parallelsugar)
library(Biostrings)
library(coRdon)
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

cu <- gRodon:::getStatisticsBatch("~/gRodon2/inst/extdata/euk_genomes/",
                         genetic_code = "1",
                         mc.cores = 6)
cui <- gRodon:::getStatisticsBatch("~/gRodon2/inst/extdata/euk_genomes/",
                         genetic_code = "1",
                         mc.cores = 8,
                         bg = "individual")
setwd("~/gRodon2/inst/extdata/")
save(cu,cui, file = "CodonStatistics_euk.rda")

# Load Growth Dataset ----------------------------------------------------------

# Codon data
load("CodonStatistics_euk.rda")
cu <- cu %>% mutate_all(unlist)
cu$Accession <- cu$File %>% gsub(pattern="[.].*",replace="")
cui <- cui %>% mutate_all(unlist)
cui$Accession <- cui$File %>% gsub(pattern="[.].*",replace="")

### MMETSP Data
setwd("~/gRodon2/inst/extdata/")
# Rose and Caron and Thomas et al.
x1 <- read.csv("euk_growth_ogt_Thomas_MMETSP.csv",stringsAsFactors = F) %>%
  subset(select=c("MMETSP.Strain","Growth.rate","OGT"))
x1$Source <- "Thomas et al."
x2 <- read.csv("euk_growth_ogt_Rose_MMETSP.csv",stringsAsFactors = F)   %>%
  subset(select=c("MMETSP.Strain","Growth.rate","OGT"))
x2$Source <- "Rose and Caron"
y <- rbind(x1,x2) %>% subset(!is.na(MMETSP.Strain))
y$Doubling.Time <- 24*log(2)/y$Growth.rate
y$Growth.rate <- NULL
# Load values matched by literature search
x3 <- read.csv("euk_growth_ogt_Weissman_MMETSP.csv",stringsAsFactors = F)   #%>%
names(x3) <- c("MMETSP.Strain","Doubling.Time","OGT","Source")
x3 <- x3[,c(1,3,4,2)]
growth_rates <- rbind(y,x3)
# Extract species names from strain names
growth_rates$Species.Name <- gsub(";.*","",growth_rates$MMETSP.Strain)
growth_rates$Species.Name[!grepl("sp[.]",growth_rates$Species.Name)] <-
  growth_rates$Species.Name[!grepl("sp[.]",growth_rates$Species.Name)] %>%
  strsplit(.,split=" ") %>%
  lapply(.,"[",1:2) %>%
  lapply(.,paste,collapse=" ") %>%
  unlist() %>%
  gsub(pattern=",",replace="")
# Deduplicate entries (one per species)
d <- growth_rates %>%
  unique() %>%
  group_by(Species.Name) %>%
  dplyr::slice(which.min(Doubling.Time)) %>%
  subset(select=c(Species.Name,Doubling.Time,OGT))
names(d) <- c("Species","d","OGT")
d$MMETSP.Strain <- NULL
#fix some misspelled species names for matching
d$Species[d$Species=="Alexandrium temarense"] <- "Alexandrium tamarense"
d$Species[d$Species=="Heterocapsa triquetra"] <- "Heterocapsa triquestra"
d$Species[d$Species=="Odontella sinensis"] <- "Odontella Sinensis"
#Match accessions
spp_id <- read.csv("mmetsp_spp.csv")
spp_id$Species.Name <- spp_id$Species
spp_id$Species.Name[!grepl("sp[.]",spp_id$Species.Name)] <-
  spp_id$Species.Name[!grepl("sp[.]",spp_id$Species.Name)] %>%
  strsplit(.,split=" ") %>%
  lapply(.,"[",1:2) %>%
  lapply(.,paste,collapse=" ") %>%
  unlist() %>%
  gsub(pattern=",",replace="")
spp_id$Species <- spp_id$Species.Name
spp_id$Species.Name <- NULL
exclude <- readLines("eukprey_list.txt")
mmetsp_d <- merge.easy(d,spp_id,key="Species") %>%
  subset(!Species %in% exclude) %>%
  subset(!is.na(OGT))
#Missing accession
mmetsp_d$Accession[mmetsp_d$Species=="Pavlova lutheri"] <- "METSP1463"

### NCBI Data
x <- read.csv("euk_growth_ogt_Weissman_RefSeq.csv",stringsAsFactors = F)
ncbi_d <- data.frame(Species=character(),
                           de=numeric(),
                           OGT=numeric(),
                           Accession=character())
for(i in 1:nrow(x)){
  ncbi_d <- rbind(ncbi_d,
                        data.frame(Species=x[i,"Species.Name"],
                                   d=x[i,"Doubling.Time"],
                                   OGT=x[i,"OGT"],
                                   Accession=x[i,"Accession"] %>%
                                     strsplit(split=";") %>%
                                     unlist() %>%
                                     gsub(pattern="[.].*",replace="")))
}

### merge and average over species ---------------------------------------------

mmetsp_da <- merge.easy(mmetsp_d,cu,key="Accession") %>%
  subset(nHE>10) %>%
  group_by(Species) %>%
  summarise_all(mean,na.rm=T) %>%
  subset(!is.na(Species)) %>%
  subset(!is.na(CUBHE))
ncbi_da <- merge.easy(ncbi_d,cu,key="Accession") %>%
  subset(nHE>10) %>%
  group_by(Species) %>%
  summarise_all(mean,na.rm=T) %>%
  subset(!is.na(Species)) %>%
  subset(!is.na(CUBHE))
stat_data <- rbind(mmetsp_da,ncbi_da)

mmetsp_di <- merge.easy(mmetsp_d,cui,key="Accession") %>%
  subset(nHE>10) %>%
  group_by(Species) %>%
  summarise_all(mean,na.rm=T) %>%
  subset(!is.na(Species)) %>%
  subset(!is.na(CUBHE))
ncbi_di <- merge.easy(ncbi_d,cui,key="Accession") %>%
  subset(nHE>10) %>%
  group_by(Species) %>%
  summarise_all(mean,na.rm=T) %>%
  subset(!is.na(Species)) %>%
  subset(!is.na(CUBHE))
stat_data_i <- rbind(mmetsp_di,ncbi_di)


# Test for outliers
test_outliers <- EnvStats::rosnerTest(stat_data$d %>% log10())
outlier_ind <- test_outliers$all.stats$Obs.Num[test_outliers$all.stats$Outlier==T]
stat_data <- stat_data[-outlier_ind,]
stat_data_i <- stat_data_i[-outlier_ind,]

# Fit Models -------------------------------------------------------------------

model_list <- fitGCModels(stat_data, stat_data)

gRodon_model_base_euk <- model_list[[3]]
gRodon_model_temp_euk <- model_list[[4]]
lambda_milc_euk <- model_list[[5]]

model_list <- fitGCModels(stat_data_i, stat_data_i)

gRodon_model_base_euk_i <- model_list[[3]]
gRodon_model_temp_euk_i <- model_list[[4]]
lambda_milc_euk_i <- model_list[[5]]

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
     gRodon_model_base_euk,
     gRodon_model_temp_euk,
     lambda_milc_euk,
     gRodon_model_base_euk_i,
     gRodon_model_temp_euk_i,
     lambda_milc_euk_i,
     file="sysdata.rda")
