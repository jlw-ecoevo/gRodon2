
#' Predict Growth Rate
#'
#' This function predicts the growth rate of prokaryotic or eukaryotic organisms based on
#' genomic codon usage patterns (Weissman et. al. 2021; Weissman et al. TBD).
#'
#' @param genes DNAStringSet object holding all in-frame coding sequences (nucleotide sequences)
#'  from a genome. See Biostrings package.
#' @param highly_expressed Logical vector describing the set of highly expressed
#' genes. Must be of same length as \code{genes}. Typically these are ribosomal proteins
#' (all models were trained using ribosomal proteins as the highly expressed set.)
#' @param mode Whether to run prediction in full, partial, metagenome_v1, metagenome_v2,
#' eukaryote, or metagenome_euk mode
#' (by default gRodon applies the full model). Mode metagenome_v2 may run slower
#' than the other prediction modes when consistency>0.6.
#' @param temperature Optimal growth temperature. By default this is set as
#' "none" for prokaryotes and we do not guarantee good results for non-mesophilic
#'  since prokaryotes few were used to fit the model. For eukaryotes, though, including
#'  an optimal growth temperature drastically improves model predictions as there is
#'  a very strong relationship between OGT and max. growth rate for these organisms.
#' @param training_set Whether to use models trained on the original Vieira-Silva et al.
#' doubling time dataset or doubling times drawn from the Madin et al. database. This
#' setting is only used for prokaryotic modes (eukaryotic models based on their own
#'  training set from Weissman et al. TBD). For metagenome_v2 mode, only the madin
#'  set is available.By default training set is now set to Madin. For AOA and NOB try AOA_NOB
#'  which includes an expanded set of these organisms (including some measurements from enrichment cultures)
#' @param depth_of_coverage When using metagenome mode, provide a vector containing
#' the coverage of your ORFs to improve your estimate
#' @param fragments Do not change (performance will suffer, publication forthcoming).
#' If using gene fragments predicted from reads, will use a
#' more permissive length filter (120bp as opposed to 240bp)
#' @param genetic_code The genetic code of the organism to be used in codon usage
#' calculations (see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). By
#' default code "1" is used for eukaryote mode and code "11" for prokaryotic modes.
#' @param n_le The number of genes to sample to generate a background estimate of
#' codon usage against which to calculate the codon usage of the highly expressed
#' genes. Only used by metagenome_v2 mode and only when consistency>0.6. Increasing
#' will make prediction slower, but may decrease variance. Decreasing is not recommended
#' but will increase speed.
#' @param bg Feature in testing, do not change
#' @return gRodon returns a list with the following elements:
#' \describe{
#'   \item{CUBHE}{Median codon usage bias of the highly expressed genes (MILC)
#'   calculated using the genome-wide codon usage as the expected bias}
#'   \item{ConsistencyHE}{Mean codon usage bias of the highly expressed genes
#'   (MILC) calculated using the codon usage of highly expressed genes as the
#'   expected bias}
#'   \item{CPB}{Genome-wide codon pair bias (Coleman et al. 2008)}
#'   \item{FilteredSequences}{Number of gene sequences filtered out during
#'   calulation (due to length and/or presence of ambiguous bases)}
#'   \item{d}{Predicted doubling time in hours}
#'   \item{LowerCI}{Lower CI of \code{d} (2.5%) from linear model}
#'   \item{UpperCI}{Upper CI of \code{d} (97.5%) from linear model}
#' }
#' @examples
#' # Load in example genome (Streptococcus pyogenes M1, downloaded from RefSeq)
#' # included with gRodon
#' library(Biostrings)
#' path_to_genome <- system.file('extdata',
#'   'GCF_000349925.2_ASM34992v2_cds_from_genomic.fna.gz',
#'   package = 'gRodon')
#' genes <- readDNAStringSet(path_to_genome)
#'
#' # Search pre-existing annotations for ribosomal proteins, which we
#' # will use as our set of highly expressed genes
#' highly_expressed <- grepl("^(?!.*(methyl|hydroxy)).*0S ribosomal protein",names(genes),ignore.case = T, perl = TRUE)
#'
#' # Run the gRodon growth prediction pipeline
#' predictGrowth(genes, highly_expressed)
#'
#' # Run gRodon with temperature option (not needed for mesophiles, gRodon not
#' # validated on extremophiles, use with care)
#' predictGrowth(genes, highly_expressed, temperature = 37)
#'
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom dplyr group_by mutate count summarise
#' @importFrom matrixStats weightedMedian
#' @import coRdon
#' @import Biostrings
predictGrowth <- function(genes,
                          highly_expressed,
                          mode = "full",
                          temperature = "none",
                          training_set = "madin",
                          depth_of_coverage = NULL,
                          fragments = FALSE,
                          genetic_code = NULL,
                          n_le = 100,
                          bg = "all"){

  if(! mode %in% c("full","partial","metagenome_v1","metagenome_v2","metagenome_euk","eukaryote","meta_testing","meta_nogc_testing")){
    stop("Invalid mode. Please pick an available prediction mode (\"full\", \"partial\", \"metagenome_v1\", \"metagenome_v2\", \"eukaryote\")")
  }

  if((! training_set  %in% c("vs","madin","AOA_NOB")) & !mode %in% c("eukaryote","metagenome_euk")){
    stop("Invalid training set. Please pick an available model (\"vs\", \"madin\", \"AOA_NOB\")")
  }

  if(sum(highly_expressed)<10){
    warning("Less than 10 highly expressed genes provided, performance may suffer")
  }

  if(!(mode %in% c("metagenome_v1","metagenome_v2","meta_testing","meta_nogc_testing","metagenome_euk")) & !is.null(depth_of_coverage)){
    warning("Ignoring depth_of_coverage because not in metagenome mode")
    depth_of_coverage <- NULL
  }

  if(mode %in% c("metagenome_v1","metagenome_v2","metagenome_euk") & is.null(depth_of_coverage)){
    warning("Provide depth_of_coverage for your ORFs for a more realistic average community growth rate")
  }

  if(is.null(genetic_code) & mode %in% c("eukaryote","metagenome_euk")){
    genetic_code <- "1"
  } else if(is.null(genetic_code)){
    genetic_code <- "11"
  } else {
    genetic_code <- as.character(genetic_code)
    warning("Models were trained with a default genetic code of '1' for eukaryotes and '11' for prokaryotes. See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi")
  }

  if(is.na(temperature) | is.null(temperature)){
    temperature="none"
    warning("Missing temperature. Flag set to `temperature=\"none\"`.")
  }

  if(temperature=="none" & mode %in% c("eukaryote","metagenome_euk")){
    warning("For best results for eukaryotes an optimal growth temperature must be provided. Much of the variation in max. growth rate between species is explained by temperature, independent of any diffferences in codon usage.")
  }

  if(mode=="metagenome_v2" & training_set %in% c("vs","AOA_NOB")){
    training_set <- "madin"
    warning("Training set automatically set to \"madin\" for metagenome_v2 mode")
  }

  i_flag <- 0

  # Calculate codon data
  if(bg=="all"){
    codon_stats <- getCodonStatistics(genes = genes,
                                      highly_expressed = highly_expressed,
                                      fragments = fragments,
                                      depth_of_coverage = depth_of_coverage,
                                      genetic_code = genetic_code)

    codon_stats$dCUB <- (codon_stats$CUB-codon_stats$CUBHE)/codon_stats$CUB
    # Predict growth rate (stored models - sysdata.rda)
    if(training_set=="vs" & mode!="eukaryote" & mode!="meta_testing" & mode!="meta_nogc_testing"){
      if(temperature == "none" & mode=="full"){
        pred <- stats::predict.lm(gRodon_model_base,
                                  newdata = codon_stats,
                                  interval = "confidence")

      } else if(temperature == "none" & mode=="metagenome_v1"){
        pred <- stats::predict.lm(gRodon_model_meta,
                                  newdata = codon_stats,
                                  interval = "confidence")

      } else if(temperature == "none" & mode=="partial"){
        pred <- stats::predict.lm(gRodon_model_partial,
                                  newdata = codon_stats,
                                  interval = "confidence")

      } else if(temperature != "none" & mode=="full"){
        codon_stats$OGT <- temperature
        pred <- stats::predict.lm(gRodon_model_temp,
                                  newdata = codon_stats,
                                  interval = "confidence")

      } else if(temperature != "none" & mode=="metagenome_v1"){
        codon_stats$OGT <- temperature
        pred <- stats::predict.lm(gRodon_model_meta_temp,
                                  newdata = codon_stats,
                                  interval = "confidence")

      } else if(temperature != "none" & mode=="partial"){
        codon_stats$OGT <- temperature
        pred <- stats::predict.lm(gRodon_model_partial_temp,
                                  newdata = codon_stats,
                                  interval = "confidence")
      }
      #Transform back from box-cox
      pred_back_transformed <- boxcoxTransform(pred,
                                               lambda_milc,
                                               back_transform = TRUE)
    } else if(training_set=="AOA_NOB" & mode!="eukaryote" & mode!="meta_testing" & mode!="meta_nogc_testing"){
      if(temperature == "none" & mode=="full"){
        pred <- stats::predict.lm(gRodon_model_base_AOANOB,
                                  newdata = codon_stats,
                                  interval = "confidence")

      } else if(temperature == "none" & mode=="metagenome_v1"){
        pred <- stats::predict.lm(gRodon_model_meta_AOANOB,
                                  newdata = codon_stats,
                                  interval = "confidence")

      } else if(temperature == "none" & mode=="partial"){
        pred <- stats::predict.lm(gRodon_model_partial_AOANOB,
                                  newdata = codon_stats,
                                  interval = "confidence")

      } else if(temperature != "none" & mode=="full"){
        codon_stats$OGT <- temperature
        pred <- stats::predict.lm(gRodon_model_temp_AOANOB,
                                  newdata = codon_stats,
                                  interval = "confidence")

      } else if(temperature != "none" & mode=="metagenome_v1"){
        codon_stats$OGT <- temperature
        pred <- stats::predict.lm(gRodon_model_meta_temp_AOANOB,
                                  newdata = codon_stats,
                                  interval = "confidence")

      } else if(temperature != "none" & mode=="partial"){
        codon_stats$OGT <- temperature
        pred <- stats::predict.lm(gRodon_model_partial_temp_AOANOB,
                                  newdata = codon_stats,
                                  interval = "confidence")
      }
      #Transform back from box-cox
      pred_back_transformed <- boxcoxTransform(pred,
                                               lambda_milc_AOANOB,
                                               back_transform = TRUE)
    } else if(mode!="eukaryote" & mode!="meta_testing" & mode!="meta_nogc_testing"){
      if(temperature == "none" & mode=="full"){
        pred <- stats::predict.lm(gRodon_model_base_madin,
                                  newdata = codon_stats,
                                  interval = "confidence")

      } else if(temperature == "none" & mode=="metagenome_v1"){
        pred <- stats::predict.lm(gRodon_model_meta_madin,
                                  newdata = codon_stats,
                                  interval = "confidence")

      } else if(temperature == "none" & mode=="metagenome_v2"){
        if(codon_stats$ConsistencyHE<0.6){
          pred <- stats::predict.lm(gRodon_model_meta_madin,
                                    newdata = codon_stats,
                                    interval = "confidence")
        } else {
          codon_stats <- getCodonStatistics_i(genes = genes,
                                              highly_expressed = highly_expressed,
                                              fragments = fragments,
                                              depth_of_coverage = depth_of_coverage,
                                              genetic_code = genetic_code,
                                              n_le = n_le)
          codon_stats$dCUB <- (codon_stats$CUB-codon_stats$CUBHE)/codon_stats$CUB
          pred <- stats::predict.lm(gRodon_model_newmeta_i,
                                    newdata = codon_stats,
                                    interval = "confidence")
          i_flag <- 1
        }
      } else if(temperature == "none" & mode=="metagenome_euk"){
        if(codon_stats$ConsistencyHE<0.6){
          pred <- stats::predict.lm(gRodon_model_base_euk,
                                    newdata = codon_stats,
                                    interval = "confidence")
        } else {
          codon_stats <- getCodonStatistics_i(genes = genes,
                                              highly_expressed = highly_expressed,
                                              fragments = fragments,
                                              depth_of_coverage = depth_of_coverage,
                                              genetic_code = genetic_code,
                                              n_le = n_le)
          codon_stats$dCUB <- (codon_stats$CUB-codon_stats$CUBHE)/codon_stats$CUB
          pred <- stats::predict.lm(gRodon_model_base_euk_i,
                                    newdata = codon_stats,
                                    interval = "confidence")
          i_flag <- 1
        }

      } else if(temperature == "none" & mode=="partial"){
        pred <- stats::predict.lm(gRodon_model_partial_madin,
                                  newdata = codon_stats,
                                  interval = "confidence")

      } else if(temperature != "none" & mode=="full"){
        codon_stats$OGT <- temperature
        pred <- stats::predict.lm(gRodon_model_temp_madin,
                                  newdata = codon_stats,
                                  interval = "confidence")

      } else if(temperature != "none" & mode=="metagenome_v1"){
        codon_stats$OGT <- temperature
        pred <- stats::predict.lm(gRodon_model_meta_temp_madin,
                                  newdata = codon_stats,
                                  interval = "confidence")

      } else if(temperature != "none" & mode=="metagenome_v2"){
        if(codon_stats$ConsistencyHE<0.6){
          codon_stats$OGT <- temperature
          pred <- stats::predict.lm(gRodon_model_meta_temp_madin,
                                    newdata = codon_stats,
                                    interval = "confidence")
        } else {
          codon_stats <- getCodonStatistics_i(genes = genes,
                                              highly_expressed = highly_expressed,
                                              fragments = fragments,
                                              depth_of_coverage = depth_of_coverage,
                                              genetic_code = genetic_code,
                                              n_le = n_le)
          codon_stats$dCUB <- (codon_stats$CUB-codon_stats$CUBHE)/codon_stats$CUB
          codon_stats$OGT <- temperature
          pred <- stats::predict.lm(gRodon_model_newmeta_temp_i,
                                    newdata = codon_stats,
                                    interval = "confidence")
          i_flag <- 1
        }
      } else if(temperature != "none" & mode=="metagenome_euk"){
        if(codon_stats$ConsistencyHE<0.6){
          codon_stats$OGT <- temperature
          pred <- stats::predict.lm(gRodon_model_temp_euk,
                                    newdata = codon_stats,
                                    interval = "confidence")
        } else {
          codon_stats <- getCodonStatistics_i(genes = genes,
                                              highly_expressed = highly_expressed,
                                              fragments = fragments,
                                              depth_of_coverage = depth_of_coverage,
                                              genetic_code = genetic_code,
                                              n_le = n_le)
          codon_stats$dCUB <- (codon_stats$CUB-codon_stats$CUBHE)/codon_stats$CUB
          codon_stats$OGT <- temperature
          pred <- stats::predict.lm(gRodon_model_temp_euk_i,
                                    newdata = codon_stats,
                                    interval = "confidence")
          i_flag <- 1
        }
      } else if(temperature != "none" & mode=="partial"){
        codon_stats$OGT <- temperature
        pred <- stats::predict.lm(gRodon_model_partial_temp_madin,
                                  newdata = codon_stats,
                                  interval = "confidence")
      }
      #Transform back from box-cox
      if(i_flag==0){
        pred_back_transformed <- boxcoxTransform(pred,
                                                 lambda_milc_madin,
                                                 back_transform = TRUE)
      } else {
        pred_back_transformed <- boxcoxTransform(pred,
                                                 lambda_newmeta_i,
                                                 back_transform = TRUE)
      }

    } else if(mode=="eukaryote"){
      codon_stats$dCUB <- (codon_stats$CUB-codon_stats$CUBHE)/codon_stats$CUB
      if(temperature == "none"){
        pred <- stats::predict.lm(gRodon_model_base_euk,
                                  newdata = codon_stats,
                                  interval = "confidence")
      } else {
        codon_stats$OGT <- temperature
        pred <- stats::predict.lm(gRodon_model_temp_euk,
                                  newdata = codon_stats,
                                  interval = "confidence")
      }
      #Transform back from box-cox
      pred_back_transformed <- boxcoxTransform(pred,
                                               lambda_milc_euk,
                                               back_transform = TRUE)
    } else if(mode=="meta_testing"){
      if(temperature == "none"){
        pred <- stats::predict.lm(gRodon_model_newmeta,
                                  newdata = codon_stats,
                                  interval = "confidence")
      } else {
        codon_stats$OGT <- temperature
        pred <- stats::predict.lm(gRodon_model_newmeta_temp,
                                  newdata = codon_stats,
                                  interval = "confidence")
      }
      #Transform back from box-cox
      pred_back_transformed <- boxcoxTransform(pred,
                                               lambda_newmeta,
                                               back_transform = TRUE)
    } else if(mode=="meta_nogc_testing"){
      if(temperature == "none"){
        pred <- stats::predict.lm(gRodon_model_newmeta_nogc,
                                  newdata = codon_stats,
                                  interval = "confidence")
      } else {
        codon_stats$OGT <- temperature
        pred <- stats::predict.lm(gRodon_model_newmeta_nogc_temp,
                                  newdata = codon_stats,
                                  interval = "confidence")
      }
      #Transform back from box-cox
      pred_back_transformed <- boxcoxTransform(pred,
                                               lambda_newmeta,
                                               back_transform = TRUE)
    }


  } else if(bg=="individual"){
    if(!(mode %in% c("meta_testing","meta_nogc_testing","metagenome","eukaryote"))){
      stop("Mode not compatible with gene-level CUB calculations")
    }
    codon_stats <- getCodonStatistics_i(genes = genes,
                                        highly_expressed = highly_expressed,
                                        fragments = fragments,
                                        depth_of_coverage = depth_of_coverage,
                                        genetic_code = genetic_code)

    codon_stats$dCUB <- (codon_stats$CUB-codon_stats$CUBHE)/codon_stats$CUB

    if(mode=="meta_testing"){
      if(temperature == "none"){
        pred <- stats::predict.lm(gRodon_model_newmeta_i,
                                  newdata = codon_stats,
                                  interval = "confidence")
      } else {
        codon_stats$OGT <- temperature
        pred <- stats::predict.lm(gRodon_model_newmeta_temp_i,
                                  newdata = codon_stats,
                                  interval = "confidence")
      }
      #Transform back from box-cox
      pred_back_transformed <- boxcoxTransform(pred,
                                               lambda_newmeta_i,
                                               back_transform = TRUE)
    } else if(mode=="meta_nogc_testing"){
      if(temperature == "none"){
        pred <- stats::predict.lm(gRodon_model_newmeta_nogc_i,
                                  newdata = codon_stats,
                                  interval = "confidence")
      } else {
        codon_stats$OGT <- temperature
        pred <- stats::predict.lm(gRodon_model_newmeta_nogc_temp_i,
                                  newdata = codon_stats,
                                  interval = "confidence")
      }
      #Transform back from box-cox
      pred_back_transformed <- boxcoxTransform(pred,
                                               lambda_newmeta_i,
                                               back_transform = TRUE)
    } else if(mode=="eukaryote"){
      if(temperature == "none"){
        pred <- stats::predict.lm(gRodon_model_base_euk_i,
                                  newdata = codon_stats,
                                  interval = "confidence")
      } else {
        codon_stats$OGT <- temperature
        pred <- stats::predict.lm(gRodon_model_temp_euk_i,
                                  newdata = codon_stats,
                                  interval = "confidence")
      }
      #Transform back from box-cox
      pred_back_transformed <- boxcoxTransform(pred,
                                               lambda_milc_euk_i,
                                               back_transform = TRUE)
    } else if(mode=="metagenome_v1"){
      if(temperature == "none"){
        pred <- stats::predict.lm(gRodon_model_meta_madin_i,
                                  newdata = codon_stats,
                                  interval = "confidence")
      } else {
        codon_stats$OGT <- temperature
        pred <- stats::predict.lm(gRodon_model_meta_temp_madin_i,
                                  newdata = codon_stats,
                                  interval = "confidence")
      }
      #Transform back from box-cox
      pred_back_transformed <- boxcoxTransform(pred,
                                               lambda_milc_madin_i,
                                               back_transform = TRUE)
    }
  } else{
    stop("Feature in testing, please set bg==\"all\" for normal gRodon behavior")
  }



  #attach prediction
  codon_stats$d <- pred_back_transformed[,"fit"]
  codon_stats$LowerCI <- pred_back_transformed[,"lwr"]
  codon_stats$UpperCI <- pred_back_transformed[,"upr"]

  #Return
  if(is.na(pred_back_transformed[,"fit"]) & pred[,"fit"]>6){
    warning("Estimated doubling time very long. Essentially goes to infinity (gives NA value after back-transforming from box-cox).")
  } else if(pred_back_transformed[,"fit"]>5 & mode!="eukaryote"){
    warning("Estimated doubling time >5 hours. CUB signal saturates at approx. 5 hrs... gRodon may underestimate doubling times above this range. Consider simply reporting as '>5hrs'. (In other words, this microbe definitely grows slowly, but we can't tell you quite how slowly).")
  } else if(pred_back_transformed[,"fit"]>40 & mode=="eukaryote"){
    warning("Estimated doubling time >40 hours. CUB signal saturates at approx. 40 hrs for eukaryotes... gRodon may underestimate doubling times above this range. Consider simply reporting as '>40hrs'. (In other words, this microbe definitely grows slowly, but we can't tell you quite how slowly).")
  }
  return(as.list(codon_stats))
}
