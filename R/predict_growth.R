
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
#' @param mode Whether to run prediction in full, partial, metagenome, or eukaryote mode
#' (by default gRodon applies the full model)
#' @param temperature Optimal growth temperature. By default this is set as
#' "none" for prokaryotes and we do not guarantee good results for non-mesophilic
#'  since prokaryotes few were used to fit the model. For eukaryotes, though, including
#'  an optimal growth temperature drastically improves model predictions as there is
#'  a very strong relationship between OGT and max. growth rate for these organisms.
#' @param training_set Whether to use models trained on the original Vieira-Silva et al.
#' doubling time dataset or doubling times drawn from the Madin et al. database. This
#' setting is only used for prokaryotic modes (eukaryotic models based on their own
#'  training set from Weissman et al. TBD)
#' @param depth_of_coverage When using metagenome mode, provide a vector containing
#' the coverage of your ORFs to improve your estimate
#' @param fragments If using gene fragments predicted from reads, will use a
#' more permissive length filter (120bp as opposed to 240bp)
#' @param genetic_code The genetic code of the organism to be used in codon usage
#' calculations (see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). By
#' default code "1" is used for eukaryote mode and code "11" for prokaryotic modes.
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
#' path_to_genome <- system.file('extdata',
#'   'GCF_000349925.2_ASM34992v2_cds_from_genomic.fna.gz',
#'   package = 'gRodon2')
#' genes <- readDNAStringSet(path_to_genome)
#'
#' # Search pre-existing annotations for ribosomal proteins, which we
#' # will use as our set of highly expressed genes
#' highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case = T)
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
                          training_set = "vs",
                          depth_of_coverage = NULL,
                          fragments = FALSE,
                          genetic_code = NULL){

  if(! mode %in% c("full","partial","metagenome","eukaryote","meta_testing","meta_nogc_testing")){
    stop("Invalid mode. Please pick an available prediction mode (\"full\", \"partial\", \"metagenome\", \"eukaryote\")")
  }

  if((! training_set  %in% c("vs","madin")) & mode!="eukaryote"){
    stop("Invalid training set. Please pick an available model (\"vs\", \"madin\")")
  }

  if(sum(highly_expressed)<10){
    warning("Less than 10 highly expressed genes provided, performance may suffer")
  }

  if((mode!="metagenome" | mode!="meta_testing" | mode!="meta_nogc_testing") & !is.null(depth_of_coverage)){
    warning("Ignoring depth_of_coverage because not in metagenome mode")
    depth_of_coverage <- NULL
  }

  if(mode=="metagenome" & is.null(depth_of_coverage)){
    warning("Provide depth_of_coverage for your ORFs for a more realistic average community growth rate")
  }

  if(is.null(genetic_code) & mode=="eukaryote"){
    genetic_code <- "1"
  } else if(is.null(genetic_code)){
    genetic_code <- "11"
  } else {
    genetic_code <- as.character(genetic_code)
    warning("Models were trained with a default genetic code of '1' for eukaryotes and '11' for prokaryotes. See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi")
  }

  if(temperature=="none" & mode=="eukaryote"){
    warning("For best results for eukaryotes an optimal growth temperature must be provided. Much of the variation in max. growth rate between species is explained by temperature, independent of any diffferences in codon usage.")
  }

  # Calculate codon data
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

    } else if(temperature == "none" & mode=="metagenome"){
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

    } else if(temperature != "none" & mode=="metagenome"){
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
  } else if(mode!="eukaryote" & mode!="meta_testing" & mode!="meta_nogc_testing"){
    if(temperature == "none" & mode=="full"){
      pred <- stats::predict.lm(gRodon_model_base_madin,
                                newdata = codon_stats,
                                interval = "confidence")

    } else if(temperature == "none" & mode=="metagenome"){
      pred <- stats::predict.lm(gRodon_model_meta_madin,
                                newdata = codon_stats,
                                interval = "confidence")

    } else if(temperature == "none" & mode=="partial"){
      pred <- stats::predict.lm(gRodon_model_partial_madin,
                                newdata = codon_stats,
                                interval = "confidence")

    } else if(temperature != "none" & mode=="full"){
      codon_stats$OGT <- temperature
      pred <- stats::predict.lm(gRodon_model_temp_madin,
                                newdata = codon_stats,
                                interval = "confidence")

    } else if(temperature != "none" & mode=="metagenome"){
      codon_stats$OGT <- temperature
      pred <- stats::predict.lm(gRodon_model_meta_temp_madin,
                                newdata = codon_stats,
                                interval = "confidence")

    } else if(temperature != "none" & mode=="partial"){
      codon_stats$OGT <- temperature
      pred <- stats::predict.lm(gRodon_model_partial_temp_madin,
                                newdata = codon_stats,
                                interval = "confidence")
    }
    #Transform back from box-cox
    pred_back_transformed <- boxcoxTransform(pred,
                                             lambda_milc_madin,
                                             back_transform = TRUE)
  } else if(mode=="eukaryote"){
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
