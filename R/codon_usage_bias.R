

getCUB <- function(fna_tab, highly_expressed, method = "MILC"){

  if(method == "MILC"){
    # MILC estimate of codon usage bias from coRdon
    #   - Bias of highly expessed genes, with genome-wide codon usage considered
    #     expected codon usage
    x <- MILC(fna_tab,id_or_name2 = "11")
    x[highly_expressed, 1] %>% median() %>% return()

  } else if(method == "MILCgenomic"){
    # MILC estimate of codon usage bias from coRdon
    #   - Bias of all genes, with genome-wide codon usage considered
    #     expected codon usage
    x <- MILC(fna_tab,id_or_name2 = "11")
    x[, 1] %>% median() %>% return()

  } else if(method == "consistency"){
    # Consistency of CUB across highly expressed genes
    #   - Bias of highly expessed genes, with codon usage of highly expressed
    #     gene set considered expected codon usage
    milc <- MILC(fna_tab,
                 subsets = list(HE = highly_expressed),
                 id_or_name2 = "11")
    milc[highly_expressed, 2] %>% mean() %>% return()

  } else {
    stop("Error: Please pick an implemented method (MILC, MILCgenomic, consistency)")
  }
}

getWeightedCUB <- function(fna_tab, highly_expressed, depth_of_coverage){
    x <- MILC(fna_tab,id_or_name2 = "11")
    x[highly_expressed, 1] %>%
      weightedMedian(., w = depth_of_coverage[highly_expressed]) %>%
      return()
}


getCodonStatistics <- function(genes, highly_expressed, fragments = FALSE, depth_of_coverage = NULL){

  if(sum(highly_expressed) == 0){
    stop("No highly expressed genes?")
  }

  #Remove short sequences and sequences that are not multiples of 3
  if(fragments == TRUE){
    filtered <- filterSeq(genes = genes,
                          highly_expressed = highly_expressed,
                          length_threshold = 120,
                          depth_of_coverage = depth_of_coverage)
  } else {
    filtered <- filterSeq(genes = genes,
                          highly_expressed = highly_expressed,
                          depth_of_coverage = depth_of_coverage)
  }
  genes <- filtered$Genes
  highly_expressed <- filtered$HE
  depth_of_coverage <- filtered$Depth

  # codon table
  codon_table <- codonTable(genes)

  # codon pair counts
  codon_pair_table <- getPairCounts(genes)

  if(!is.null(depth_of_coverage)){
    return(data.frame(CUBHE = getWeightedCUB(codon_table,
                                     highly_expressed,
                                     depth_of_coverage),
                      ConsistencyHE = NA,
                      CPB = NA,
                      FilteredSequences = filtered$Filtered,
                      stringsAsFactors = FALSE))
  } else {
    return(data.frame(CUBHE = getCUB(codon_table,
                                     highly_expressed,
                                     method = "MILC"),
                      ConsistencyHE = getCUB(codon_table,
                                             highly_expressed,
                                             method = "consistency"),
                      CPB = getCPB(codon_pair_table),
                      FilteredSequences = filtered$Filtered,
                      stringsAsFactors = FALSE))
  }
}
