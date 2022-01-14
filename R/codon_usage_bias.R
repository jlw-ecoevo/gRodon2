


getCUB <- function(fna_tab, highly_expressed, method = "MILC", genetic_code = "11"){

  if(method == "MILC"){
    # MILC estimate of codon usage bias from coRdon
    #   - Bias of highly expessed genes, with genome-wide codon usage considered
    #     expected codon usage
    x <- MILC(fna_tab,id_or_name2 = genetic_code)
    x[highly_expressed, 1] %>% median() %>% return()

  } else if(method == "MILCgenomic"){
    # MILC estimate of codon usage bias from coRdon
    #   - Bias of all genes, with genome-wide codon usage considered
    #     expected codon usage
    x <- MILC(fna_tab,id_or_name2 = genetic_code)
    x[, 1] %>% median() %>% return()

  } else if(method == "MILC_i"){
    # MILC estimate of codon usage bias from coRdon
    #   - Bias of all genes, with genome-wide codon usage considered
    #     expected codon usage
    y <- numeric()
    for(i in 1:sum(highly_expressed)){
      y <- c(y,MILC(fna_tab[highly_expressed][i]))
    }
    y <- MILC(fna_tab,id_or_name2 = genetic_code)
    y[highly_expressed, 1] %>% mean() %>% return()

  } else if(method == "MILCgenomic_i"){
    # MILC estimate of codon usage bias from coRdon
    #   - Bias of all genes, with genome-wide codon usage considered
    #     expected codon usage
    y <- numeric()
    for(i in 1:length(fna_tab)){
      y <- c(y,MILC(fna_tab[i]))
    }
    y <- MILC(fna_tab,id_or_name2 = genetic_code)
    y[, 1] %>% mean() %>% return()

  } else if(method == "MILCgenomic"){
    # MILC estimate of codon usage bias from coRdon
    #   - Bias of all genes, with genome-wide codon usage considered
    #     expected codon usage
    x <- MILC(fna_tab,id_or_name2 = genetic_code)
    x[, 1] %>% median() %>% return()


  } else if(method == "consistency"){
    # Consistency of CUB across highly expressed genes
    #   - Bias of highly expessed genes, with codon usage of highly expressed
    #     gene set considered expected codon usage
    milc <- MILC(fna_tab,
                 subsets = list(HE = highly_expressed),
                 id_or_name2 = genetic_code)
    milc[highly_expressed, 2] %>% mean() %>% return()

  } else {
    stop("Error: Please pick an implemented method (MILC, MILCgenomic, consistency)")
  }
}

getWeightedCUB <- function(fna_tab, highly_expressed, depth_of_coverage, genetic_code = "11"){
    x <- MILC(fna_tab,id_or_name2 = genetic_code)
    x[, 1] %>%
      weightedMedian(., w = depth_of_coverage) %>%
      return()
}

getWeightedCUBHE <- function(fna_tab, highly_expressed, depth_of_coverage, genetic_code = "11"){
  x <- MILC(fna_tab,id_or_name2 = genetic_code)
  x[highly_expressed, 1] %>%
    weightedMedian(., w = depth_of_coverage[highly_expressed]) %>%
    return()
}

getWeightedConsistency <- function(fna_tab, highly_expressed, depth_of_coverage, genetic_code = "11"){
  x <- MILC(fna_tab,
            subsets = list(HE = highly_expressed),
            id_or_name2 = genetic_code)
  sum(x[highly_expressed, 2]*(depth_of_coverage[highly_expressed]/sum(depth_of_coverage[highly_expressed]))) %>%
    return()
}

getCodonStatistics <- function(genes, highly_expressed, fragments = FALSE, depth_of_coverage = NULL, genetic_code = "11"){

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
  codon_pair_table <- getPairCounts(genes, genetic_code = genetic_code)

  if(!is.null(depth_of_coverage)){
    gc <- sum(alphabetFrequency(genes)[,2:3]*depth_of_coverage)/sum(alphabetFrequency(genes)*depth_of_coverage)
    return(data.frame(CUBHE = getWeightedCUBHE(codon_table,
                                     highly_expressed,
                                     depth_of_coverage,
                                     genetic_code = genetic_code),
                      GC = gc,
                      GCdiv = abs(0.5-gc),
                      ConsistencyHE = getWeightedConsistency(codon_table,
                                                     highly_expressed,
                                                     depth_of_coverage,
                                                     genetic_code = genetic_code),
                      CUB = getWeightedCUB(codon_table,
                                   highly_expressed,
                                   depth_of_coverage,
                                   genetic_code = genetic_code),
                      CPB = NA,
                      FilteredSequences = filtered$Filtered,
                      nHE = sum(highly_expressed),
                      stringsAsFactors = FALSE))
  } else {
    gc <- sum(alphabetFrequency(genes)[,2:3])/sum(alphabetFrequency(genes))
    return(data.frame(CUBHE = getCUB(codon_table,
                                     highly_expressed,
                                     method = "MILC",
                                     genetic_code = genetic_code),
                      GC = gc,
                      GCdiv = abs(0.5-gc),
                      ConsistencyHE = getCUB(codon_table,
                                             highly_expressed,
                                             method = "consistency",
                                             genetic_code = genetic_code),
                      CUB = getCUB(codon_table,
                                   highly_expressed,
                                   method = "MILCgenomic",
                                   genetic_code = genetic_code),
                      CPB = getCPB(codon_pair_table),
                      FilteredSequences = filtered$Filtered,
                      nHE = sum(highly_expressed),
                      stringsAsFactors = FALSE))
  }
}
