
shuffleGenes <- function(gene){
  DNAString(stringi::stri_rand_shuffle(as.character(gene))) %>%
    return()
}

singleMILC <- function(gene,n=100){
  bg <- replicate(n,shuffleGenes(gene)) %>% DNAStringSet()
  # print(bg)
  # print(gene)
  genes_bg <- c(DNAStringSet(gene),bg)
  l1 <- logical(n+1)
  l1[1] <- T
  MILC(codonTable(genes_bg))[1,1] %>%
    return()
}


CUBi <- function(genes,highly_expressed,n_le=100){
  genes_list <- as.list(as.character(genes)) %>%
    lapply(DNAString)
  genes_list_HE <- genes_list[highly_expressed]
  genes_list_LE <- genes_list[!highly_expressed][sample(1:sum(!highly_expressed),n_le,replace=T)]

  x_HE <- lapply(genes_list_HE,singleMILC) %>% unlist()
  x_LE <- lapply(genes_list_LE,singleMILC) %>% unlist()

  return(c(median(x_LE),median(x_HE)))
}

getCUB <- function(fna_tab, highly_expressed, method = "MILC", genetic_code = "11"){

  if(method == "MILC"){
    # MILC estimate of codon usage bias from coRdon
    #   - Bias of highly expessed genes, with genome-wide codon usage considered
    #     expected codon usage
    x <- MILC(fna_tab,id_or_name2 = genetic_code)
    # print(x[highly_expressed, 1])
    x[highly_expressed, 1] %>% median() %>% return()

  } else if(method == "MILCgenomic"){
    # MILC estimate of codon usage bias from coRdon
    #   - Bias of all genes, with genome-wide codon usage considered
    #     expected codon usage
    x <- MILC(fna_tab,id_or_name2 = genetic_code)
    x[, 1] %>% median() %>% return()

  # } else if(method == "MILC_i"){
  #   # MILC estimate of codon usage bias from coRdon
  #   #   - Bias of all genes, with genome-wide codon usage considered
  #   #     expected codon usage
  #   y <- numeric()
  #   for(i in 1:sum(highly_expressed)){
  #     print(MILC(fna_tab[highly_expressed][i]))
  #     y <- c(y,MILC(fna_tab[highly_expressed][i]))
  #   }
  #   print(y)
  #   y %>% unlist() %>% mean() %>% return()
  #
  # } else if(method == "MILCgenomic_i"){
  #   # MILC estimate of codon usage bias from coRdon
  #   #   - Bias of all genes, with genome-wide codon usage considered
  #   #     expected codon usage
  #   y <- numeric()
  #   for(i in 1:length(fna_tab)){
  #     y <- c(y,MILC(fna_tab[i]))
  #   }
  #   y %>% unlist() %>% mean() %>% return()

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


getWeightedCUBi <- function(genes,highly_expressed,depth_of_coverage,n_le=100){
  genes_list <- as.list(as.character(genes)) %>%
    lapply(DNAString)
  samp_le <- sample(1:sum(!highly_expressed),n_le,replace=T)
  genes_list_HE <- genes_list[highly_expressed]
  genes_list_LE <- genes_list[!highly_expressed][samp_le]

  x_HE <- lapply(genes_list_HE,singleMILC) %>% unlist()
  x_LE <- lapply(genes_list_LE,singleMILC) %>% unlist()

  wHE <- x_HE %>%
    weightedMedian(., w = depth_of_coverage[highly_expressed])
  wLE <- x_LE %>%
    weightedMedian(., w = depth_of_coverage[!highly_expressed][samp_le])

  return(c(median(wLE),median(wHE)))
}


# getWeightedCUBHE_i <- function(fna_tab, highly_expressed, depth_of_coverage, genetic_code = "11"){
#   y <- numeric()
#   for(i in 1:sum(highly_expressed)){
#     y <- c(y,MILC(fna_tab[highly_expressed][i]))
#   }
#   sum(unlist(y)*depth_of_coverage[highly_expressed]/sum(depth_of_coverage[highly_expressed])) %>%
#     return()
# }
#
# getWeightedCUB_i <- function(fna_tab, highly_expressed, depth_of_coverage, genetic_code = "11"){
#   y <- numeric()
#   for(i in 1:length(fna_tab)){
#     y <- c(y,MILC(fna_tab[i]))
#   }
#   sum(unlist(y)*depth_of_coverage/sum(depth_of_coverage)) %>%
#     return()
# }

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



getCodonStatistics_i <- function(genes, highly_expressed, fragments = FALSE, depth_of_coverage = NULL, genetic_code = "11",n_le=100){

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
    cubi <- getWeightedCUBi(genes,highly_expressed,depth_of_coverage,n_le=n_le)
    return(data.frame(CUBHE = cubi[2],
                      GC = gc,
                      GCdiv = abs(0.5-gc),
                      ConsistencyHE = getWeightedConsistency(codon_table,
                                                             highly_expressed,
                                                             depth_of_coverage,
                                                             genetic_code = genetic_code),
                      CUB = cubi[1],
                      CPB = NA,
                      FilteredSequences = filtered$Filtered,
                      nHE = sum(highly_expressed),
                      stringsAsFactors = FALSE))
  } else {
    cubi <- CUBi(genes,highly_expressed,n_le=n_le)
    gc <- sum(alphabetFrequency(genes)[,2:3])/sum(alphabetFrequency(genes))
    return(data.frame(CUBHE = cubi[2],
                      GC = gc,
                      GCdiv = abs(0.5-gc),
                      ConsistencyHE = getCUB(codon_table,
                                             highly_expressed,
                                             method = "consistency",
                                             genetic_code = genetic_code),
                      CUB = cubi[1],
                      CPB = NA,
                      FilteredSequences = filtered$Filtered,
                      nHE = sum(highly_expressed),
                      stringsAsFactors = FALSE))
  }
}
