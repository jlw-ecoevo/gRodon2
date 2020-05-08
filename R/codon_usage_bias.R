

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

getCodonStatistics <- function(genes, highly_expressed, metagenome){

  if(sum(highly_expressed) == 0){
    stop("No highly expressed genes?")
  }

  #Remove short sequences and sequences that are not multiples of 3
  filtered <- filterSeq(genes, highly_expressed)
  genes <- filtered$Genes
  highly_expressed <- filtered$HE

  # codon table
  codon_table <- codonTable(genes)

  # codon pair counts
  codon_pair_table <- getPairCounts(genes)

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
