# `gRodon`

`gRodon` is an R package to estimate maximal growth rates of prokaryotic organisms from genome-wide codon usage statistics. You can find a detailed tutorial (vignette) on how to use the `gRodon` package [here](https://jlw-ecoevo.github.io/gRodon-vignette).

You can use `gRodon` to get maximal growth rate predictions from **genomes**, as well-as bulk community-wide average growth rates from **metagenomes**. 

To run `gRodon` you will need a fasta file with your coding sequence (ORFs), as well as a list of highly expressed proteins (typically ribosomal proteins). If you would like to run abundance-weighted metagenome mode you will also need mean depth of coverage estimates for each of your ORFs.

## Installation

The easiest way to install `gRodon` is with [`devtools`](https://github.com/r-lib/devtools).

`devtools::install_github("jlw-ecoevo/gRodon")`

gRodon has a few dependencies - namely the Biostrings, coRdon, and matrixStats packages which are bioconductor packages and cannot be installed via CRAN. To install them run the following:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")
BiocManager::install("coRdon")
install.packages("matrixStats")
```

## A simple example

Currently `gRodon` only has a single function available to users: `predictGrowth`. 

To see the details of how this function works type `?predictGrowth()`.

A minimal example with data included in the package is:

```
# Load in example genome (Streptococcus pyogenes M1, downloaded from RefSeq)
# included with gRodon
path_to_genome <- system.file('extdata',
  'GCF_000349925.2_ASM34992v2_cds_from_genomic.fna',
  package = 'gRodon')
genes <- readDNAStringSet(path_to_genome)

# Search pre-existing annotations for ribosomal proteins, which we
# will use as our set of highly expressed genes
highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case = T)

# Run the gRodon growth prediction pipeline
predictGrowth(genes, highly_expressed)
```

## Using `gRodon` with [`docker`](https://www.docker.com/)

We have compiled a docker image for `gRodon` to ease the installation process. You can pull it to your local computer and run it like this:

```bash
# pull the image
$ docker pull shengwei/grodon:latest

# start an interactive container
$ docker run -ti --rm shengwei/grodon:latest
```

Now you're inside of docker container, let's start an `R` session
```bash
$ root@5218b31cd695:/mnt# R
```

Now you're inside of `R` REPL of the docker container, let's test `gRodon`: 
```
> library(gRodon)
> library(Biostrings)
> path_to_genome <- system.file('extdata',
  'GCF_000349925.2_ASM34992v2_cds_from_genomic.fna.gz',
  package = 'gRodon')
> genes <- readDNAStringSet(path_to_genome)
> highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case = T)
> predictGrowth(genes, highly_expressed)
```

To mount your own data volume and run in non-interactive mode, please refer to [this](https://hub.docker.com/r/shengwei/das_tool/) example.

If you want to modify and build your own docker image, the source code can be found [here](https://github.com/housw/Bioinfo-Dockerfiles/blob/master/gRodon).

## Citation
If you find `gRodon` is useful to your study, please cite: 

> Jake L. Weissman, Shengwei Hou, Jed A. Fuhrman. Estimating maximal microbial growth rates from cultures, metagenomes, and single cells via codon usage patterns. bioRxiv 2020.07.25.221176; doi: https://doi.org/10.1101/2020.07.25.221176




