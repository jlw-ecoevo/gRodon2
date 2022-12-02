**Do you have growth rates matched to genomes for cultured isolates? Email JL <jw4336@terpmail.umd.edu> and she will be happy to incorporate that data into the next version of gRodon in development! They are always on the lookout for more data, and always happy to have more collaborators on board.**

# `gRodon`

`gRodon` is an R package to estimate maximal growth rates of prokaryotes and microbial eukaryotes (**new in v2.0.0**) from genome-wide codon usage statistics. You can find a detailed tutorial (vignette) on how to use the `gRodon` package [here](https://jlw-ecoevo.github.io/gRodon-vignette).

You can use `gRodon` to get maximal growth rate predictions from individual **genomes**, as well-as bulk community-wide average growth rates from **metagenomes**. 

To run `gRodon` you will need a fasta file with your coding sequence (ORFs), as well as a list of highly expressed proteins (typically ribosomal proteins). If you would like to run abundance-weighted metagenome mode you will also need mean depth of coverage estimates for each of your ORFs.

## Installation

The easiest way to install `gRodon` is with [`devtools`](https://github.com/r-lib/devtools).

`devtools::install_github("jlw-ecoevo/gRodon2")`

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
library(gRodon)
library(Biostrings)

# Load in example genome (Streptococcus pyogenes M1, downloaded from RefSeq)
# included with gRodon
path_to_genome <- system.file('extdata',
  'GCF_000349925.2_ASM34992v2_cds_from_genomic.fna.gz',
  package = 'gRodon')
genes <- readDNAStringSet(path_to_genome)

# Search pre-existing annotations for ribosomal proteins, which we
# will use as our set of highly expressed genes
highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case = T)

# Run the gRodon growth prediction pipeline
predictGrowth(genes, highly_expressed)
```

## Using `gRodon` with [`docker`](https://www.docker.com/)

We have compiled a docker image for `gRodon` v1.0.0 (**no eukaryotes, no metagenome_v2 mode**) to ease the installation process. You can pull it to your local computer and run it like this:

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

## `gRodon` with conda

[Susheel Busi has setup a conda environment to run `gRodon`](https://github.com/susheelbhanu/gRodon)  v1.0.0 (**no eukaryotes, no metagenome_v2 mode**) with some helper scripts that made be useful to some users. The `gRodon` developers take no responsibility for the functioning of this code though and all questions should be submitted directly to the author.

## Citation
If you find `gRodon` is useful to your study, please cite us!

#### For prokaryotic prediction: the [`gRodon` paper](https://doi.org/10.1073/pnas.2016810118): 

> JL Weissman, Shengwei Hou, Jed A. Fuhrman. Estimating maximal microbial growth rates from cultures, metagenomes, and single cells via codon usage patterns. Proceedings of the National Academy of Sciences Mar 2021, 118 (12) e2016810118; DOI: 10.1073/pnas.2016810118

```
@article {Weissmane2016810118,
	author = {Weissman, JL and Hou, Shengwei and Fuhrman, Jed A.},
	title = {Estimating maximal microbial growth rates from cultures, metagenomes, and single cells via codon usage patterns},
	volume = {118},
	number = {12},
	elocation-id = {e2016810118},
	year = {2021},
	doi = {10.1073/pnas.2016810118},
	publisher = {National Academy of Sciences},
	issn = {0027-8424},
	URL = {https://www.pnas.org/content/118/12/e2016810118},
	eprint = {https://www.pnas.org/content/118/12/e2016810118.full.pdf},
	journal = {Proceedings of the National Academy of Sciences}
}
```

#### For eukaryotic prediction: the [`gRodon2` paper](https://doi.org/10.1101/2021.10.15.464604): 

> JL Weissman, Edward-Robert O Dimbo, Arianna I Krinos, Christopher Neely, Yuniba Yagues, Delaney Nolin, Shengwei Hou, Sarah Laperriere, David A Caron, Benjamin L Tully, Harriet Alexander, Jed A Fuhrman. Estimating the maximal growth rates of eukaryotic microbes from cultures and metagenomes via codon usage patterns. bioRxiv 2021.10.15.464604; DOI: https://doi.org/10.1101/2021.10.15.464604

```
@article {Weissman2021.10.15.464604,
	author = {Weissman, JL and Dimbo, Edward-Robert O and Krinos, Arianna I and Neely, Christopher and Yagues, Yuniba and Nolin, Delaney and Hou, Shengwei and Laperriere, Sarah and Caron, David A and Tully, Benjamin L and Alexander, Harriet and Fuhrman, Jed A},
	title = {Estimating the maximal growth rates of eukaryotic microbes from cultures and metagenomes via codon usage patterns},
	elocation-id = {2021.10.15.464604},
	year = {2021},
	doi = {10.1101/2021.10.15.464604},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2021/10/16/2021.10.15.464604},
	eprint = {https://www.biorxiv.org/content/early/2021/10/16/2021.10.15.464604.full.pdf},
	journal = {bioRxiv}
}
```

#### For metagenomic prediction: the [`Metagenome Mode v2` paper](https://doi.org/10.1101/2022.04.12.488109): 

> JL Weissman, Marie Peras, Tyler P Barnum, Jed A Fuhrman. Benchmarking community-wide estimates of growth potential from metagenomes using codon usage statistics. bioRxiv 2022.04.12.488109; DOI: https://doi.org/10.1101/2022.04.12.488109

```
@article {Weissman2022.04.12.488109,
	author = {Weissman, JL and Peras, Marie and Barnum, Tyler P. and Fuhrman, Jed A.},
	title = {Benchmarking community-wide estimates of growth potential from metagenomes using codon usage statistics},
	elocation-id = {2022.04.12.488109},
	year = {2022},
	doi = {10.1101/2022.04.12.488109},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2022/04/13/2022.04.12.488109},
	eprint = {https://www.biorxiv.org/content/early/2022/04/13/2022.04.12.488109.full.pdf},
	journal = {bioRxiv}
}
```

#### We also encourage you to cite `gRodon`'s dependencies:

> Elek A, Kuzman M, Vlahovicek K (2020). coRdon: Codon Usage Analysis and Prediction of Gene Expressivity. R package version 1.8.0, https://github.com/BioinfoHR/coRdon

> Pagès H, Aboyoun P, Gentleman R, DebRoy S (2020). Biostrings: Efficient manipulation of biological strings. R package version 2.58.0, https://bioconductor.org/packages/Biostrings.

> Henrik Bengtsson (2021). matrixStats: Functions that Apply to Rows and Columns of  Matrices (and to Vectors). R package version 0.58.0. https://CRAN.R-project.org/package=matrixStats

As well as the original paper describing the MILC statistic:

> Supek, Fran, and Kristian Vlahovicek. “Comparison of codon usage measures and their applicability in prediction of microbial gene expressivity.” BMC bioinformatics vol. 6 182. 19 Jul. 2005, doi:10.1186/1471-2105-6-182

