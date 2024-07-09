FROM shengwei/ubuntu2204py312r441:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN apt-get -qq update && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install devtools, biocmanager, matrixstats, biostrings, cordon
RUN conda install --yes conda-forge::r-devtools conda-forge::r-biocmanager  conda-forge::r-matrixstats && \
    conda install --yes bioconda::bioconductor-biostrings bioconda::bioconductor-cordon && \
    conda clean -ya && \
    export TAR=/bin/tar

RUN Rscript -e 'devtools::install_github("jlw-ecoevo/gRodon2")' 

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]
