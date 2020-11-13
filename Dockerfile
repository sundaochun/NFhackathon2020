## Start from this Docker image
FROM rocker/tidyverse
## use rocker as a base image

## install synapser reqs
RUN apt-get update -y
RUN apt-get install -y dpkg-dev zlib1g-dev libssl-dev libffi-dev
RUN apt-get install -y curl libcurl4-openssl-dev

## install synapser
RUN R -e "install.packages('synapser', repos=c('http://ran.synapse.org', 'http://cran.fhcrc.org'))"
RUN R -e "install.packages('synapserutils', repos=c('http://ran.synapse.org', 'http://cran.fhcrc.org'))"


## install bioconductor packages
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install(c('GSVA', 'GSEABase', 'org.Hs.eg.db', 'limma', 'GOsummaries', 'GSVAdata', 'biomaRt', 'maftools', 'Biostrings','dendextend','WGCNA','ggdendro','GOSemSim','Biobase'))"

## install cran packages
RUN R -e "install.packages(c('gProfileR', 'umap', 'dbscan', 'ggfortify', 'pheatmap', 'ggpubr', 'DT', 'here', 'reshape2', 'RColorBrewer','ggplot2','gplots'))"

RUN mkdir /home/rstudio/output
COPY 0_setup.Rmd /home/rstudio/0_setup.Rmd
COPY 1_FinalCode.Hackathon2020.R /home/rstudio/1_FinalCode.Hackathon2020.R