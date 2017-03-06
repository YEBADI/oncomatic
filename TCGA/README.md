#The Cancer Genome Atlas

This is the TCGA readme.

## How to use bioconductor TCGA biolinks to download and analyse data
http://bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/tcgaBiolinks.html


## Paper summaries

### TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4856967/

TCGA is a large collection of clinical and molecular phenotypes of more than 10 000 tumor patients across 33 different tumor types. 

From this, TCGA has published over 20 marker papers detailing the genomic and epigenomic alterations associated with these tumor types. 

Implementing novel methods can elucidate new biological pathways and diagnostic markers. 

Mining the TCGA data presents several bioinformatics challenges, like data retrieval and integration with clinical data and other molecular data types (e.g. RNA and DNA methylation).

An R/Bioconductor package was developed (TCGAbiolinks) to address these challenges and offer bioinformatics solutions by using a guided workflow to allow users to query, download and perform integrative analyses of TCGA data.

Four different TCGA tumor types (Kidney, Brain, Breast and Colon) were ued as case studies to illustrate examples of reproducibility, integrative analysis and utilization of different Bioconductor packages to advance and accelerate novel discoveries.

# initial investigation

Obtain mutations data from TCGA (breast or pancreatic)

play around with R analysis

need to subset data to obtain top hit dif gene expression + genes of interest to end-user
