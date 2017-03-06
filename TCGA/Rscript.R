# setting up TCGA R environment if not yet present

source("https://bioconductor.org/biocLite.R")

biocLite("TCGAbiolinks")
biocLite("dplyr")
biocLite("DT")

# setting up TCGA R libraries

library(TCGAbiolinks)
library(dplyr)
library(DT)

