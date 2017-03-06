# setting up TCGA R environment if not yet present

source("https://bioconductor.org/biocLite.R")

biocLite("TCGAbiolinks")
biocLite("dplyr")
biocLite("DT")

# setting up TCGA R libraries

library(TCGAbiolinks)
library(dplyr)
library(DT)


# downloading XML-indexed data
clinical <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")

datatable(clinical, filter = 'top', options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), rownames = FALSE)
