# setting up TCGA R environment if not yet present

source("https://bioconductor.org/biocLite.R")

biocLite("TCGAbiolinks")
biocLite("dplyr")
biocLite("DT")

# setting up TCGA R libraries

library(TCGAbiolinks)
library(dplyr)
library(DT)

###################################
#### ignore this section
#### datatable needs to run on apocrita
#### can't run on PC
#
##### downloading XML-indexed data & setting up datatable#
# clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
# # ^this works
# datatable(clinical, filter = 'top', options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), rownames = FALSE)
# # ^ this says can't run on laptop as data too big 
##################################


# We can easily search GDC data for data to download
# using the GDCquery function.
# Using filters akin to the filters in the TCGA portal

# we won't need to look at methylations
# also primary concern is breast and pancreatic tissues

# direct query
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")


query.mirna <- GDCquery(project = "TCGA-BRCA", 
                        data.category = "Transcriptome Profiling", 
                        data.type = "miRNA Expression Quantification",
                        sample.type = c("Primary solid Tumor","Solid Tissue Normal"))

query.mirna <- GDCquery(project = "TCGA-ACC", 
                        data.category = "Transcriptome Profiling", 
                        data.type = "miRNA Expression Quantification",
                        sample.type = c("Primary solid Tumor","Solid Tissue Normal"))