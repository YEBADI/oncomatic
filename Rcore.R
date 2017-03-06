##### Analysis of datasets 2017 Yusef Badi ###########

rm (list =ls())  # clear the R workspace

######### loading the required packages ###############


source("https://bioconductor.org/biocLite.R")

biocLite("ggplot2")
biocLite("limma")
biocLite("edgeR")
biocLite("dendextend")
biocLite("calibrate")
biocLite("survival")

###########   setting up libraries  #########

library(ggplot2)

# we will use ggplot2 to generate graphs, then use 
# plotly javascript to make interactive graphs

library(limma)
library(edgeR)
library(dendextend)
library(calibrate)
library(survival)

#### read in gene names from table data input ###
