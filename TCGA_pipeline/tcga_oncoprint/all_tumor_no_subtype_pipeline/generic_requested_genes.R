### HISTORY ####################################################################
# Version   Date          Coder               Comments
# 1.0       20/03/2017    Yusef               Functionally produces oncoprint
#                                             for user picked mutated genes.
### DESCRIPTION ################################################################
#  R pipeline to download TCGAbiolinks data for tumor type
#  and produce an oncoprint for user picked mutated genes.
### PARAMETERS #################################################################
# 
### FUNCTIONS ##################################################################

################################ USER INPUT ####################################
# To place in README.txt.
# Please give two params:
# 1. Tumor type, in capital letters; either "BRCA" or "PAAD".
#
# 2. The genes of interest in capital letters, separated by " - " with NO SPACE 
# (e.g. " ATM-BRCA1-CHEK2 "): ')
# 3. The pipeline of choice; "mutect2", "varscan2", "muse", or "somaticsniper".
########################## PARAMETER INPUT #####################################
cat("-- reading arguments\n", sep = "");
cmd_args = commandArgs(trailingOnly=TRUE);
for (arg in cmd_args) cat("  ", arg, "\n", sep="");

args1 <- cmd_args[1] # "BRCA" or "PAAD"

args2 <- cmd_args[2] #select genes
# "ATM-P53-BRCA1-BRCA2-PTEN-CHEK2-PALB2-STK11-BARD1-BRIP1-CASP8-CDH1-CHEK2"

# Pipeline_options are "muse", "varscan2", "somaticsniper", or "mutect2".
args3 <- cmd_args[3] # "mutect2"

if(args1 == "BRCA" ){
  ProjectID <- "TCGA-BRCA"
  # Define tumor type according to TCGA format e.g. BRCA (breast), PAAD (Pancreas)
  tumor.type <-  "BRCA"
} else{
  ProjectID <- "TCGA-PAAD"
  tumor.type <-  "PAAD"
}

user.gene.request <- args2
user.gene.list <- unlist(strsplit(user.gene.request,"-"))

################################### MAIN #######################################
#===============================================================================
#    Preparation of environment
#===============================================================================
# Identify and install missing libraries.
TCGA.libs <- c("TCGAbiolinks","SummarizedExperiment");
new.libs <- TCGA.libs[!(TCGA.libs %in% installed.packages() [,"Package"])]
if(length(new.libs)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite(new.libs)
}

# loads in the libraries to R
library("TCGAbiolinks")
library("SummarizedExperiment")

#############################  MAIN  ###########################################
#===============================================================================
#    Preparation of environment
#===============================================================================
# Define variables.
ProjectDir <- "mutations";
ProjectName <-  "mutations";

# create mutations directory and go into it
if (file.exists(ProjectDir)){
  cat( paste("\nFolder with project [", ProjectName, "] already exists. 
         The requested data will be stored in [", ProjectDir, "]\n", sep=" ") );
} else {
  dir.create(ProjectDir);
}
setwd(ProjectDir);

#===============================================================================
#    Data query and download clincal data
#    Note: tumor can be softcoded as param input from command line.
#    Note: can develop program and let user pick clin.data terms.
#===============================================================================
# Biolinks command to download the general clinical data for all patient samples
clin.data <- GDCquery_clinic( ProjectID, "clinical" );
# Check dimensions
dim(clin.data)

# Subsetting the clinical data to covariates of interest to new data matrix
clin.covariates.for.oncodata <- c("race", "gender", "vital_status", 
                                          "tumor_stage", "bcr_patient_barcode");
clin.forvisual <- clin.data[ , clin.covariates.for.oncodata]

# Check dimensions of slimmed data matrix
dim(clin.forvisual) #1097 5

#===============================================================================
#    Access and download the mutation data for tumor type and
#    generate an oncoprint for each pipeline
#    Note: As mentioned before, should be softcoded to let user decide cancer
#    type as well as pick pipeline option (passed in as parameter).
#===============================================================================
# Download mutations data & visualise clinical & mutations data in an oncoprint
# Oncoprint of user requested genes
mut.data <- GDCquery_Maf(tumor = tumor.type, pipelines = args3, 
                         save.csv = TRUE);
mut.data.matrix <- matrix(data=0, nrow=length(user.gene.list), ncol=1, 
        dimnames=list(rownames=user.gene.list, colnames="number_of_reports"));
all.gene.positions <- NULL
for(gene in user.gene.list) {
  user.gene.list.items <- which(mut.data$Hugo_Symbol %in% user.gene.list)
  mut.data.matrix[gene, "number_of_reports"] <- length(user.gene.list.items)
  all.gene.positions <- c(all.gene.positions, user.gene.list.items)
}
mut.data.matrix.ordered <- mut.data.matrix[ order(mut.data.matrix[, 
                                    "number_of_reports"], decreasing=TRUE),];
user.choice.gene.names <- names( mut.data.matrix.ordered[ 1:20 ] );
user.choice.mut.data <- mut.data[ all.gene.positions, ];
TCGAvisualize_oncoprint(
    mut = user.choice.mut.data,
    genes = user.choice.gene.names,
    filename = paste("oncoprint_", tumor.type, "-", user.gene.request, 
                                    "_", args3, ".pdf", sep=""),
    annotation = clin.forvisual,
    color=c("background"="#CCCCCC","DEL"="purple",
                  "INS"="yellow","SNP"="brown"),
    rows.font.size= 8,
    width = 5,
    heatmap.legend.side = "right",
    dist.col = 0,
    label.font.size = 6
    );
}
  dev.off();