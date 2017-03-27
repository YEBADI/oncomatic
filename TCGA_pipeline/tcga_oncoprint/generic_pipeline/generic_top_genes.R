### HISTORY ####################################################################
# Version		Date			    Coder			          Comments
# 5.0 			23/03/2017		Yusef + Emanuela		Functionally produces oncoprint
#                                             for top number of mutated genes.
### DESCRIPTION ################################################################
#  R pipeline to download TCGAbiolinks data for tumor type
#  and produce an oncoprint for the top 20 most mutated genes.
### PARAMETERS #################################################################
# User gives param either "BRCA" or "PAAD" to pick tumor type.
# User gives param to select number of top genes to show.
# User gives param to choose pipeline.
### FUNCTIONS ##################################################################

################################ USER INPUT ####################################
# To place in README.txt.
# Please give one param:
# 1. Tumor type, in capital letters; either "BRCA" or "PAAD".
#
########################## PARAMETER INPUT #####################################
cat("-- reading arguments\n", sep = "");
cmd_args = commandArgs(trailingOnly=TRUE);
for (arg in cmd_args) cat("  ", arg, "\n", sep="");

args1 <- cmd_args[1] # "BRCA" or "PAAD"

args2 <- cmd_args[2] # 20 # top number of genes to display (e.g. 20)

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

#############################  MAIN  ###########################################
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

# Load in the libraries to R.
library("TCGAbiolinks")
library("SummarizedExperiment")

# Define variables.
ProjectDir <- "mutations";
# The dataset name as defined by TCGAbiolinks.
ProjectName <-  "mutations";

# Create and enter "mutations" directory.
if (file.exists(ProjectDir)){
	cat( paste("\nFolder with project [", ProjectName, "] already exists. The 
             requested data will be stored in [", ProjectDir, "]\n", sep=" ") );
} else {
	dir.create(ProjectDir);
}
setwd(ProjectDir);

#===============================================================================
#    Data query and download clincal data
#    Note: tumor can be softcoded as param input from command line.
#    Note: can develop program and let user pick clin.data terms.
#===============================================================================
# Biolinks command to download the general clinical data for all tumor samples
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
#    generate an oncoprint using chosen pipeline.
#===============================================================================

# Download mutations data and generate oncoprint of top 20 genes for all pipes.
mut.data <- GDCquery_Maf(tumor = tumor.type, pipelines = args3, 
                         save.csv = TRUE);
genes.names <- unique(mut.data$Hugo_Symbol);
all.mut <- matrix(data=0, nrow=length(genes.names), ncol=1, 
          dimnames=list(rownames=genes.names, colnames="number_of_reports"));
# Identify position of each gene and populate the matrix
all.positions <- NULL;
for(gene in genes.names) {
	gene.position <- which(mut.data$Hugo_Symbol == gene);
	all.mut[gene, "number_of_reports"] <- length(gene.position);
	all.positions <- c(all.positions, gene.position );
}
# Order the whole all.mut matrix by the most reported gene highest to lowest.
all.mut.ordered <- all.mut[ order(all.mut[, "number_of_reports"], 
                                  decreasing=TRUE),];
# Isolate top 20 gene names.
top.mut.genes <- names( all.mut.ordered[ 1:args2 ] );
# Subset mutation data
top.mut.data <- mut.data[ all.positions, ];
# Oncoprint of top 20 genes.
TCGAvisualize_oncoprint(
    mut = top.mut.data,
    genes = top.mut.genes,
    filename = paste("top", args2,"_oncoprint_", tumor.type, "_", 
                     args3, ".pdf", sep=""),
    annotation = clin.forvisual,
    color=c("background"="#CCCCCC","DEL"="purple",
            "INS"="yellow","SNP"="brown"),
    rows.font.size= 8,
    width = 5,
    heatmap.legend.side = "right",
    dist.col = 0,
    label.font.size = 6
  );

dev.off();