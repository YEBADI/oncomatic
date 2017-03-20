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

# Set working directory to home.
setwd(home)

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
ProjectID <- "TCGA-BRCA"; # The dataset name as defined by TCGAbiolinks.
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
# Biolinks command to download the general clinical data for all BRCA samples
clin.data <- GDCquery_clinic( "TCGA-BRCA", "clinical" );
# Check dimensions
dim(clin.data)
# Define tumor type according to TCGA format e.g. BRCA (breast), PAAD (Pancreas)
tumor.type <-  "BRCA"

# Write clinical data matrix to target file
write.table( clin.data, file = paste(tumor.type,"_", ProjectName, 
                          "_clinical.txt", sep=""), sep="\t", row.names=FALSE );

# Subsetting the clinical data to covariates of interest to new data matrix
clin.covariates.for.oncodata <- c("race", "gender", "vital_status", 
                                          "tumor_stage", "bcr_patient_barcode");
clin.forvisual <- clin.data[ , clin.covariates.for.oncodata]
clin.data.slimmed <- clin.data[ , clin.covariates ];

# Check dimensions of slimmed data matrix
dim(clin.data.slimmed) #1097 9

#===============================================================================
#    Access and download the mutation data for tumor type and
#    generate an oncoprint for each pipeline
#    Note: As mentioned before, should be softcoded to let user decide cancer
#    type as well as pick pipeline option (passed in as parameter).
#===============================================================================

# [1] MuTect2
# * Applies Bayesian classifier to detect somatic mutations with very low allele 
# fractions and requires only a few supporting reads however is confounded by 
# low tumor purity.
# * Has high sensitivity and calls mutations with allelic fractions as low as 
# 0.1 and below.
# * Low specifiticity.
# * Applies carefully tuned hard filters to compensate for specificity issues.

# [2] Varscan2
# * Is generally outperformed by MuSE and MuTect2.
# * Low sensitivity and fails to pick up somatic SNVs of low allelic fraction
# as supresses mutations below allelic threshold.
# * Sensitivity can be improved but this drastically drops specificitiy and 
# returns high levels of false positives.
# * Outperformed MuTect to identify variants present at 10%, 5%, 2.5% and 1%
# at sequencing depths of 100x, 250x, 500x and 1000x respectively.
# * See (Stead et al, 2013).
# * However, Varscan2 circumvents confounding factor of tumor purity and extreme
# read depth as does not use probabilistic framework to detect variants and
# assess confidence in them but uses a robust heuristic/statistic approach to 
# call variants that meet desired thresholds for read depth, base quality, 
# variant allele frequency, and statistical significance.

# [3] MuSE
# * Has outperformed MuTect2 on calling variants from ACC TCGA data.
# * Maximum likelihood or the Markov chain Monte Carlo (MCMC) method estimates 
# Model parameters. 
# * Variants are classified into somatic, germ-line, and reversal to the 
# homozygous reference by comparing the somatic variant allele fraction (π) 
# between the paired tumor–normal samples.
# * A sample-specific error model to account for tumor heterogeneity and 
# identify cutoffs is built.
# * Filters reduce the number of false positives by considering the sequence 
# context surrounding the point mutations. 

# [4] SomaticSniper
# Returns high level of false positives and many of these are not in agreement 
# with other InDel callers.
# Most of the literature shows it is outperformed by most other variant callers.

# Verdict
#* TCGA-BRCA analysis is suitable with either Mutect2 or MuSE. Ideally both 
# should be run and the results correlated / contrasted.
#* TCGA-PAAD analysis would only be advisable with Varscan2 due to significantly
# low tumor purity (circa 40%).

# Download mutations data & visualise clinical & mutations data in an oncoprint

# Oncoprint of user requested genes
user.gene.request = readline('Please type the genes of interest in capital 
       letters, separated by " ; " with NO SPACE (e.g. " ATM;BRCA1;CHEK2 "): ')
# COMMONLY MUTATED BREAST CANCER GENES
# ATM;p53;BRCA1;BRCA2;PTEN;CHEK2;PALB2;STK11;BARD1;BRIP1;CASP8;CDH1;CHEK2
user.gene.list <- unlist(strsplit(user.gene.request,";"))
pipeline_options <- c("muse", "varscan2", "somaticsniper", "mutect2");
for(pipe in pipeline_options){
  mut.data <- GDCquery_Maf(tumor = tumor.type, pipelines = "mutect2", 
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
  for(pipe in pipeline_options){
    TCGAvisualize_oncoprint(
        mut = user.choice.mut.data,
        genes = user.choice.gene.names,
        filename = paste("oncoprint_", tumor.type, user.gene.request, 
                                        "_", pipe, ".pdf", sep=""),
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
}