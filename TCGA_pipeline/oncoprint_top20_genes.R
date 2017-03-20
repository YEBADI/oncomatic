### HISTORY ####################################################################
# Version		Date			    Coder			          Comments
# 3.0 			20/03/2017		Yusef + Emanuela		Functionally produces oncoprint
#                                             for top 20 mutated genes.
### DESCRIPTION ################################################################
#  R pipeline to download TCGAbiolinks data for tumor type
#  and produce an oncoprint for the top 20 most mutated genes.
### PARAMETERS #################################################################
#  I need to pass in paramaters of gene names from the console
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

# Load in the libraries to R.
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

# Download mutations data and generate oncoprint of top 20 genes for all pipes.
pipeline_options <- c("muse", "varscan2", "somaticsniper", "mutect2");
for(pipe in pipeline_options){
	mut.data <- GDCquery_Maf(tumor = tumor.type, pipelines = pipe, 
                           save.csv = TRUE);
	genes.names <- unique(mut.data$Hugo_Symbol);
	all.mut <- matrix(data=0, nrow=length(user.gene.list), ncol=1, 
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
	top.mut.genes <- names( all.mut.ordered[ 1:20 ] );
  # Subset mutation data
	top.mut.data <- mut.data[ all.positions, ];
}
# Oncoprint of top 20 genes.
for(pipe in pipeline_options){
  TCGAvisualize_oncoprint(
      mut = top.mut.data,
      genes = top.mut.genes,
      filename = paste("oncoprint_", tumor.type, "_", pipe, ".pdf", sep=""),
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
}