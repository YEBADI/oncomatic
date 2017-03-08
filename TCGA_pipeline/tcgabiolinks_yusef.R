### HISTORY ##############################################################################
# Version		Date			Coder			Comments
# 1.0 			07/03/2017		Emanuela		init
#
### DESCRIPTION ##########################################################################
# script for Yusef - become familiar with TCGAbiolinks
#
### PARAMETERS ###########################################################################
# 
### FUNCTIONS ############################################################################
# install libraries if not present
TCGA.libs <- c("TCGAbiolinks","SummarizedExperiment");

if (length(setdiff(TCGA.libs, rownames(installed.packages()))) > 0) {
   source("http://bioconductor.org/biocLite.R")
   biocLite(setdiff(TCGA.libs, rownames(installed.packages())))
}

# load libraries
library("TCGAbiolinks")
library("SummarizedExperiment")


### MAIN #################################################################################
# define variables
ProjectDir <- "mutations"; # in this instance, I am already in Desktop and want to work from /Desktop/mutations/
ProjectID <- "TCGA-BRCA"; # as defined by TCGAbiolinks
ProjectName <-  "mutations";

# set/create the project directory i.e. "/Desktop/mutations/"
if (file.exists(ProjectDir)){

	cat( paste("\nFolder with project [", ProjectName, "] already exists. The requested data will be stored in [", ProjectDir, "]\n", sep=" ") );
} else {
	dir.create(ProjectDir);

}

# change working directory to the project workspace
setwd(ProjectDir);


#===============================================================================
#    Data query and download
#===============================================================================
# ACCESS AND DOWNLOAD THE CLINICAL DATA
clin.data <- GDCquery_clinic( "TCGA-BRCA", "clinical" );
dim(clin.data) #1097 39

# Define type of tumour e.g. BRCA (breast)
tumour.type <-  "BRCA"

# Write clinical data to outfile - this is your target file
#write.table( clin.data, file = paste(tumour.type,"_", ProjectName, "_clinical.txt", sep=""), sep="\t", row.names=FALSE );

# what clinical information do you want to appear in your figure? view you options by:
colnames(clin.data);

# Subset clinical data to covariates of interest
clin.covariates <- c("bcr_patient_barcode","disease","gender","race","vital_status");
clin.data.slimmed <- clin.data[ , clin.covariates ];
dim(clin.data.slimmed) #1097 5


# ACCESS AND DOWNLOAD THE MUTATION DATA
## [Yusef: briefly read a little on each pipeline i.e. an overview of mutect2 / varscan etc]
### can you please write a couple of sentences on each? I don't mind if you do this as an overview of all the methods or a couple of lines about each
### in your opinion, which method would be most applicable to breast/pancreatic cancer? i.e. look at literature
### should we give users the option to choose what method they want to apply?

pipelines <- c("muse", "varscan2", "somaticsniper", "mutect2");

# for each methods, create a file of all the mutations
for(pipe in pipelines){

	# for each method (pipe) obtain the mutation data
	mut.data <- GDCquery_Maf( tumour = tumour.type, save.csv = TRUE, pipelines = pipe );

	# identify all genes (symbol) reported as mutated
	unique.genes <- unique(mut.data$Hugo_Symbol);

	# create input data matrix to collect information about the most reported genes
	all.mut <- matrix(data=0, nrow=length(unique.genes), ncol=1, dimnames=list(rownames=unique.genes, colnames="no_reports"));

	# identify position of each gene and populate the matrix
	all.positions <- NULL;

	for(gene in unique.genes) {
		gene.position <- which(mut.data$Hugo_Symbol == gene);
		all.mut[gene, "no_reports"] <- length(gene.position);
		all.positions <- c(all.positions, gene.position );
	}

	# order by most reported gene
	all.mut.ordered <- all.mut[ order(all.mut[, "no_reports"], decreasing=TRUE),];

	# define the "top" 20 genes to be used in the visualisation process and subset mutation data
	top.mut.genes <- names( all.mut.ordered[ 1:20 ] );
	top.mut.data <- mut.data[ all.positions, ];

	# make sure barcodes in the clinical data correspond to those in the mutation data
	top.mut.bar <- NULL;

	# if you look at what I am printing 
	for( i in 1:length(top.mut.data$Tumor_Sample_Barcode) ){
		tmp.bar.split <- unlist(strsplit(top.mut.data$Tumor_Sample_Barcode[i], split="-", fixed=TRUE));
		tmp.bar.join <- paste(tmp.bar.split[1], "-", tmp.bar.split[2], "-", tmp.bar.split[3], sep = "")
		cat("\n", tmp.bar.split, "\t", tmp.bar.join);
		top.mut.bar <- rbind(top.mut.bar, tmp.bar.join)
	}

	rownames(top.mut.bar)<- top.mut.bar[,1];
	unique.bar <- unique(rownames(top.mut.bar));
	top.clin.data <- clin.data.slimmed[ , unique.bar]

	# visualise mutation data using oncoprint using ComplexHeatmap package
	TCGAvisualize_oncoprint(
		mut = top.mut.data,
		genes = top.mut.data$Hugo_Symbol, 
		filename = paste("oncoprint_", tumour.type, "_", pipe,".pdf", sep=""),
		annotation = top.clin.data,
		color = c("background"="#CCCCCC","DEL"="purple","INS"="yellow","SNP"="brown"),
		rows.font.size = 8,
		width = 5,
		heatmap.legend.side = "right",
		dist.col = 0,
		label.font.size = 6
	);
	dev.off();

}

## [Yusef: now imagine that I am a researcher interested in specific genes]
## Look up commonly mutated genes in breast cancer and amend the code so that it:
### reads in a user-supplied gene list
### generates the plot for these genes of interest
