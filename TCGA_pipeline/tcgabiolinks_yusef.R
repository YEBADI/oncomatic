### HISTORY ##############################################################################
# Version		Date			Coder			Comments
# 2.0 			13/03/2017		Yusef + Emanuela		init
#
### DESCRIPTION ##########################################################################
# Edited script by Yusef (source v1.0 by Emanuela) to become familiar with TCGAbiolinks
#
### PARAMETERS ###########################################################################
# 
### FUNCTIONS ############################################################################
# install libraries if not present
# TCGA.libs <- c("TCGAbiolinks","SummarizedExperiment");

# if (length(setdiff(TCGA.libs, rownames(installed.packages()))) > 0) {
#   source("http://bioconductor.org/biocLite.R")
#   biocLite(setdiff(TCGA.libs, rownames(installed.packages())))
# }
# ^ this uses setdiff to compare our libraries and see if they are installed
# it then installs it if they are missing
# using setdiff to do this is bad because setdiff will fail if anything else installed
# and if it's not installed it then installs 
# should use %in% instead 
# also install.packages doesn't work for all R versions
# but biocLite() function installs for all R versions


# this sets working directory initially to Desktop (softcoded this)
input = readline('Enter Desktop directory in full (case sensitive) e.g. /home/<username>/Desktop: ')
# e.g. can write /home/samimonk/Desktop
setwd(input)

# this creats a vector with two strings in it (the library names)
TCGA.libs <- c("TCGAbiolinks","SummarizedExperiment");

# this will identify if library missing from out installed libraries
new.libs <- TCGA.libs[!(TCGA.libs %in% installed.packages() [,"Package"])]

# this detects and installs the specific libraries if one or both are missing
if(length(new.libs)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite(new.libs)
}

# loads in the libraries to R
library("TCGAbiolinks")
library("SummarizedExperiment")

### MAIN #################################################################################
# define variables
ProjectDir <- "mutations"; # Currently in Desktop directory and want to work from /Desktop/mutations/
ProjectID <- "TCGA-BRCA"; # the dataset name as defined by TCGAbiolinks
ProjectName <-  "mutations";

# set/create the project directory i.e. "/Desktop/mutations/"
if (file.exists(ProjectDir)){

	cat( paste("\nFolder with project [", ProjectName, "] already exists. The requested data will be stored in [", ProjectDir, "]\n", sep=" ") );
# ^ this is just telling the user that the directory exists and we will be working in it 

} else {
	dir.create(ProjectDir);
# ^ this makes the directory if it didn't exist already (and is called mutations since we assigned that earlier)
}

# change working directory to the project workspace (mutations folder)
setwd(ProjectDir);


#===============================================================================
#    Data query and download
#===============================================================================
# ACCESS AND DOWNLOAD THE CLINICAL DATA
# this is a biolinks command to download the data, we assign it into the object clin.data
clin.data <- GDCquery_clinic( "TCGA-BRCA", "clinical" );


# this command just tells us the dimensions of the datafile
# 1097 number of rows, 39 collumns
dim(clin.data)

# Define type of tumour e.g. BRCA (breast)
tumour.type <-  "BRCA"

# Write clinical data to outfile - this is your target file
# write.table writes the data matrix to a file (called the target file)
# using "\t" creates a delimited output (the default is "")
# so this is a "tab delimited text file"
# row.names as FALSE means that theres no extra collumn with row names in our data matrix
# clin.data is the name of the object that represents the data matrix written into the file
# so we can handle the data from the object data matrix by using the object name "clin.data"

write.table( clin.data, file = paste(tumour.type,"_", ProjectName, "_clinical.txt", sep=""), sep="\t", row.names=FALSE );

# what clinical information do you want to appear in your figure?
# can view the options by seeing the collumn names of our data matrix object "clin.data"

colnames(clin.data);

# after checking this I identified the following as useful:

# race = useful states if "black or african american" or "white" or "asian"
# year_of_birth = very useful
# bcr_patient_barcode = very useful (but is it the same as submitter Id?)
# disease = it's all just BRCA, i guess it confirms you picked the right group since this will be soft coded
# gender = very useful, mostly female but some male (can see data trends possibly)
# tissue_or_organ_of_origin = seems a better indictator than site of biopsy (mainly c50.9)
# morphology = gives cancer morphology codes (e.g. 8520/3), seems seems quite useful
# vital_status = seems very useful, states if patient is alive or dead
# tumor_stage = very useful (stage i, ii, iii etc.)

# Subsetting the clinical data to covariates of interest
# we assign the collumn names to the object "clin.covariates"
clin.covariates <- c("race", "year_of_birth", "bcr_patient_barcode", "disease", "gender", "tissue_or_organ_of_origin", "morphology", "vital_status", "tumor_stage");
# now we produce a new data matrix with only the collumns that we choose, leaving the rows as they are
clin.data.slimmed <- clin.data[ , clin.covariates ];
# checking the dimensions of our new "slimmed down" data matrix to confirm that we only have 9 rows
dim(clin.data.slimmed) #1097 9

# we now have our slimmed down data matrix


# ACCESS AND DOWNLOAD THE MUTATION DATA
# overview of each pipeline
# this is a very useful paper on InDel calling software http://humgenomics.biomedcentral.com/articles/10.1186/s40246-015-0042-2
# [1] mutect2
# MuTect2 is a somatic SNP and indel caller that uses MuTect 1 cosomatic genotyping engine (Cibulskis et al., 2013) 
# with the assembly-based machinery of HaplotypeCaller. 
# The basic operation of MuTect2 proceeds similarly to that of the HaplotypeCaller.
# BUT MuTect2 doesn't use ploidy assumption but informs its genotype likelihood and variant quality calculations 
# by varying allelic fraction for each variant (as is often seen in tumors with purity less than 100%)
# and also for multiple subclones, and/or copy number variation (either local or aneuploidy). 

# MuTect2 applies some hard filters to variants before producing output.
# Gene Variant Call Formaat (GVCF) generation is not available in MuTect2.

# [2] varscan
Most of the published variant callers for next-generation sequencing data employ a probabilistic framework, 
such as Bayesian statistics, to detect variants and assess confidence in them. 
These approaches generally work quite well, but can be confounded by numerous factors such as extreme read depth, 
pooled samples, and contaminated or impure samples. 

In contrast, VarScan employs a robust heuristic/statistic approach to call variants that meet desired thresholds for read depth, base quality, variant allele frequency, and statistical significance.

VarScan is under continued development and improvement at a leading genome center with early access to new sequencing technologies, substantial computing resources, immense public/private datasets, and established expertise in sequencing, genetics, and genomics.

Detecting Subclonal Mutations
A 2013 study by Stead et al evaluated several somatic mutation callers including MuTect, Strelka, and VarScan2. They found that VarScan2 performed best overall with sequencing depths of 100x, 250x, 500x and 1000x required to accurately identify variants present at 10%, 5%, 2.5% and 1% respectively. 

# [3] muse


# [4] somaticsniper


### in your opinion, which method would be most applicable to breast/pancreatic cancer? i.e. look at literature

### should we give users the option to choose what method they want to apply?
# input = readline('Please key in the pipline of choice [1] mutect2, [2] varscan, [3] muse, or [4] somaticsniper: ')
#
# if input == 1
#  x
#
# else if input == 2

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
