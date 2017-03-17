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
# this is a biolinks command to download the general clinical data for all BRCA samples
# we assign this data it into the object clin.data
clin.data <- GDCquery_clinic( "TCGA-BRCA", "clinical" );


# this command just tells us the dimensions of the datafile
# 1097 number of rows, 39 collumns
dim(clin.data)

# Define type of tumor e.g. BRCA (breast)
tumor.type <-  "BRCA"

# Write clinical data to outfile - this is your target file
# write.table writes the data matrix to a file (called the target file)
# using "\t" creates a delimited output (the default is "")
# so this is a "tab delimited text file"
# row.names as FALSE means that theres no extra collumn with row names in our data matrix
# clin.data is the name of the object that represents the data matrix written into the file
# so we can handle the data from the object data matrix by using the object name "clin.data"

write.table( clin.data, file = paste(tumor.type,"_", ProjectName, "_clinical.txt", sep=""), sep="\t", row.names=FALSE );

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

# we now have our slimmed down clinical data matrix


# ACCESS AND DOWNLOAD THE MUTATION DATA
# now that we have the clinical data for the BRCA samples, we need to obtain the mutation data

# overview of each pipeline
# The TCGAbiolinks package has four inbuilt InDel pipelines provided
# below, we provide a summary of each one to inform and let users pick which
# pipeline will be most effective for their research purposes.
# all pipelines generate VCF output

# [1] MuTect2
# * Applies Bayesian classifier to detect somatic mutations with very low allele fractions
# requires only a few supporting reads however is confounded by low tumor purity.
# * Has high sensitivity and calls mutations with allelic fractions as low as 0.1 and below.
# * Low specifiticity.
# * Applies carefully tuned hard filters to compensate partially specificity issues.

# [2] Varscan2
# * Is generally outperformed by MuSE and MuTect2.
# * Low sensitivity and fails to pick up somatic SNVs of low allelic fraction
# as supresses mutations below allelic threshold.
# * Sensitivity can be improved but this drastically drops specificitiy and returns high levels of false positives.
# * Outperformed MuTect to identify variants present at 10%, 5%, 2.5% and 1%
# at sequencing depths of 100x, 250x, 500x and 1000x respectively (Stead et al, 2013).

# * However, Varscan2 circumvents confounding factor of tumor purity and extreme read depth
# as does not use probabilistic framework to detect variants and assess confidence in them
# but uses a robust heuristic/statistic approach to call variants that meet desired thresholds 
# for read depth, base quality, variant allele frequency, and statistical significance.

# [3] MuSE
# * Has outperformed MuTect2 on calling variants from ACC TCGA data.
# * Maximum likelihood or the Markov chain Monte Carlo (MCMC) method estimates Model parameters. 
# * Variants are classified into somatic, germ-line, and reversal to the homozygous reference
# by comparing the somatic variant allele fraction (π) between the paired tumor–normal samples.
# * A sample-specific error model to account for tumor heterogeneity and identify cutoffs is built.
# * Filters reduce the number of false positives by considering the sequence context surrounding the point mutations. 


# [4] SomaticSniper
# Returns high level of false positives and many of these are not in agreement with other InDel callers.
# Most of the literature shows it is outperformed by most other variant callers.

# Verdict
#* TCGA-BRCA analysis is suitable with either Mutect2 or MuSE. Ideally both should be run and the results correlated / contrasted.
#* TCGA-PAAD analysis would only be advisable with Varscan2 due to significantly low tumor purity (circa 40%).

# GIVING USERS ABILITY TO CHOOSE PIPELINE OPTION
### To give users the option to choose what method they want to apply:
# input = readline('Please key in the pipline of choice [1] mutect2, [2] varscan, [3] muse, or [4] somaticsniper: ')
#
# if input == 1
#  x
#
# else if input == 2

#establishing the pipline options object
pipeline_options <- c("muse", "varscan2", "somaticsniper", "mutect2");

# for each methods, create a file of all the mutations
# I will consider using R specific vectorised functions instead of "for" loops
# as these run faster. For example:
# vapply(pipeline_options, function(pipe)){
# I've tried to do this but haven't gotten it to work yet

for(pipe in pipeline_options){

	# for each method (pipe) obtain the mutation data
  # we create an object called mut.data which downloads the mutation data for each pipeline option
  # tumor is picking BRCA (as we defined tumor.type to be BRCA)
  # save.csv means it is saved as a csv format
  # pipelines = pipe, means for each "pipeline option" it will download the whole BRCA data
	mut.data <- GDCquery_Maf(tumor = tumor.type, pipelines = pipe, save.csv = TRUE);

	# identify all genes (symbol) reported as mutated
  # unique means it will return each gene name WITHOUT ANY OF THE DUPLICATES
  # so to clarify Hugo_symbol has all of the genes including duplicates whereas gene.names only has the unique genes
  # hugo symbol means the standardised gene name by human genome organsaiton (HUGO)
	genes.names <- unique(mut.data$Hugo_Symbol);

	# create an empty data matrix to serve as our "input data matrix"
  # i.e. this is where we will collect information about the most reported genes
  # when we run the analysis pipline on the BRCA dataset we just downloaded for each pipline
  # now we just have a list of gene names and then one column on its right with the "number_of_reports" listed in it
  # there is no title for the gene names and this isn't 'counted' as a column
  # so ncol = 1 means theres one column (number of reports) to the right of the names
	all.mut <- matrix(data=0, nrow=length(genes.names), ncol=1, dimnames=list(rownames=genes.names, colnames="number_of_reports"));

	# identify position of each gene and populate the matrix
  # first, in order to assign the gene position values into a list
  # we make the object all.positions
  # then when we get the index values (in gene.positions) it comes in a non-list format (what is this called?)
  # so we assign it to our all.positions obejct to make it into a list
	all.positions <- NULL;

  # now we've defined gene as the individual parts of gene.names
  # now we are collecting the indicies of where the unique genes (from gene.names)  are within the total list of genes (from hugo_symbol)
  # the which function returns indicides for the events where the event is TRUE
  # as in, where the gene name in hugo_symbol equals the unique gene name, we get an index value
  # then we define the length of our "gene.names" and "number of reports" columns to that of the length of the indicies list we just obtained
  for(gene in genes.names) {
		gene.position <- which(mut.data$Hugo_Symbol == gene);
		all.mut[gene, "number_of_reports"] <- length(gene.position);
		all.positions <- c(all.positions, gene.position );
	}

	# order the whole all.mut matrix by the most reported gene highest to lowest
	all.mut.ordered <- all.mut[ order(all.mut[, "number_of_reports"], decreasing=TRUE),];

	# We isolate the character strings (i.e. the gene "names") from the top 20 genes
  # these will be used in the visualisation process
	top.mut.genes <- names( all.mut.ordered[ 1:20 ] );

  # Subset mutation data; we are extracting only the genes we indexed from the whole mut.data vector
  # (question: are we losing any data since we are removing duplicates?)
	top.mut.data <- mut.data[ all.positions, ];





################################## I couldn't follow what's specifically happening from here ###################



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
		filename = paste("oncoprint_", tumor.type, "_", pipe,".pdf", sep=""),
		annotation = top.clin.data,
		color = c("background"="#CCCCCC","DEL"="purple","INS"="yellow","SNP"="brown"),
		rows.font.size = 8,
		width = 5,
		heatmap.legend.side = "right",
		dist.col = 0,
		label.font.size = 6
	);



user.gene.request = readline('Please type the genes of interest in capital letters, separated by a comma (e.g. ATM, BRCA1 ): ')
user.gene.list <- user.gene.request(as.vector(unlist(strsplit(str,",")),mode="list"))

  TCGAvisualize_oncoprint(
    mut = top.mut.data,
    genes = top.mut.data$user.gene.list, 
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



################################## I couldn't follow what's specifically happening to here ###################


# COMMONLY MUTATED BREAST CANCER GENES
#    ATM.
#    p53
#    BRCA1.
#    BRCA2.
#    PTEN
#    CHEK2
#    PALB2
#    STK11
#    BARD1.
#    BRIP1.
#    CASP8.
#    CDH1.
#    CHEK2.

## [Yusef: now imagine that I am a researcher interested in specific genes]
## Look up commonly mutated genes in breast cancer and amend the code so that it:
### reads in a user-supplied gene list
### generates the plot for these genes of interest