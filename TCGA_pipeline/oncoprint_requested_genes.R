# ONCOPRINT OF USER REQUESTED GENES

# set wkdir to my personal Desktop
setwd(/home/samimonk/Desktop)

# detect and install missing libraries
TCGA.libs <- c("TCGAbiolinks","SummarizedExperiment");
new.libs <- TCGA.libs[!(TCGA.libs %in% installed.packages() [,"Package"])]
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

# create mutations directory and go into it
if (file.exists(ProjectDir)){
  cat( paste("\nFolder with project [", ProjectName, "] already exists. The requested data will be stored in [", ProjectDir, "]\n", sep=" ") );
} else {
  dir.create(ProjectDir);
}
setwd(ProjectDir);

# Download clinical data
clin.data <- GDCquery_clinic( "TCGA-BRCA", "clinical" );
dim(clin.data)
tumor.type <-  "BRCA"
write.table( clin.data, file = paste(tumor.type,"_", ProjectName, "_clinical.txt", sep=""), sep="\t", row.names=FALSE );
clin.covariates.for.oncodata <- c("race", "gender", "vital_status", "tumor_stage", "bcr_patient_barcode");

#setting up display of clinical data for oncoprint
clin.forvisual <- clin.data[ , clin.covariates.for.oncodata]
clin.data.slimmed <- clin.data[ , clin.covariates ];


######## Download mutations data & visualising clinical and mutations data in an oncoprint

# oncoprint of user requested genes
user.gene.request = readline('Please type the genes of interest in capital letters, separated by " ; " with NO SPACE (e.g. " ATM;BRCA1;CHEK2 "): ')
# COMMONLY MUTATED BREAST CANCER GENES
# ATM;p53;BRCA1;BRCA2;PTEN;CHEK2;PALB2;STK11;BARD1;BRIP1;CASP8;CDH1;CHEK2
user.gene.list <- unlist(strsplit(user.gene.request,";"))

pipeline_options <- c("muse", "varscan2", "somaticsniper", "mutect2");

for(pipe in pipeline_options){
  mut.data <- GDCquery_Maf(tumor = tumor.type, pipelines = "mutect2", save.csv = TRUE);
  mut.data.matrix <- matrix(data=0, nrow=length(user.gene.list), ncol=1, dimnames=list(rownames=user.gene.list, colnames="number_of_reports"));
  all.gene.positions <- NULL
  for(gene in user.gene.list) {
    user.gene.list.items <- which(mut.data$Hugo_Symbol %in% user.gene.list)
    mut.data.matrix[gene, "number_of_reports"] <- length(user.gene.list.items)
    all.gene.positions <- c(all.gene.positions, user.gene.list.items)
  }
  mut.data.matrix.ordered <- mut.data.matrix[ order(mut.data.matrix[, "number_of_reports"], decreasing=TRUE),];
  user.choice.gene.names <- names( mut.data.matrix.ordered[ 1:20 ] );
  user.choice.mut.data <- mut.data[ all.gene.positions, ];
  for(pipe in pipeline_options){
    TCGAvisualize_oncoprint(
        mut = user.choice.mut.data,
        genes = user.choice.gene.names,
        filename = paste("oncoprint_", tumor.type, user.gene.request, "_", pipe, ".pdf", sep=""),
        annotation = clin.forvisual,
        color=c("background"="#CCCCCC","DEL"="purple","INS"="yellow","SNP"="brown"),
        rows.font.size= 8,
        width = 5,
        heatmap.legend.side = "right",
        dist.col = 0,
        label.font.size = 6
      );
  }

  dev.off();

}
