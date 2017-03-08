#The Cancer Genome Atlas

This is the TCGA readme.

## How to use bioconductor TCGA biolinks to download and analyse data
http://bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/tcgaBiolinks.html


## Paper summaries

### TCGA Workflow: Analyze cancer genomics and epigenomics data using Bioconductor packages

https://f1000research.com/articles/5-1542/v1

### TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4856967/

TCGA is a large collection of clinical and molecular phenotypes of more than 10 000 tumor patients across 33 different tumor types. 

From this, TCGA has published over 20 marker papers detailing the genomic and epigenomic alterations associated with these tumor types. 

Implementing novel methods can elucidate new biological pathways and diagnostic markers. 

Mining the TCGA data presents several bioinformatics challenges, like data retrieval and integration with clinical data and other molecular data types (e.g. RNA and DNA methylation).

An R/Bioconductor package was developed (TCGAbiolinks) to address these challenges and offer bioinformatics solutions by using a guided workflow to allow users to query, download and perform integrative analyses of TCGA data.

Four different TCGA tumor types (Kidney, Brain, Breast and Colon) were ued as case studies to illustrate examples of reproducibility, integrative analysis and utilization of different Bioconductor packages to advance and accelerate novel discoveries.

# Cancer Genomics from Discovery science to Personalised Medicine

http://www.nature.com/nm/journal/v17/n3/full/nm.2323.html

development of genomics has turned cancer medicine from morphology-based to genetics-based taxonomy.
increasingly customised care for patient depending on proteomics and genomics of tumour.
personalised medicine is a maturing reality and no longer fantasy
level of genomic complexity is far beyond expectations: practical utilisation of genetic information is unrealised in medicine
2 main challenges causing slow down in translation from science to medicine:
1. incomplete catalogue of genomic alterations
2. lack of understanding of how alterations affect pathways / signal transduction / what the biological consequences are etc.
- hard to obtain this as specific assays are not always available + cancer context is very complex and hard to identify which driver events cause what (as intertwined)...

this paper highlights importance of biological significance and how personalised medicine is achivable.

first cancer mutation discovered HRAS gene gly12 to val12, then others (KRAS, NRAS) discovered: this ushered in cancer genomics study

RAS proteins are important for transduction of proliferation and survival
- after 30 years of trying to inhibit it unsuccessfully, only NOW is being used as some kind of marker...
mutant cancer genes can be targets for therapy

most early efforts for cancer genetics research focused on druggable targets (protein/lipid kinases)
V600 mutation linked to BRAF mutation (kinase central to MAPK cascade, transduces growth factors via RAS proteins to phosphorylate MEK proteins)
BRAF inhibition was developed after 8 years; effective in tumours with V600E mutation
PIK3 is most commonly mutated gene in breast cancer: drives AKT signalling for growth and stability


# initial investigation

Obtain mutations data from TCGA (breast or pancreatic)

play around with R analysis

need to subset data to obtain top hit dif gene expression + genes of interest to end-user

## using this tutorial

https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html

