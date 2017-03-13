#The Cancer Genome Atlas

This is the TCGA readme.

## How to use bioconductor TCGA biolinks to download and analyse data
http://bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/tcgaBiolinks.html
## another useful tutorial
https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html

# initial investigation

Obtain mutations data from TCGA (breast or pancreatic)

play around with R analysis

need to subset data to obtain top hit dif gene expression + genes of interest to end-user

## Identifying useful subsets
## all the subsets available are:

### useful
- race = useful states if "black or african american" or "white" or "asian"
- year_of_birth = very useful
- bcr_patient_barcode = very useful (but is it the same as submitter Id?)
- disease = it's all just BRCA, i guess it confirms you picked the right group since this will be soft coded
- gender = very useful, mostly female but some male (can see data trends possibly)
- tissue_or_organ_of_origin = seems a better indictator than site of biopsy (mainly c50.9)
- morphology = gives cancer morphology codes (e.g. 8520/3), seems seems quite useful
- vital_status = seems very useful, states if patient is alive or dead
- tumor_stage = very useful (stage i, ii, iii etc.)


### Not useful
- submitter_id = appears to be identical to the patient barcode value
- last_known_disease_status = all not reported
- classification_of_tumor = not reported
- primary_diagnosis = not sure why these are all c50.9 since tissue of origin is varied
- updated_datetime = seems irrelevent
- age_at_diagnosis = seems to be not useful in this context
- days_to_death = seems to actually be days SINCE death, mostly NA, not useful
- days_to_last_follow_up = seems a bit irrelevent especially since many are dead
- state = all NA
- days_to_last_known_disease_status = all NA
- days_to_recurrence = all NA
- tumour grade = all not reported
- diagnosis_id = seems a bit irrelevent
- treatments (Ids) = patient barcode with another code, seems a bit irrelevent
- site_of_resection_or_biopsy = possibly useful (majority are c50.9) but tissue of origin is better
- days_to_birth = irrelevent
- ethnicity = it's all either "(not) hispanic or latino" or "not reported"
- year_of_death = kind of useful
- demographic_id = seems a bit irrelevent
- progression_or_recurrence = all not reported
- prior_malignancy = all not reported
- created_datetime = all NA
- alcohol_intensity = it's a shame this is all NA as would be relevent
- alcohol_history = it's a shame this is all NA as it would be relevent
- bmi = it's a shame this is all NA as would be relevent
- weight = it's a shame this is all NA as would be relevent
- cigarettes_per_day = it's a shame this is all NA as would be relevent
- height = all NA so useless here
- years_smoked = all NA so also useless for BRCA data
- exposure_id = unsure of the significance of this

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
