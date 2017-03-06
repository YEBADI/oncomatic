# BCCTBbp: Breast Cancer Campaign Tissue Bank bioinformatics portal

This is a readme of BCNTB.

## Paper Summaries

### BCCTBbp: Breast Cancer Campaign Tissue Bank bioinformatics portal

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4384036/pdf/gku984.pdf

BCCTBbp  is  dedicated  to  maximising  research  on  patient  tissues

- stores  genomics,  methylomics,  transcriptomics,  proteomics  and  microRNA  data

- these have been mined from literature
- allows researchers to pose new questions and link them to  pathways  and mechanisms  involved  in  breast  cancer.

Currently, the  portal  holds  146  datasets and contains over  227, 795  expression genomic  measurements  from  various breast tissues (e.g. normal, malignant or benign lesions), cell lines and body fluids.

- Has a dedicated analytical layer allowing researchers to further analyse stored datasets

Aim of this project is to work on this aspect and further develop and enhance the integrated analytical layers to facilitate tailored analysis of key datasets, thus allowing users to pose questions not addressed in the original publications simply and efficiently.


### Analytical Options

Four analysis categories can be performed: 

#### Molecular Classification
Each of the included dataset samples can be classified based on
- the PAM50 set of markers
- OR the hormonal receptor status can be inferred from the transcriptome using the MCLUST R package (http://cran.r-project.org/web/packages/mclust/)
#### Tumour Purity
Cancer samples frequently contain a small proportion of normal adjacent tissue that might confuse sample analysis. A method to infer sample cancer purity is implemented using ESTIMATE
#### Gene Expression
Gene expression plots can be obtained for a gene of interest
#### Survival
For datasets that contain patient survival information, survival analysis can be performed using the ‘survival’ package in R (http://cran.r-project.org/web/packages/survival/)

Two modes are available
- survival based on sample sub-groups within the dataset
- OR survival based on the expression value of a gene of interest

In the gene-based survival analysis, two classes of gene expression values (high and low) are determined based on the median expression value.
