## Instructions: 
## Variables:
## The logical variables in lines 35-36 are set to whether ancestry/phenotype PCA files should be generated. We recommend to set them at FALSE as these files are provided with the data, unless one desires to regenerate them (propably due to updates in the data/scripts). In such case, replace the condition "!file.exists(...)" with TRUE. 
## Similar considerations are to be taken with regard to generating updated GDS files and LD-pruning in the script "ancestry_pcs.R" andbam files in the script "generate_bam.sh" . For now, these files have been generated and provided in the data folder. 

## setup the environment
library(GENESIS)
library(GWASTools)
library(SNPRelate)
library(PCAtools)
library(tidyverse) 
library(magrittr)
library(here)
library(corral)
library(BiocSingular)
library(RNOmni)
library(vroom)
library(rtracklayer)

setwd(here()) # set up the main folder as directory, e.g. ~/user/.../pipe_repo 
dir.create(file.path('./data', 'pop_standardized_data'), showWarnings = F)
genes_all_pops <- NULL # stores intersection genes between all populations

## This script imports 
# a. Ancestry PCs
# b. Phenotype PCs
# c. Population-specific raw counts data
## Cross IDs in these files, rewrite the PCs into a single file where IDs (cols) order matches the order in the .fam file (expression matrix)
## Output: 
## ancestry PCs, phenotype PCs and GDS files (ONLY IF required/not provided)
## Combined covariates file, rearranged such that included samples re only those shared between both sets of PCs and with the gene-expression file. This  will be used as inlut for TensorQTL.
## Subset of the normalized gene-expression data, including samples shared with the covariates above. 
## Updated .bed file with name format "PopulationName_phenotype.bed" in ./data/ 
# Remeber: ----------------------------------------------------------------------
# Need only to correct for ancestry PCs for MKK; we found that for all other populations, non of the eigen values passed the Gavish Dunnoho threshold, except in MKK.    
populations <- c("ESN", "GWD", "LWK", "MSL", "YRI", "MKK")
counts_aggregated <- read_tsv('./data/counts_aggregated.txt')
for( pop_name in populations){
   
    print(pop_name) 
    
    source('./scripts/get_pop_data.R')
    
    source('./scripts/filter_normalize.R')
} 
# genes_all_pops contains the intersection genes across all 6 populations. 
## Note for later: make the script removr the directory 'pop_standarized_data' when done 



for( pop_name in populations){
    counts_population <- read_rds(str_c('./data/pop_standardized_data/',pop_name, 'counts_population.rds' ))

    keep_rows <- which(counts_population$Name %in% genes_all_pops)

    counts_population <- counts_population[keep_rows,]

    standardized_data <- read_rds(str_c('./data/pop_standardized_data/',pop_name, 'standardized_data.rds' ))[keep_rows,]
    
    compute_pheno_pcs <- T # !file.exists(str_c('./data/pheno_pcs/',pop_name,'_pcs_pheno.rds'))
    compute_ances_pcs <- !file.exists(str_c('./data/ancestry_pcs/',pop_name,'_pcs_ances.rds'))

    if (compute_pheno_pcs) source('./scripts/pheno_pcs.R')
    
    if (compute_ances_pcs) source('./scripts/ancestry_pcs.R')
    
   source( './scripts/cross_files.R')
}




