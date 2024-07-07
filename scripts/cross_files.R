#library(SNPRelate)
#library(PCAtools)
#library(tidyverse) 
#library(magrittr)
#library(here)
##library(readr)
#library(corral)
##library(scran)
#library(BiocSingular)
#library(RNOmni)
#library(vroom)
#library(rtracklayer)

#setwd(here()) # set up the main folder as directory, e.g. ~/user/.../pipe_repo
#pop_name <- 'GWD'
#source('./scripts/get_pop_data.R')
#source('./scripts/filter_normalize.R')
#source('./scripts/pheno_pcs.R')
#source('./scripts/ancestry_pcs.R')



dir.create(file.path('./data', 'covariates'), showWarnings = F)

# At this stage, the file ./scripts/get_pop_data.R' must have been sourced and object "counts_population" has been created. This object is the intersection of .fam and raw counts for id matching where the columns have been filtered and normalized.
# Now we filter gene annotation GTF down to only gene entries - this is the same for every population
annotation_gtf <- readGFF("./data/Homo_sapiens.GRCh38.110.collapse.gtf")
annotation_genes <- filter(annotation_gtf, type == "gene") %>%
    select(seqid, start, gene_id)
annotation_genes$end <- annotation_genes$start + 1
annotation_genes %<>% select(seqid, start, end, gene_id) # Use %<>% as the double pipe to send an object to a command and replace it with the output
annotation_genes <- annotation_genes %>% rename( "seqid" = "#chr", "gene_id" = "phenotype_id")

#  Note: counts_population is already crossed with the .fam file. Now we need to find the samples shared with the pcs files.
process_pop <- function(population_name, population_counts) {
   (pheno_pcs_file <- read_rds(str_c('./data/pheno_pcs/',population_name, '_pcs_pheno.rds')))
   (ances_pcs_file <- read_rds(str_c('./data/ancestry_pcs/',population_name, '_pcs_ances.rds')) )
 
    shared_columns_ances <- intersect(colnames(counts_population), colnames(ances_pcs_file))
    shared_columns_pheno <- intersect(colnames(counts_population), colnames(pheno_pcs_file))
    if(!any(c(is.null(pheno_pcs_file), is.null(ances_pcs_file))) ){
        (shared_columns <- intersect(shared_columns_pheno, shared_columns_ances))
    }else{
        if(is.null(pheno_pcs_file) ){
            shared_columns <- shared_columns_ances
        }else{
            shared_columns <- shared_columns_pheno
        }
    }
   
    counts_ids <- colnames(counts_population %>% select(-c('Name','Description')))
    counts_population <- counts_population[,c(TRUE, TRUE, counts_ids %in% shared_columns)]  
    # The following lines are to make sure that the selectin is done in the right order, i.e. covariates file's columns match the counts_population order.
    all_pcs_file <- rbind(pheno_pcs_file[,c('variable', colnames(counts_population %>% select(-c('Name','Description'))))],
                          ances_pcs_file[,c('variable',colnames(counts_population %>% select(-c('Name','Description'))))]) %>%
                    as.data.frame
    write_tsv(all_pcs_file, str_c('./data/covariates/',population_name, "_covariates.txt") )
    

    counts_population <- counts_population %>%
        select(-c('Description'))%>%
        rename("Name" = "phenotype_id" )
    # Join the gene annotoation
    phenotype_bed <- inner_join(annotation_genes, counts_population) %>% as_tibble %>% arrange(`#chr`, start) 
    write_tsv(phenotype_bed, str_c('./data/bed_bim_fam/',population_name, "_phenotype.bed"))
}

process_pop(pop_name,counts_population)
