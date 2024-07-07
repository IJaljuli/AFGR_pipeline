# This script is meant to be sourced to get the data of a selected population 
# 'pop' is a variable given at source
pop_counts <- function(population_name){ 
    if(population_name != 'MKK'){
        fam_name <- str_c("./data/bed_bim_fam/all", population_name, "30x.ID.fam", sep = ".")
    }else{
        fam_name <- './data/bed_bim_fam/MKK.filtered.biallelic.maf05.autosomal.chr.id.fam'
    }
    fam_file <- read_delim(fam_name, col_names = FALSE)
    fam_ids <- fam_file$X2
    if( population_name=='MKK'){
        id_key_path <- "./data/all_MKK_166_newIDS.csv"
        id_key <- read_csv(id_key_path, col_names = FALSE) %>% filter(is.na(X3)) %>% select(-X3)
        txt1 <- paste(id_key$X1, id_key$X2, sep = "' ~ '") # creating coding in txt1 
        txt1 <- paste0("'",txt1)
        txt1 <- paste0(txt1,"'")
        txt2 <- paste0(c('case_match(fam_ids', txt1, ' .default = fam_ids)'), collapse = ', ')
        fam_ids <- eval(parse(text=txt2))
    }
    counts_ids <- colnames(counts_aggregated)[-c(1,2)]
    shared_samples <-  intersect(fam_ids,counts_ids)
    counts_aggregated[,c(TRUE, TRUE, counts_ids %in% shared_samples )]
}

#counts_aggregated <- read_tsv('./data/counts_aggregated.txt')
counts_population <- pop_counts(pop_name)
