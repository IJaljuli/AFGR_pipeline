# This script is meant to be sourced to get the data of a selected population 
# 'pop' is a variable given at source
pop_counts <- function(population_name){
    # the following is not neaded after the latest update. 
    # if(population_name != 'MKK'){
    #     fam_name <- str_c("./data/bed_bim_fam/collapsed.hwepass01.withgnomad.var.merged.processedMKK.split1kG30x", population_name, "maf05.geno01.fam", sep = ".")
    # }else{
    #    fam_name <- './data/bed_bim_fam/MKK.filtered.biallelic.maf05.autosomal.chr.id.fam'
    #}
    
    fam_name <- str_c("./data/bed_bim_fam/collapsed.hwepass01.withgnomad.var.merged.processedMKK.split1kG30x", population_name, "maf05.geno01.fam", sep = ".")
    
   
    
    fam_file <- read_delim(fam_name, col_names = FALSE)
    fam_ids <- fam_file$X2
    # the following is not neaded after the latest update. 
    # if( population_name=='MKK'){
    #     id_key_path <- "./data/MKK.dictionary.txt"
    #     id_key <- read_tsv(id_key_path, col_names = FALSE) %>% drop_na()
    #     txt1 <- paste(id_key$X1, id_key$X2, sep = "' ~ '") # creating coding in txt1 
    #     txt1 <- paste0("'",txt1)
    #     txt1 <- paste0(txt1,"'")
    #     txt2 <- paste0(c('case_match(fam_ids', txt1, ' .default = fam_ids)'), collapse = ', ')
    #     fam_ids <- eval(parse(text=txt2))
    # }
    counts_ids <- colnames(counts_aggregated)[-c(1,2)]
    shared_samples <-  intersect(fam_ids,counts_ids)
    counts_aggregated_selected <-  counts_aggregated[,c(TRUE, TRUE, counts_ids %in% shared_samples )]
   return(counts_aggregated_selected)
}

counts_aggregated <- read_tsv('./data/counts_aggregated_GM2NA.txt')

counts_population <- pop_counts(pop_name)
