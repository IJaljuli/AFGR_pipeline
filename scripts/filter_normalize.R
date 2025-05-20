# This script is intended to import raw gene counts >> filter out genes >> normalize. Output: a normalized phenotype matrix. 
dir.create(file.path('./data', 'pop_standardized_data'), showWarnings = F)
counts_raw <- counts_population %>% select(-c('Name','Description')) %>% as.matrix

# Use two methods from scran to identify sufficiently variable genes for downstream analysis and subset the data.
# https://bioconductor.org/books/3.18/OSCA.basic/feature-selection.html#quantifying-per-gene-variation
# ONLY valid for integer counts data (Poisson, geometric, negative binomial etc.)
# Other omes that are proportions or log-normal do not have similarly prinicpled methods at this time
set.seed(12345L)
counts_gene_var <- scran::modelGeneVarByPoisson(counts_raw)
counts_gene_var_top <- scran::getTopHVGs(counts_gene_var)
counts_gene_cv2 <- scran::modelGeneCV2(counts_raw)
counts_gene_cv2_top <- scran::getTopHVGs(counts_gene_cv2, var.field = "ratio", var.threshold = 1L)
counts_gene_shared <- intersect(counts_gene_var_top, counts_gene_cv2_top)

counts_variable <- counts_raw[counts_gene_shared, ]

# Use corral to normalize counts for library depth and gene length, and use RNOmni for the inverse normal transformation (on the ranks).
# Use the standardized_data for QTL mapping directly.
# Also use corral for counts data
# For other proportions/lognormal data types just use inverse normal, because there are not "library sizes"
standardized_data <- corral::corral_preproc(counts_variable, rtype = "freemantukey") %>% apply(1, RNOmni::RankNorm) %>% t()
counts_population  <- data.frame(Name = counts_population$Name[counts_gene_shared],
                                Description = counts_population$Description[counts_gene_shared],
                                standardized_data)
write_rds(standardized_data ,str_c('./data/pop_standardized_data/',pop_name, 'standardized_data.rds' ))
write_rds(counts_population, str_c('./data/pop_standardized_data/',pop_name, 'counts_population.rds' ))
if(is.null(genes_all_pops)){
    genes_all_pops <- counts_population$Name
}else{
 genes_all_pops <- intersect(genes_all_pops, counts_population$Name)
}

dir.create(file.path('./data', 'pop_standardized_data_splicing'), showWarnings = F)
splice_junctions_population_raw <- read_tsv(file.path("./data/pop_splicing", str_c(pop_name, "_splice_junctions.tsv"))) %>%
    column_to_rownames("cluster") %>%
    as.matrix
standardized_data_splicing <- apply(1, RNOmni::RankNorm) %>% t
splice_junctions_population <- data.frame(cluster = rownames(splice_junctions_population_raw), standardized_data_splicing)
write_rds(standardized_data_splicing, str_c('./data/pop_standardized_data_splicing/', pop_name, 'standardized_data_splicing.rds' ))
