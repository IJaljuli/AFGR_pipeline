# This script is intended to import raw gene counts >> filter out genes >> normalize. Output: a normalized phenotype matrix. 
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
standardized_data <- corral::corral_preproc(counts_variable) %>% apply(1, RNOmni::RankNorm) %>% t()
counts_population  <- data.frame(Name = counts_population$Name[counts_gene_shared],
                                Description = counts_population$Description[counts_gene_shared],
                                standardized_data)




