# This script is intended to import normalized phenotype data "standardized_data" via "filter_normalize.R", and calculate phenotype PCs.
dir.create(file.path('./data', 'pheno_pcs'), showWarnings = F)

# Run PCA for hidden factor estimation.
# Use the Gavish-Donoho threshold to choose the number of PCs.
# Alternatively, can use chooseMarcenkoPastur() for the less conservative Marcenko-Pastur method.
pca_standardized <- PCAtools::pca(standardized_data, BSPARAM = BiocSingular::ExactParam())
(n_pcs <- PCAtools::chooseGavishDonoho(standardized_data, var.explained = pca_standardized$sdev^2, noise = 1))

# Return to this after dan responds ------------------------------------------------------------------------------------
## If the number of select PCs is 1, it may actually be 0 or 1, so we need to manually check how many eigenvalues 
 if (n_pcs == 1) {
     eigenval_limit <- attr(n_pcs, "limit")
     n_pcs <- sum((pca_standardized$sdev^2) > eigenval_limit)
 }

pca_eigenvectors <- NULL

if( n_pcs > 0 ){
# You will probably use this object for your hidden factors.
pca_eigenvectors <- pca_standardized$rotated[,1:n_pcs] %>% t %>% as.data.frame
## Optionally, run rank normal transformation on PCs - this is not normally done but is provided here as an example.
## Iman doesn't like this step, will be using pca_eigenvectors object. 
# pca_eigenvectors_standardized <- apply(pca_standardized$rotated[,1:n_pcs], 2, RNOmni::RankNorm)}
pca_eigenvectors <- cbind(variable = paste0('phenotype_PC', 1:n_pcs), pca_eigenvectors)
}

write_rds(pca_eigenvectors, str_c('./data/pheno_pcs/',pop_name, '_pcs_pheno.rds'))
