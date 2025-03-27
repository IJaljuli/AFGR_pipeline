# Read part I instructions
# In this script we create GDS files for faster processing of the bam files. 
## Check if GDS/ancestry_pcs directories exsist, create if not and generate missing files  
#get_bam <- file.exists(str_c("./bed_bim_fam/all", pop_name, "30x.ID.bed", sep = ".") ) # needs to be moved to sourcing file 
#if( !get_bam ){
#}


my.bam.path <- './data/bed_bim_fam/'
my.PCs.path <- './data/ancestry_pcs/'
# pop_name must be provided in sourcing script



dir.create(file.path('./data', 'gds'), showWarnings = F)
(get_gds <- !file.exists(str_c("./data/gds/collapsed.hwepass01.withgnomad.var.merged.processedMKK.split1kG30x", pop_name, "maf05.geno01.gds", sep = ".")))
get_gds <- T

dir.create(file.path('./data', 'ancestry_pcs'), showWarnings = F)
(do_ldpruning <- !(file.exists('./data/ancestry_pcs/pruned_variants.rds') | file.exists(str_c('./ancestry_pcs/', pop_name,  "_pruned_variants.rds", sep = '')) ))
do_ldpruning <- T


## Functions -------------------------------------------------------------------------
make_gds <- function(population_name) {
    bed_name <- str_c("collapsed.hwepass01.withgnomad.var.merged.processedMKK.split1kG30x", population_name, "maf05.geno01.bed", sep = ".")
    bed_name <- str_c(my.bam.path,bed_name)
   
    bim_name <- str_c("collapsed.hwepass01.withgnomad.var.merged.processedMKK.split1kG30x", population_name, "maf05.geno01.bim", sep = ".")
    bim_name <- str_c(my.bam.path, bim_name) 

    fam_name <- str_c("collapsed.hwepass01.withgnomad.var.merged.processedMKK.split1kG30x", population_name, "maf05.geno01.fam", sep = ".") 
    fam_name <- str_c(my.bam.path, fam_name) 
    
    gds_name <- str_c("./data/gds/collapsed.hwepass01.withgnomad.var.merged.processedMKK.split1kG30x", population_name, "maf05.geno01.gds", sep = ".")
    snpgdsBED2GDS(bed.fn = bed_name, 
                  bim.fn = bim_name, 
                  fam.fn = fam_name, 
                  out.gdsfn = gds_name)
}


ld_prune <- function(population_name) {
    gds_name <- str_c("./data/gds/collapsed.hwepass01.withgnomad.var.merged.processedMKK.split1kG30x", population_name, "maf05.geno01.gds", sep = ".")
    gds_file <<- snpgdsOpen(gds_name)
    pruned_snps <- snpgdsLDpruning(gds_file, method = "corr", slide.max.bp = 10e6, ld.threshold = sqrt(0.1), verbose = TRUE, num.thread = 1L)
    snps_vector <- unlist(pruned_snps, use.names = FALSE) 
    snps_vector
}

choose_npcs <- function(pca_estimate) {
    n_pcs <- chooseGavishDonoho(.dim = c(length(pca_estimate$eigenval), length(pca_estimate$snp.id)), var.explained = pca_estimate$eigenval, noise = 1)
    n_pcs
}

genotype_pcs <- function(population_name, pruned_variants) {
    gds_name <- str_c("./data/gds/collapsed.hwepass01.withgnomad.var.merged.processedMKK.split1kG30x", population_name, "maf05.geno01.gds", sep = ".")
    if (is.null(gds_file)) {
        gds_file <- snpgdsOpen(gds_name)
    }
    gds_pca <- snpgdsPCA(gds_file, snp.id = pruned_variants[[population_name]], algorithm = "exact", eigen.method = "DSPEV", eigen.cnt = -1)
    snpgdsClose(gds_file)
    n_pcs <- choose_npcs(gds_pca)
    sample_names <- gds_pca$sample.id
    # If the number of select PCs is 1, it may actually be 0 or 1, so we need to manually check how many eigenvalues 
    if (n_pcs == 1) {
        eigenval_limit <- attr(n_pcs, "limit")
        n_pcs <- sum(gds_pca$eigenval > eigenval_limit)
    }
    if(n_pcs > 0){
       gds_pca <-  gds_pca$eigenvect[,1:n_pcs] %>% t %>% as.data.frame
       colnames(gds_pca)<- sample_names
       return(gds_pca)
    }else{
        return(NULL)
    }
}

extract_eigenvalues <- function(pca_estimate, population_name) {
    num_eigenval <- length(pca_estimate$eigenval)
    tibble(population = population_name, number = 1:(num_eigenval-1), eigenvals = pca_estimate$eigenval[-num_eigenval])
}


## Analysis ------------------------------------------------------------------------------------------------------------------
### Part I: This part creates GDS files, performs LD-pruning on the genotype matrices. 
### These files need to be produced once prior to analysis. Due to the long time to produce such files, they are provided in the main directory. Do not run these lines unless the described files were not initially provided or the provided files were changed.

## Create GDS files for all populations
gds_file <- NULL
if( get_gds){
    # if (pop_name == 'MKK'){
    #     snpgdsBED2GDS(bed.fn = str_c(my.bam.path,"/MKK.filtered.biallelic.maf05.autosomal.chr.id.bed"),
    #               bim.fn = str_c(my.bam.path,"/MKK.filtered.biallelic.maf05.autosomal.chr.id.bim"),
    #               fam.fn =str_c(my.bam.path,"/MKK.filtered.biallelic.maf05.autosomal.chr.id.fam"),
    #               out.gdsfn = "./data/gds/all.MKK.30x.ID.gds")
    # }else{
        make_gds(pop_name)
    # }
}

# Variant pruning: Pruned variants are variants which are independent under LD up to a specific LD cutoff (here we allow pairwise LD < sqrt(0.1).
if( do_ldpruning ){
    message('LD prining in process')
    pruned_variants <- list(ld_prune(pop_name))
    names(pruned_variants) <- pop_name # populations
    write_rds(pruned_variants,str_c(pop_name,  "pruned_variants.rds", sep = '_'))
}else{
    if(file.exists("./data/ancestry_pcs/pruned_variants.rds")){
         pruned_variants <- read_rds("./data/ancestry_pcs/pruned_variants.rds")
    }else{
         pruned_variants <- read_rds(str_c("./data/ancestry_pcs/", pop_name, "_pruned_variants.rds"))
    }
}

# Find SVD decomposition >>  select top npcs eigenvectors based on their eigen values where npcs is calculated via the Gavish-Dunnoho formula.
(pop_pca <- genotype_pcs(pop_name, pruned_variants))

# if( pop_name == 'MKK'){
#         id_key_path <- './data/MKK.dictionary.txt' #"./data/all_MKK_166_newIDS.csv"
#         id_key <- read_tsv(id_key_path, col_names = FALSE) %>% drop_na()
#         txt1 <- paste(id_key$X1, id_key$X2, sep = "' ~ '") # creating coding in txt1 
#         txt1 <- paste0("'",txt1)
#         txt1 <- paste0(txt1,"'")
#         txt2 <- paste0(c('case_match(colnames(pop_pca)', txt1, ' .default = colnames(pop_pca))'), collapse = ', ')
#         colnames(pop_pca) <- eval(parse(text=txt2)) 
# }

if(!is.null(pop_pca)){
    pop_pca <- cbind(variable = paste0('ancestry_PC', 1:nrow(pop_pca)), pop_pca)
}

write_rds(pop_pca, str_c('./data/ancestry_pcs/',pop_name, "_pcs_ances.rds"))

