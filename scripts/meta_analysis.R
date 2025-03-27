setwd('/labs/smontgom/grps_smontgom/jaljuli/afgr_eqtls/pipe_repo/scripts')

# Import tensorQTL results, combine in one data set.  
library(tidyverse) 
library(meta)
library(metarep)

population_names <- c('YRI','MSL','LWK', 'ESN','GWD', 'MKK' )
all_pop_data <- NULL

## check if phenotypes are indeed shared across all populations: 
all_genes <- NULL
for(pop_name in population_names){
  tmp <- read_tsv(str_c('/labs/smontgom/grps_smontgom/jaljuli/afgr_eqtls/pipe_repo/data/bed_bim_fam/',
                        pop_name, '_phenotype.bed') ) %>% 
    select(phenotype_id) %>% distinct(.) %>% mutate(population = pop_name) 
  all_genes <- rbind(tmp, all_genes)
}

tmp <- all_genes %>% group_by(phenotype_id) %>% summarise(n_pops = n())

# all selected phenotypes are reported in the 6 populations. However, I later find out that the tensorQTL outputs do not report them for all populations. 
## I suspect that only significant p-values were outputed. 
table(tmp$n_pops)

all_genes <- tmp$phenotype_id


### random-effects meta analysis:



# ! MKK is not available yet!!
population_names <- c('YRI','MSL','LWK', 'ESN','GWD', 'MKK' )
all_pop_data <- NULL

for(pop_name in population_names){
  path_name <- str_c("../output/cis/",
                     pop_name, "_.cis_qtl.txt.gz")
  all_pop_data <-  read.table( gzfile(path_name) , sep="\t",header = T) %>%
    mutate(population = pop_name) %>%
    rbind(., all_pop_data)
}
dim(all_pop_data)

tensor_output_genes <- all_pop_data %>% pull(phenotype_id) %>% unique()
# Genes not outputed by tensorqtl
sum( ! (all_genes %in% tensor_output_genes))

all( tensor_output_genes %in% all_genes)


keep_anal <- all_pop_data %>% 
  group_by(variant_id, phenotype_id) %>% 
  summarise(n=n()) %>% 
  filter(n==5)

dim(keep_anal)

all_pop_data  <- all_pop_data %>%
  right_join(keep_anal)

tmp <- all_pop_data %>% 
  group_by(variant_id, phenotype_id) %>% 
  summarise(n=n()) %>% 
  group_by(phenotype_id) %>% 
  summarise(n_max = max(n))

summary(all_pop_data$pval_beta)
summary(all_pop_data$pval_nominal)



all_pop_data  <- all_pop_data %>%
  right_join(keep_anal)





metagen_summary <- function(slope, slope_se, population){
  mar <- metagen(TE = slope, seTE = slope_se, studlab = population, common = F, random = T) %>%
    metarep(x = .,u = 2, report.u.max = T)
  return(data.frame(pooled_effect = mar$TE.random,
                    pooled_se = mar$seTE.random,
                    pooled_pvalue = mar$pval.random,
                    r_value = mar$r.value,
                    u_twosided = max(mar$u_L, mar$u_R),
                    u_max_L = mar$u_L, 
                    u_max_R = mar$u_R,
                    inconsistency = (mar$u_R > 0)&(mar$u_L>0) ))
}

meta_repli_results <- all_pop_data %>% 
  group_by(phenotype_id, variant_id) %>%
  reframe( meta_analysis =  metagen_summary(slope, slope_se, population)) %>% 
  unnest(cols = c(meta_analysis))

View(meta_repli_results)


# associations with replicated effects (not inclusing inconsistent)
meta_repli_results %>% 
  filter(!inconsistency) %>% 
  filter((u_twosided>2)) %>% 
  arrange(desc(u_twosided))

# assossiations with established inconsistency at 95%
meta_repli_results %>% filter(inconsistency)
