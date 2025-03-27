#!/usr/bin/env bash 

#### Collapse genes (only if needed, file is provided!)

# bash /oak/stanford/groups/jaljuli/afgr_eqtls/pipe_repo/scripts/collapse_genes.sh

#### Generate PLINK 1 files  (only if needed, files are provided!)

# bash /oak/stanford/groups/smontgom/jaljuli/afgr_eqtls/pipe_repo/scripts/generate_bam.sh



#### Run analyses in R to silter/normalize phenotype data, create covariates and phenotype .bed files 
## Note: The scripts assume that pruned_variantes.rds and .gds files are provided in designated directories. If need to regenrate, open ./sript/ancestry_pcs.R file and set 'get_gsd' and 'do_ldpruning' to TRUE 



Rscript  /oak/stanford/groups/smontgom/jaljuli/afgr_eqtls/pipe_repo/scripts/source_Rscripts.R

#### Submit an interactive job to run TensorQTL in cis and susie modes 

tmux

srun -c 1 --mem 32G -t 10:00:00 --gres=gpu:1 -A smontgom -p nih_s10 --pty zsh

bash  /oak/stanford/groups/smontgom/jaljuli/afgr_eqtls/pipe_repo/scripts/tensorqtl_allpop.sh


