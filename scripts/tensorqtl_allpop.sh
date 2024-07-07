#!/usr/bin/env bash
#ls
## Allow aliases to work in a non-interactive script - only needed for bash
#SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
#cd SCRIPT_DIR
#cd .. 
#cd .. 
#SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
#SCRIPT_DIR

shopt -s expand_aliases

alias tensorqtl="micromamba run --prefix=/oak/stanford/groups/smontgom/jaljuli/micromamba/envs/tensorqtl tensorqtl"

# From here, you can just run tensorqtl
# Need to make plink files from VCF
echo popName
# For loop runs tensorQTL in cis mode (result: eQTLs)  on each population separately
for popName in ESN; do   # GWD LWK MSL YRI ; do # MKK    
    
    ba_path="/oak/stanford/groups/smontgom/jaljuli/afgr_eqtls/pipe_repo/data/bed_bim_fam/all.${popName}.30x.ID"

    bed_path="/oak/stanford/groups/smontgom/jaljuli/afgr_eqtls/pipe_repo/data/bed_bim_fam/${popName}_phenotype.bed"

    covariates_path="/oak/stanford/groups/smontgom/jaljuli/afgr_eqtls/pipe_repo/data/covariates/${popName}_covariates.txt"

    #tensorqtl \
        #"${ba_path}" \
        #"${bed_path}" \
        #${popName}_\
        #--covariates "${covariates_path}" \
        #--mode cis \
        #--maf_threshold 0.05 \
        #--seed 1 \
        #--output_dir /oak/stanford/groups/smontgom/jaljuli/afgr_eqtls/pipe_repo/output/cis

    tensorqtl \
        "${ba_path}" \
        "${bed_path}" \
        ${popName}_\
        --covariates "${covariates_path}" \
        --mode cis_susie \
        --maf_threshold 0.05 \
        --seed 1 \
        --cis_output /oak/stanford/groups/smontgom/jaljuli/afgr_eqtls/pipe_repo/output/cis/${popName}_.cis_qtl.txt.gz \
        --output_dir /oak/stanford/groups/smontgom/jaljuli/afgr_eqtls/pipe_repo/output/cis_susie

done

