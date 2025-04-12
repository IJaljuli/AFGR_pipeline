#!/usr/bin/env bash
#ls
## Allow aliases to work in a non-interactive script - only needed for bash
#SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
#cd SCRIPT_DIR
#cd .. 
#cd .. 
#SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
#SCRIPT_DIR

# shopt -s expand_aliases

# From here, you can just run tensorqtl
# Need to make plink files from VCF
# For loop runs tensorQTL in cis mode (result: eQTLs)  on each population separately
for popName in ESN GWD LWK MSL YRI MKK; do     
    
    echo $popName

    ba_path="/oak/stanford/groups/smontgom/jaljuli/afgr_eqtls/pipe_repo/data/bed_bim_fam/collapsed.hwepass01.withgnomad.var.merged.processedMKK.split1kG30x.${popName}.maf05.geno01"


    bed_path="/oak/stanford/groups/smontgom/jaljuli/afgr_eqtls/pipe_repo/data/bed_bim_fam/${popName}_phenotype.bed"

    covariates_path="/oak/stanford/groups/smontgom/jaljuli/afgr_eqtls/pipe_repo/data/covariates/${popName}_covariates.txt"

    export PATH="${HOME}/.pixi/bin:/oak/stanford/groups/smontgom/jaljuli/pixi/bin:${PATH}"
    tensorqtl \
        "${ba_path}" \
        "${bed_path}" \
        ${popName}_ \
        --covariates "${covariates_path}" \
        --mode cis_nominal \
        --maf_threshold 0.05 \
        --pval_threshold 1 \
        --seed 1 \
        --output_dir /oak/stanford/groups/smontgom/jaljuli/afgr_eqtls/pipe_repo/output/cis_parquet

    bed_splicing_path="/oak/stanford/groups/smontgom/jaljuli/afgr_eqtls/pipe_repo/data/bed_bim_fam/${popName}_phenotype_splicing.bed"
    covariates_splicing_path="/oak/stanford/groups/smontgom/jaljuli/afgr_eqtls/pipe_repo/data/covariates/${popName}_covariates_splicing.txt"
    tensorqtl \
        "${ba_path}" \
        "${bed_path}" \
        ${popName}_ \
        --covariates "${covariates_path}" \
        --mode cis_nominal \
        --maf_threshold 0.05 \
        --pval_threshold 1 \
        --seed 1 \
        --windows 100000 \
        --output_dir /oak/stanford/groups/smontgom/jaljuli/afgr_eqtls/pipe_repo/output/cis_splicing_parquet
done

