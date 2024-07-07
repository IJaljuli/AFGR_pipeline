#!/usr/bin/env bash

fd "all.(GWD|LWK|MSL|YRI|ESN).30x.ID.vcf.gz$" vcfs | sed -e 's/vcfs\///' | sed -e 's/\.vcf\.gz//' | xargs -I % plink --make-bed --maf 0.05 --vcf vcfs/%.vcf.gz --out ./plink/%

fd "MKK.filtered.biallelicAndRare.autosomal.chr.id.vcf.gz$" vcfs | sed -e 's/vcfs\///' | sed -e 's/\.vcf\.gz//' | xargs -I % plink --make-bed --maf 0.05 --vcf vcfs/%.vcf.gz --out ./plink/%
