#!/usr/bin/env bash

mamba run -n gtex-pipelines collapse_annotation.py ./oak/stanford/groups/smontgom/dnachun/data/AFGR/rnaseqc/Homo_sapiens.GRCh38.110.genes.gtf ./Homo_sapiens.GRCh38.110.collapse.gtf
