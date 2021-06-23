#!/bin/bash

printf "\nSETTING UP ENV: rascl\n"

conda create -n rascl biopython=1.70 hyphy=2.5.31 python-bioext=0.19.7 raxml-ng=1.02 snakemake=5.3 -y 

printf "\nSUCCESSFULL INSTALL\n"

conda activate rascl

exit 0
