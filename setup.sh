#!/bin/bash

printf "\nSETTING UP ENV: rascl\n"

#python3 -m venv sc2-clade
#source ./sc2-clade/bin/activate

conda create -n rascl biopython=1.70 hyphy=2.5.31 python-bioext=0.19.7 raxml-ng=1.02 snakemake=5.3 -y 


#pip3 install numpy 
#pip3 install Cython 
#pip3 install bio 
#pip3 install bioext

printf "\nSUCCESSFULL INSTALL\n"

conda activate rascl

#printf "\nSUCCESSFULL INSTALL\n\nRUN --> source ./sc2-clade/bin/activate\n"
