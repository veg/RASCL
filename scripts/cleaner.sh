#!/bin/bash

FASTA=$1
OUTPUT=$2

echo "Cleaning: "$FASTA


sed 's, ,_,g' -i $FASTA

awk '{ if ($0 ~ "^>") {b=gensub(/>(.+)\|(EPI_ISL_|epi_isl_)([0-9]+)\|(.+)/, ">epi_isl_\\3/\\1","g"); print b;} else print;}' $FASTA > $OUTPUT





exit 0
