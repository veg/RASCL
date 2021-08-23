#!/bin/bash

FASTA=$1
OUTPUT=$2

echo "Cleaning: "$FASTA

sed 's, ,_,g' -i $FASTA

awk '{ if ($0 ~ "^>") {b=gensub(/>(.+)\|(EPI_ISL_|epi_isl_)([0-9]+)\|(.+)/, ">epi_isl_\\3/\\1","g"); print b;} else print;}' $FASTA > $OUTPUT

# For 452 MUTS
#epi_isl_2870330,AY.3

sed "s/\,/\_/g" -i $OUTPUT
sed "s/\./\_/g" -i $OUTPUT

exit 0
