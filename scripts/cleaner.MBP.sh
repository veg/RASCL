#!/bin/bash

FASTA=$1
OUTPUT=$2

mkdir -p results
#mkdir -p results/Example1

echo "# Input FASTA: "$FASTA
echo "# Output FASTA: "$OUTPUT

# Creating a copy of the data in the results folder
echo ""
echo "Creating a copy of the data, placed in the results folder"
echo "cp $FASTA ${OUTPUT}.bk"
cp $FASTA ${OUTPUT}.bk


echo ""
echo "Cleaning our FASTA header information in: "${OUTPUT}.bk
echo "Replacing illegal characters with underscores"

## Backup here will be in the data folder
## inplace modifications
sed -i '' 's/ /_/g'  ${OUTPUT}.bk
sed -i '' 's/:/_/g'  ${OUTPUT}.bk
sed -i '' 's/|/_/g'  ${OUTPUT}.bk
sed -i '' 's/-/_/g'  ${OUTPUT}.bk
sed -i '' 's/\//_/g' ${OUTPUT}.bk
sed -i '' 's/\./_/g' ${OUTPUT}.bk

## Run GAWK
## This is for GISAID data.
#echo ""
#echo "Running awk to clean up header information"
##gawk '{ if ($0 ~ "^>") {b=gensub(/>(.+)\|(EPI_ISL_|epi_isl_)([0-9]+)\|(.+)/, ">epi_isl_\\3/\\1","g"); print b;} else print;}' ${OUTPUT}.bk > $OUTPUT
awk '{ if ($0 ~ "^>") {b=gensub(/>(.+)\|(EPI_ISL_|epi_isl_)([0-9]+)\|(.+)/, ">epi_isl_\\3/\\1","g"); print b;} else print;}' ${OUTPUT}.bk 

mv ${OUTPUT}.bk ${OUTPUT}

#echo ""
#echo "Removing backup .bk files"
## Delete backup file
#rm ${OUTPUT}.bk


exit 0
