#!/bin/bash
#@Usage: bash process_json.sh

# Description: Really runs on the master node for now. Creates summary and annotation jsons

## Declares
#BASEDIR="/home/aglucaci/SARS-CoV-2/clades"

BASEDIR=/home/aglucaci/SARS-CoV-2/clades/SNAKEMAKE/results

#TAG="B-1-617-1"
TAG="B-1-617-2"

DATA_DIR="$BASEDIR"/"$TAG"

# Static settings

REF_TAG="REFERENCE"
ANNOTATION_JSON="$TAG"_annotation.json
SUMMARY_JSON="$TAG"_summary.json

for file in "$DATA_DIR"/*.combined.fas; do
   #echo $file
   echo python3 generate-report.py -f $file -A $ANNOTATION_JSON -S $SUMMARY_JSON -r "$REF_TAG"
   python3 generate-report.py -f $file -A $ANNOTATION_JSON -S $SUMMARY_JSON -r "$REF_TAG"
done

exit 0












# end of file
