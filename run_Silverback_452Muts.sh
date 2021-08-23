#!/bin/bash

set -euo pipefail

#printf "Running snakemake...\n"

#snakemake --forceall --dag | dot -Tpdf > dag.pdf

mkdir -p logs

#snakemake \
#      -s Snakefile \
#      --cluster-config cluster.json \
#      --cluster "qsub -V -l nodes={cluster.nodes}:ppn={cluster.ppn} -q {cluster.name} -l walltime=72:00:00 -e logs -o logs" \
#      --jobs 25 all \
#      --rerun-incomplete \
#      --keep-going \
#      --reason \
#      --latency-wait 60 \
#      --use-conda


DATADIR=/data/shares/veg/SARS-CoV-2/452Mutations/export

# Prioritize VOIs/VOCs
# https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/

#VOIs=(B.1.525.452.fas B.1.526.452.fas B.1.617.1.452.fas C.37.452.fas)
VOIs=(B.1.526.452.fas B.1.617.1.452.fas C.37.452.fas)

VOCs=(B.1.1.7.452.fas B.1.351.452.fas P.1.452.fas P.1.1.452.fas B.1.617.2.452.fas AY.1.452.fas AY.3.452.fas)
#VOCs=(AY.3.452.fas)


# Some VOCs have multiple lineage designations e.g. 617.2, AY.1, AY.3 or P.1. and P.1.1

#for clade in "${VOCs[@]}"
for clade in "${VOIs[@]}"
do
   echo "# Processing: "$clade
   wg_input="$DATADIR"/$clade
   echo "# Input whole genome file: "$wg_input

   # Pass in the clade label
   label=${clade%".fas"}
   echo "# Label (clade): "$label

   # Write out alternative JSON (snakemake_config.json)
   echo "# Writing json output for label"
   jq -n --arg label "$label" --arg clade "$clade" --arg DATADIR "$DATADIR" '{"LABEL": $label, "GISAID_WG": $clade, "DATADIR": $DATADIR}' > snakemake_config_452Muts.json
   
   # Command
   snakemake \
       -s Snakefile_452Muts \
       --cluster-config cluster.json \
       --cluster "qsub -V -l nodes={cluster.nodes}:ppn={cluster.ppn} -q {cluster.name} -l walltime=72:00:00 -e logs -o logs" \
       --jobs 30 all \
       --rerun-incomplete \
       --keep-going \
       --reason \
       --latency-wait 60 \
       --use-conda

   #if [ "$?" -eq "0" ]
   #then
   #    echo "ok"
   #else
   #    echo "Fail"
   #    continue
   #fi

      
   #continue
   echo "# Continuing to the next clade"
done




#exit 0 

