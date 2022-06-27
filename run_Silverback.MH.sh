#!/bin/bash

set -euo pipefail

printf "Running snakemake...\n"

# Uncomment this command to create the pipeline DAG
#snakemake --forceall --dag | dot -Tpdf > dag.pdf

mkdir -p logs

snakemake \
      -s Snakefile.MH \
      --cluster-config cluster.json \
      --cluster "qsub -V -l nodes={cluster.nodes}:ppn={cluster.ppn} -q {cluster.name} -l walltime=120:00:00 -e logs -o logs" \
      --jobs 23 all \
      --rerun-incomplete \
      --keep-going \
      --reason \
      --latency-wait 60 


# --cluster "qsub -V -N RASCL -l nodes={cluster.nodes}:ppn={cluster.ppn} -q {cluster.name} -l walltime=120:00:00 -e logs -o logs" \

exit 0
