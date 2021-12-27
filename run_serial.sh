#!/bin/bash

CORES="${1:-1}"

set -euo pipefail

printf "Running snakemake on ${CORES} cores\n"

# Uncomment this command to create the pipeline DAG
#snakemake --forceall --dag | dot -Tpdf > dag.pdf

mkdir -p logs

snakemake \
      -s Snakefile.serial \
      --cores $CORES \
      --jobs 23 all \
      --rerun-incomplete \
      --keep-going \
      --reason \
      --latency-wait 60 \
      -p 