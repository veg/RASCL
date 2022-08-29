#!/bin/bash

clear

echo "Version: v0.1 --- "
echo "2022, RAPID ASSESSMENT OF SELECTION IN CLADES (RASCL)."
echo ""

# Set up the pipeline failure expectations.
set -euo pipefail

echo "Initialized --- "

# Uncomment if you want to generate an analysis DAG file.
#snakemake --forceall --dag | dot -Tpdf > RASCL_DAG.pdf

echo "Creating 'logs' directory"
mkdir -p logs

echo "Executing LOCAL Snakemake command"

# Execute the Snakemake command
snakemake \
      -s Snakefile \
      --cluster-config cluster.json \
      --jobs 1 all \
      --rerun-incomplete \
      --keep-going \
      --reason \
      --latency-wait 60 

# End Snakemake command

exit 0

# End of file.
