## RASCL: RAPID ASSESSMENT OF SARS-COV-2 CLADES THROUGH MOLECULAR SEQUENCE ANALYSIS

### Overview
This application is designed to use molecular sequence data from genotypically distinct viral lineages of SARS-CoV-2 to identify distinguishing features and evolution within lineages. Using whole genomes, a "query" set of sequences will be compared to against a globally diverse set of "background" sequences. The background data set contains globally circulating SARS-CoV-2 sequences, and the query data set is the set of sequences you want to compare. The application uses a number of open-source tools, as well as selection analysis tools from [HyPhy](hyphy.org), and assembles the results from the analysis into JSON files which can then be visualized with our full feature [Observable notebook](https://observablehq.com/@aglucaci/sars-cov-2-clades)

### Installation and dependencies

This application is currently designed to run in an HPC environment.

There is an assumption that the freely available [Anaconda](https://anaconda.org/) software is installed on your machine.

#### To install -- Steps necessary to complete before running
1. `git clone https://github.com/veg/SARS-CoV-2_Clades.git`
2. `conda env create -f environment.yml`.  This will create a virtual environment called (RASCL) with the necessary dependencies.
3. At this point, run `conda activate RASCL` and your environment will be ready to go.

#### Configuration -- Steps necessary to complete before running

The user input data (which consists of the clade of interest downloaded from GISAID as a whole genome fasta) should be stored in the `./data/{LABEL}` subdirectory, where the LABEL variable is a folder which corresponds to your clade of interest (i.e. B.1.1.7). 

1. Use `mkdir data/{LABEL}` to create the directory to house your data, then place your clade of interest fasta file within this directory.
2. In the `snakemake_config.json` change the `LABEL` variable to point to the label of your data set folder (Important: they need to match). So we will use 'B.1.1.7' in both cases for this example' (i.e. `"LABEL":"B.1.1.7"`)
3. Additionally, in the `snakemake_config.json` change the `GISAID_WG` variable to your file name (i.e. `"GISAID_WG":"gisaid_hcov-19_2021_05_11_19.fasta"`). This fasta file should already be placed within the `./data/{LABEL}` folder (i.e. `/data/B.1.1.7/gisaid_hcov-19_2021_05_11_19.fasta`)
4. `cluster.json` can be modified for your HPC environment. If you want to use more cores, adjust the values in this file. This can be used to distribute jobs to run across the cluster and to specify a queue.

The results of running this application will be placed in the `./results/{LABEL}` subdirectory. This will contain a new folder with the name of of your clade i.e. the `"LABEL"` variable from the `snakemake_config.json`. We will store all intermediate files and JSON results in this subdirectory. However, they are not tracked by this GitHub repository.

At the conclusion of the run, the selection output files (BGM, MEME, FEL, SLAC, BUSTED[S], PRIME, FADE, aBSREL, RELAX, and Contrast-FEL) will be aggregated into two JSON files (Summary.json and Annotation.json) for an [Observable notebook](https://observablehq.com/@aglucaci/sars-cov-2-clades) to ingest. At this point, the user can use our visualizations to investigate the nature and extent of selective forces acting on SARS-CoV-2 genes within the clade of interest.

### Running the analysis

We provide an example HPC bash script to run the analysis in `run_Silverback.sh` which is designed to run on the Temple University computing cluster. This file can be modified to run in your own computing environment. In the `cluster.json` specify the name of the queue on your system, along with the computing resources to be used.

### Visualization

At the completion of the pipeline, the JSON outputs (Summary.json and Annotation.json) will be generated. These can be ingested into our full feature [Observable Notebook](https://observablehq.com/@aglucaci/sars-cov-2-clades). 

#### Exploring results with our interactive notebook

We provide visualizations, an alignment viewer, site-level phylogenetic trees, and summary results and full tables in our interactive notebooks. You can explore all of the results for a particular gene through the dropdown box (see below). Or review full results for a particular site of interest (see below).
![](https://i.imgur.com/7UrADgu.gif)

