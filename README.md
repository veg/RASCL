## RASCL: RAPID ASSESSMENT OF SARS-COV-2 CLADES THROUGH MOLECULAR SEQUENCE ANALYSIS

### Overview
This application is designed to use molecular sequence data from genotypically distinct viral lineages of SARS-CoV-2 to identify distinguishing features and evolution within lineages.
Using whole genomes, a "query" set of sequences will be compared to against a globally diverse set of "background" sequences. The background data set contains globally circulating SARS-CoV-2 sequences, and the query data set is the set of sequences you want to compare.
The application uses a number of open-source tools, as well as selection analysis tools from [HyPhy](hyphy.org), and assembles the results from the analysis into JSON files which can then be visualized with our full feature [Observable notebook](https://observablehq.com/@aglucaci/sars-cov-2-clades)

### Installation
*There is an assumption that [ANACONDA](https://anaconda.org/) is installed on your machine*

1. `git clone` this repo 
2. in the base directory, run `bash ./setup.sh`. This will create a virtual env (rascl) with the necessary dependencies.

At this point, run `conda activate rascl` and your env will be ready to go. Next step is putting your data in the correct place.

### Data Input -- Steps necessary to complete before running
1. in the `snakemake_config.json` change the `LABEL` to the label of your data set (i.e. `"LABEL":"B.1.1.7"`)
2. (optional) add your own reference file, in the `snakemake_config.json` change the `GISAID_WG` to your file name (i.e. `"GISAID_WG":"your_file.fasta"`) 


### Directory Structure 


### How to use the pipeline
1. Clone the repository
2. Download [GISAID](https://www.gisaid.org/) 
3. use `mkdir` to put the GISAID data into the `/data/{LABEL}` folder, where `LABEL` is the clade you will analyse: for example, "B-1-617" or "B-1-1-7" (`mkdir data/B-1-617`)
4. Configure the `env` to set up the environment -- THIS ASSUMES YOU HAVE -- python3 `bash ./setup.sh` then `source ./sc2-clade/bin/activate` 
5. edit the `LABEL` variable in the `snakemake_config.json`
6. edit the `GISAID_WG` variable in the `snakemake_config.json`

`conda env create -f environment.yml`

[BIOCONDA](https://bioconda.github.io/user/install.html)

5. execute the snakefile (see below if on server) 

#### If you are running on the server then run:
```
snakemake -s Snakefile --cluster "qsub -V -l nodes=1:ppn=8 -q epyc2 -l walltime=999:00:00" --jobs 50 all --rerun-incomplete --keep-going -use-conda
```


