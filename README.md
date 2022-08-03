## RASCL: RAPID ASSESSMENT OF SELECTION IN CLADES THROUGH MOLECULAR SEQUENCE ANALYSIS

### Overview
This application is designed to use molecular sequence data from genotypically distinct viral lineages to identify distinguishing features and evolution within lineages. Using whole genome sequences, a "query" set of sequences will be compared to against a globally diverse set of "background" sequences. The background data set contains globally circulating viral sequences, and the query data set is the set of sequences you want to compare. The application uses a number of open-source tools, as well as selection analysis tools from [HyPhy](hyphy.org), and assembles the results from the analysis into JSON files which can then be visualized with our full feature [Observable notebook]. We provide a list of selected results for several SARS-CoV-2 clades at (https://observablehq.com/@aglucaci/rascl)

### Installation and dependencies

This application is currently designed to run in an HPC environment.

There is an assumption that the freely available [Anaconda](https://anaconda.org/) software is installed on your machine.

#### To install -- Steps necessary to complete before running
1. `git clone https://github.com/veg/RASCL.git`
2. `cd RASCL`
3. `conda env create -f environment.yml`.  This will create a conda environment called (RASCL) with the necessary dependencies.
4. At this point, run `conda activate RASCL` and your environment will be ready to go.

#### Minimum Configuration -- Steps necessary to complete before running

The user input data (which consists of the clade of interest downloaded as a FASTA file of viral whole genome's) should be stored in the `./data}` subdirectory. We provide demo data for an test-run using the sequences in `data/Example`. These correspond to a set of "background" pre-Alpha variant set of sequences `data/Example/Background-preAlpha.fasta` and a "query" set of sequences corresponding to Alpha variant sequences `data/Example/Query-Alpha.fasta`

The LABEL variable corresponds to your viral clade of interest (e.g. "B.1.1.7") and will be used for annotation. 

1. Place your viral clade of interest fasta file within the "data" directory.
2. In the `config.json` change the following:
       The `Label` variable corresponds to a tag for your clade of interest (e.g. "B.1.1.7").
       The `Background_WholeGenomeSeqs` variable to correspond to your query whole genome sequences (e.g. "Example1/Background-preAlpha.fasta")
       The `Query_WholeGenomeSeqs` variable to correspond to your query whole genome sequences (e.g. "Example1/Query-Alpha.fasta")
  
3. `cluster.json` can be modified for your HPC environment. If you want to use more cores, adjust the values in this file. This can be used to distribute jobs to run across the cluster and to specify a queue.

The results of running this application will be placed in the `./results/{LABEL}` subdirectory. This will contain a new folder with the name of of your clade i.e. the `"LABEL"` variable from the `config.json`. We will store all intermediate files and JSON results in this subdirectory. However, they are not tracked by this GitHub repository.

At the conclusion of the run, the selection output files (BGM, MEME, FEL, SLAC, BUSTED[S], PRIME, FADE, RELAX, and Contrast-FEL, etc) will be aggregated into two JSON files (Summary.json and Annotation.json) for an Observable notebook to ingest. At this point, the user can use our visualizations to investigate the nature and extent of selective forces acting on viral genes within the clade of interest.

#### Advanced Configuration

The `config.json` file also contains a number of advanced features corresponding to the parameters we use for downsampling viral gene sequences. Specificially, from the total number of query or background sequences we aim to downsample to the `max_background` and `max_query` sequences. These can be modified by the user in order to capature an additional number of sequences from their input. We also have two additional values for `threshold_query` and `threshold_background` which correspond to the initial genetic distance threshold we apply during downsampling. 

### Running the analysis

We provide an example HPC bash script to run the analysis in `run_Silverback.sh` which is designed to run on the Temple University computing cluster. This file can be modified to run in your own computing environment. In the `cluster.json` specify the name of the queue on your system, along with the computing resources to be used.

Note, that in some cases not all of the pipeline steps will complete (e.g. insufficient sequences to run analyses on all gene segments). In this case please run, from the top RASCL directory, (with the value of `LABEL` from `config.json`, and `BASEDIR` corresponding to the main directory of your analysis.

```
bash scripts/process_json.sh {BASEDIR} {LABEL}
```

The results of the analysis will be placed into the `results/LABEL` directory as `LABEL_summary.json` and `LABEL_annotation.json`.

### Visualization

At the completion of the pipeline, the JSON outputs (Summary.json and Annotation.json) will be generated. These can be ingested into our full feature [Observable Notebook](https://observablehq.com/@aglucaci/rascl_latest). We suggest that users make a free account on ObservableHQ and fork this notebook, which allows the user to point the notebook to their data.

The version of the notebook at https://observablehq.com/@spond/sars-cov-2-clades allows one to upload summary and annotation JSON files.

#### Exploring results with our interactive notebook

We provide visualizations, an alignment viewer, site-level phylogenetic trees, and summary results and full tables in our interactive notebooks. You can explore all of the results for a particular gene through the dropdown box. Or review full results for a particular site of interest (see below).

![](https://i.imgur.com/Da3p3x0.gif)

### Galaxy workflow

We also provide an alternative way to use RASCL within the [Galaxy](https://galaxy.hyphy.org) ecosystem. User accounts are free to sign up.

[https://galaxy.hyphy.org/u/hyphy/w/rapid-assessment-of-selection-on-clades-and-lineages](https://galaxy.hyphy.org/u/hyphy/w/rapid-assessment-of-selection-on-clades-and-lineages).

 

