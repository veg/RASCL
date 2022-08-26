## RASCL: RAPID ASSESSMENT OF SELECTION IN CLADES THROUGH MOLECULAR SEQUENCE ANALYSIS

### Overview
This application is designed to use molecular sequence data from genotypically distinct viral lineages to identify distinguishing features and evolution within lineages. Using whole genome sequences, a "query" set of sequences will be compared to against a globally diverse set of "background" sequences. The background data set contains globally circulating viral sequences, and the query data set is the set of sequences you want to compare. The application uses a number of open-source tools, as well as selection analysis tools from [HyPhy](hyphy.org), and assembles the results from the analysis into JSON files which can then be visualized with our full feature [Observable notebook]. We provide a list of selected results for several SARS-CoV-2 clades at (https://observablehq.com/@aglucaci/rascl)

### Installation and dependencies (Conda-based)
#### Environment dependencies
This application is currently designed to run in an HPC environment.

There is an assumption that the freely available [Anaconda](https://anaconda.org/) software is installed on your machine.

#### To install (Conda-based) -- Steps necessary to complete before running
1. `git clone https://github.com/veg/RASCL.git`
2. `cd RASCL`
3. `conda env create -f environment.yml`.  This will create a conda environment called (RASCL) with the necessary dependencies.
4. At this point, run `conda activate RASCL` and your environment will be ready to go.

#### Configuration settings -- Steps necessary to complete before running

The user input data (which consists of the clade of interest downloaded as a FASTA file of viral whole genome's) should be stored in the `./data` subdirectory. We provide demo data for a test-run using the sequences in `data/Example1`. These correspond to a set of "background" pre-Alpha variant set of sequences `data/Example1/Background-preAlpha.fasta` and a "query" set of sequences corresponding to Alpha variant sequences `data/Example1/Query-Alpha.fasta`

The `Label` variable within the `config.json` corresponds to the label of your viral clade of interest (e.g. "B.1.1.7") and will be used for annotation downstream. 

From your working directory:

1. `mkdir data/{name}`
2. Place your viral clade of interest fasta file within the `data/{name}` subdirectory.
- The `{name}` you pick does not matter. Recall, in our test-run our `Label` is `B.1.1.7`, and we put our input files in the `./data/Example1` directory.
3. In the `config.json` change the following:
- The `Label` variable to correspond to a tag for your clade of interest (e.g. "B.1.1.7"). Make sure to include the `"` around your label. 
- The `Background_WholeGenomeSeqs` variable to correspond to your query whole genome sequences (e.g. "Example1/Background-preAlpha.fasta"). Include the relative path as if you were within the `./data/` directory.
- The `Query_WholeGenomeSeqs` variable to correspond to your query whole genome sequences (e.g. "Example1/Query-Alpha.fasta"). Include the relative path as if you were within the `./data/` directory.
- The `cleaner` variable corresponds to bash script containing commands necessary to parse the headers of the fasta. Within the `./scripts` directory, there are two files, one with `linux` (defualt), and one with `MBP` in the name. change this file according to your OS.


4. The `cluster.json` file can be modified for your computing environment. If you want to use more cores, adjust the values in this file. This can be used to distribute jobs to run across the cluster and to specify a queue. The `cluster` variable refers to the workload manager. The `nodes` variable is a request for resource allocation from the server, in this case it refers to the number of nodes. The `ppn` variable is a request for resource allocation from the server, in this case it refers to the number of processors per node. The `name` variable is a specification to submit the jobs for the RASCL application to a specific queue. These have different names and priorities, please refer to your local system administrator for more information. We have added an additional variable `walltime` which is a request for a certain period of time for resource allocation from the server.

*Notes for running on your local machine.* If you have the computational ability to parallelize jobs across more than one processor, then you can adjust the `ppn` variable (processors per node). Starting with 1 `ppn` is a good place to start, as some computations can take up quite a bit of compute power. The remaining variables within this file are not relevant to a local run, so they can be left alone (see **Advanced Configurations** for more information).

5. Within the `run_LOCAL.sh` file (within the RASCL directory), change the `--jobs` parameter, to reflect the number of jobs you want to parallelize at one time. Again, much like the `ppn` variable, starting with 1 is a good idea.

At this point, your configuration settings are set.

#### Running the analysis.

When in the `RASCL` working directory:
1. `bash ./run_HPC.sh` for running on a high performance computing server
2. `bash ./run_LOCAL.sh` for running on your local machine

The results of running this application will be placed in the `results/{Label}` subdirectory. This will contain a new folder with the name of of your clade i.e. the `"Label"` variable from the `config.json`. We will store all intermediate files and JSON results in this subdirectory. However, they are not tracked by this GitHub repository.

At the conclusion of the run, the selection output files (BGM, MEME, FEL, SLAC, BUSTED[S], PRIME, FADE, RELAX, and Contrast-FEL, etc) will be aggregated into two JSON files (Summary.json and Annotation.json) for an Observable notebook to ingest. At this point, the user can use our visualizations to investigate the nature and extent of selective forces acting on viral genes within the clade of interest.

### Conda-independent Installation (utilizes a python `venv`, unless you want to install system-wide, then omit step 1)
#### Environment dependencies
1. `git clone https://github.com/veg/RASCL.git`
2. `cd RASCL`
2. `python3.8 -m venv tester` Note: (requires python3.8 to run, even if you dont want to use a `venv` you need python 3.8)
3. `source tester/bin/activate
4. `pip install biopython==1.77`
5. `pip install snakemake==7.9.0`
6. `pip install Cython==0.29.32`
7. `pip install bioext=0.19.7`
8. `git clone --recursive https://github.com/amkozlov/raxml-ng` (follow developer install directions)
  - requies version `1.1.0`
9. `git clone https://github.com/veg/tn93.git` (follow developer install directions)
  - requires version `1.0.9`
10. `git clone https://github.com/veg/hyphy.git` (follow developer install directions)
  - requies version `2.5.41`
11. `brew install gnu-sed`

Note: Change the `raxml_ng` variable within the `config.json`, which corresponds to the full path to the `raxml-ng` executable that was installed, it should look something like: `/usr/path/to/raxml-ng/bin/raxml-ng` for local runs, as well as the `cleaner` variable to either the file containing `linux` (defualt), or to the one with `MBP`, depending on your OS.

### Running the analysis

We provide an example HPC bash script to run the analysis in `run_Silverback.sh` which is designed to run on the Temple University computing cluster. This file can be modified to run in your own computing environment. In the `cluster.json` specify the name of the queue on your system, along with the computing resources to be used.

Note, that in some cases not all of the pipeline steps will complete (e.g. insufficient sequences to run analyses on all gene segments). In this case please run, from the top RASCL directory, (with the value of `Label` from `config.json`, and `WD` corresponding to the working directory.

```
bash scripts/process_json.sh {WD} {Label}
```

The results of the analysis will be placed into the `results/Label` directory as `{Label}_summary.json` and `{Label}_annotation.json`.

See `Visualization` section for next steps on how to use our interactive notebook.

### Visualization

At the completion of the pipeline, the JSON outputs (Summary.json and Annotation.json) will be generated. These can be ingested into our full feature [Observable Notebook](https://observablehq.com/@aglucaci/rascl_latest). We suggest that users make a free account on ObservableHQ and fork this notebook, which allows the user to point the notebook to their data.

The version of the notebook at https://observablehq.com/@spond/sars-cov-2-clades allows one to upload summary and annotation JSON files.

If you get an error about `tree_tags`, go to the `species` variable line and change it to the gene you analyzed.

#### Exploring results with our interactive notebook

We provide visualizations, an alignment viewer, site-level phylogenetic trees, and summary results and full tables in our interactive notebooks. You can explore all of the results for a particular gene through the dropdown box. Or review full results for a particular site of interest (see below).

![](https://i.imgur.com/Da3p3x0.gif)

### Galaxy workflow

We also provide an alternative way to use RASCL within the [Galaxy](https://galaxy.hyphy.org) ecosystem. User accounts are free to sign up.

[https://galaxy.hyphy.org/u/hyphy/w/rapid-assessment-of-selection-on-clades-and-lineages](https://galaxy.hyphy.org/u/hyphy/w/rapid-assessment-of-selection-on-clades-and-lineages).

#### Advanced Configuration

The `config.json` file also contains a number of advanced features corresponding to the parameters we use for downsampling viral gene sequences. Specificially, from the total number of query or background sequences we aim to downsample to the `max_background` and `max_query` sequences. These can be modified by the user in order to capature an additional number of sequences from their input. We also have two additional values for `threshold_query` and `threshold_background` which correspond to the initial genetic distance threshold we apply during downsampling. 

#### Advanced Configuration for HPC environment and downsampling

If you want to use more cores, adjust the values in the `cluster.json` file. This can be used to distribute jobs to run across the cluster and to specify a queue. The `cluster` variable refers to the workload manager. The `nodes` variable is a request for resource allocation from the server, in this case it refers to the number of nodes. The `ppn` variable is a request for resource allocation from the server, in this case it refers to the number of processors per node. The `name` variable is a specification to submit the jobs for the RASCL application to a specific queue. These have different names and priorities, please refer to your local system administrator for more information. We have added an additional variable `walltime` which is a request for a certain period of time for resource allocation from the server. We provide an example HPC bash script to run the analysis in `run_HPC.sh` which is designed to run on the Temple University computing cluster. This file can be modified to run in your own computing environment. In the `cluster.json` specify the name of the queue on your system, along with the computing resources to be used.

#### Testing on a singular gene 
If you want to run the script on a singular gene instead of the entire SARS-CoV-2 genome, you can go into the Snakemake file (`Snakefile`), comment out line 58, 59, 60, and 61. Uncomment out line 64 and input the gene that you want to run (from the `gene` list).
 

