# SARS-CoV-2 Clades analysis pipeline

### Overview.

The pipeline uses several open-source tools to prepare SARS-CoV-2 full length genome data for clade selection analysis using [HyPhy](hyphy.org), assemble the results of these analysis using a Python3 script into a monolithic JSON file which can then be visualized with an [Observable notebook](https://observablehq.com/@aglucaci/sars-cov-2-clades). The results are then visualized to pick out patterns at each site, such as positive and negative selection, the history of selection at that site (when that site started becoming selected for), haplotype information, global distribution at that site, etc.  

### How to use the pipeline
1. Clone the repository
2. Download GISAID data into the data/ folder
3. Make sure that you name the folder after the clade you will analyse: for example, "B-1-617" or "B-1-1-7"
4. Configure and execute the snakefile 


