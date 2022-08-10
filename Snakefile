#==============================================================================
# Snakefile for RAPID ASSESSMENT OF SELECTION IN CLADES (RASCL).
# 2022
#==============================================================================

#==============================================================================
# Imports
#==============================================================================
import os
import sys
import json
import csv
from pathlib import Path

#==============================================================================
# Declares
#==============================================================================

with open("config.json", "r") as input_sc:
  config = json.load(input_sc)
#end with

with open("cluster.json", "r") as input_c:
  cluster = json.load(input_c)
#end with

#==============================================================================
# User settings 
#==============================================================================
BASEDIR = os.getcwd()

print("We are operating out of base directory:", BASEDIR)

# Which clades are you analyzing?
LABEL = config["Label"]

# Get this information from the config.
Query_WholeGenomeSeqs = config["Query_WholeGenomeSeqs"]
Background_WholeGenomeSeqs = config["Background_WholeGenomeSeqs"]

# A different set of variables to specify full-paths
Input_Query_WholeGenomeSeqs = os.path.join(BASEDIR, "data", Query_WholeGenomeSeqs )
Input_Background_WholeGenomeSeqs = os.path.join(BASEDIR, "data", Background_WholeGenomeSeqs )

# Report to user
print("We are using the following input files")
print("Query FASTA file:", Input_Query_WholeGenomeSeqs)
print("Background FASTA file:", Input_Background_WholeGenomeSeqs)
print()

#==============================================================================
# End -- User defined settings
#==============================================================================

genes = ["leader", "nsp2", "nsp3", "nsp4", "3C", "nsp6", "nsp7", "nsp8", \ 
         "nsp9", "nsp10", "helicase", "exonuclease", "endornase", "S",   \
         "E", "M", "N", "ORF3a", "ORF6", "ORF7a", "ORF8" ,"RdRp",        \
         "methyltransferase"]

# For debugging or single gene analyses
#genes = ["S"]

# Basename of INPUT_WG (used in rule clean)
Input_Query_WholeGenomeSeqs_Basename = os.path.basename(Input_Query_WholeGenomeSeqs)
Input_Background_WholeGenomeSeqs_Basename = os.path.basename(Input_Background_WholeGenomeSeqs)

# Reference NCBI sequence dirs 
# This is hard-coded to the Wuhan reference from NCBI.
REF_SEQ_DIR = os.path.join(BASEDIR, "data", "ReferenceSeq")

# Set output directory
OUTDIR = os.path.join(BASEDIR, "results", LABEL)

# Create output directories
Path(os.path.join(BASEDIR, "results")).mkdir(parents=True, exist_ok=True)
Path(OUTDIR).mkdir(parents=True, exist_ok=True)

# Settings, these can be passed in or set in a config.json type file
PPN = cluster["__default__"]["ppn"] 

#==============================================================================
# Rule All 
#==============================================================================
rule all:
    input:
        os.path.join(OUTDIR, Input_Query_WholeGenomeSeqs_Basename + ".fa"),
        os.path.join(OUTDIR, Input_Background_WholeGenomeSeqs_Basename + ".fa"),
        expand(os.path.join(OUTDIR, "{GENE}.query.bam"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.query.msa.OG"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.query.msa.SA"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.query.compressed.fas"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.query.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.reference.bam"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.reference.msa.OG"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.reference.msa.SA"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.reference.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.reference.compressed.fas"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.combined.fas"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.AA.fas"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.combined.fas.raxml.bestTree"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.int.nwk"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.clade.nwk"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.full.nwk"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.combined.fas.BGM.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.SLAC.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.FEL.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.MEME.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.MEME-full.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.PRIME.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.FADE.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.BUSTEDS.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.RELAX.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.CFEL.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.FMM.json"), GENE=genes),
        os.path.join(OUTDIR, LABEL + "_summary.json"),
        os.path.join(OUTDIR, LABEL + "_annotation.json")
#end rule

#==============================================================================
# Rules -- Main analysis 
#==============================================================================
rule cleaner_query:
    input:
        input = Input_Query_WholeGenomeSeqs
    output:
        output = os.path.join(OUTDIR, Input_Query_WholeGenomeSeqs_Basename + ".fa")
    shell:
       "bash scripts/cleaner.sh {input.input} {output.output}"
#end rule

rule cleaner_background:
    input:
        input = Input_Background_WholeGenomeSeqs
    output:
        output = os.path.join(OUTDIR, Input_Background_WholeGenomeSeqs_Basename + ".fa")
    shell:
       "bash scripts/cleaner.sh {input.input} {output.output}"
#end rule

#==============================================================================
# PROCESS QUERY SEQUENCES
#==============================================================================
rule bealign_query:
    input:
        input_genome = rules.cleaner_query.output.output,
        input_gene_ref_seq = os.path.join(REF_SEQ_DIR, "{GENE}.fas")
    output:
        output = os.path.join(OUTDIR, "{GENE}.query.bam")
    shell:
        "bealign -r {input.input_gene_ref_seq} -m HIV_BETWEEN_F {input.input_genome} {output.output}"
#end rule

rule bam2msa_query:
    input:
        input_bam = rules.bealign_query.output.output
    output:
        output_msa = os.path.join(OUTDIR, "{GENE}.query.msa.OG")
    shell:
        "bam2msa {input.input_bam} {output.output_msa}"        
#end rule

rule strike_ambigs_query:
   input:
       input_msa = rules.bam2msa_query.output.output_msa
   output:
       output_strike_ambigs = os.path.join(OUTDIR, "{GENE}.query.msa.SA")
   shell:
      "hyphy scripts/strike-ambigs.bf --alignment {input.input_msa} --output {output.output_strike_ambigs}"
#end rule

rule tn93_cluster_query:
    params:
        THRESHOLD_QUERY = config["threshold_query"],
        MAX_QUERY = config["max_query"] 
    input:
        input_msa = rules.strike_ambigs_query.output.output_strike_ambigs
    output:
        output_fasta = os.path.join(OUTDIR, "{GENE}.query.compressed.fas"),
        output_json = os.path.join(OUTDIR, "{GENE}.query.json")
    shell:
        "python3 scripts/tn93_cluster.py --input {input.input_msa} --output_fasta {output.output_fasta} --output_json {output.output_json} --threshold {params.THRESHOLD_QUERY} --max_retain {params.MAX_QUERY}"
#end rule

#==============================================================================
# PROCESS BACKGROUND SEQUENCES
#==============================================================================
rule bealign_ref:
    input:
        input_genome = rules.cleaner_background.output.output,
        input_gene_ref_seq = os.path.join(REF_SEQ_DIR, "{GENE}.fas")
    output:
        output = os.path.join(OUTDIR, "{GENE}.reference.bam")
    shell:
        "bealign -r {input.input_gene_ref_seq} -m HIV_BETWEEN_F -K {input.input_genome} {output.output}"
#end rule

rule bam2msa_ref:
    input:
        input_bam = rules.bealign_ref.output.output
    output:
        output_msa = os.path.join(OUTDIR, "{GENE}.reference.msa.OG")
    shell:
        "bam2msa {input.input_bam} {output.output_msa}"
#end rule 

rule strike_ambigs_ref:
   input:
       input_msa = rules.bam2msa_ref.output.output_msa
   output:
       output_strike_ambigs = os.path.join(OUTDIR, "{GENE}.reference.msa.SA")
   shell:
      "hyphy scripts/strike-ambigs.bf --alignment {input.input_msa} --output {output.output_strike_ambigs}"
#end rule

rule tn93_cluster_ref:
    params:
        THRESHOLD_REF = config["threshold_background"],
        MAX_REF = config["max_background"],
    input:
        input_msa = rules.strike_ambigs_ref.output.output_strike_ambigs,
        input_gene_ref_seq = os.path.join(REF_SEQ_DIR, "{GENE}.fas")
    output:
        output_fasta = os.path.join(OUTDIR, "{GENE}.reference.compressed.fas"),
        output_json = os.path.join(OUTDIR, "{GENE}.reference.json")
    shell:
        "python3 scripts/tn93_cluster.py --input {input.input_msa} --output_fasta {output.output_fasta} --output_json {output.output_json} --threshold {params.THRESHOLD_REF} --max_retain {params.MAX_REF} --reference_seq {input.input_gene_ref_seq}"
#end rule 

#==============================================================================
# Combine them, the alignments
#==============================================================================
rule combine:
    params:
        THRESHOLD_QUERY = config["threshold_query"]
    input:
        input_compressed_fas = rules.tn93_cluster_query.output.output_fasta,
        input_msa = rules.tn93_cluster_ref.output.output_fasta,
	input_gene_ref_seq = os.path.join(REF_SEQ_DIR, "{GENE}.fas")
    output:
        output = os.path.join(OUTDIR, "{GENE}.combined.fas")
        #output_csv = os.path.join(OUTDIR, "{GENE}.combined.fas.csv")
    shell:
        "python3 scripts/combine.py --input {input.input_compressed_fas} -o {output.output} --threshold {params.THRESHOLD_QUERY} --msa {input.input_msa} --reference_seq {input.input_gene_ref_seq}"
#end rule

#==============================================================================
# Convert to protein
#==============================================================================
rule convert_to_protein:
    input:
        combined_fas = rules.combine.output.output
    output:
        proteinput_fas = os.path.join(OUTDIR, "{GENE}.AA.fas")
    shell:
        "hyphy conv Universal 'Keep Deletions' {input.combined_fas} {output.proteinput_fas}"
#end rule 

#==============================================================================
# Combined ML Tree
#==============================================================================
rule raxml:
    params:
        THREADS = PPN
    input:
        combined_fas = rules.combine.output.output
    output:
        combined_tree = os.path.join(OUTDIR, "{GENE}.combined.fas.raxml.bestTree")
    shell:
        "raxml-ng --model GTR --msa {input.combined_fas} --threads {params.THREADS} --tree pars{{3}} --force"
#end rule 

rule annotate:
    input:
       input_tree = rules.raxml.output.combined_tree,
       input_compressed_fas = rules.tn93_cluster_query.output.output_fasta
    output:
       output_int_tree = os.path.join(OUTDIR, "{GENE}.int.nwk"),
       output_clade_tree = os.path.join(OUTDIR, "{GENE}.clade.nwk"),
       output_full_tree = os.path.join(OUTDIR, "{GENE}.full.nwk")
    shell:
       "bash scripts/annotate.sh {input.input_tree} 'REFERENCE' {input.input_compressed_fas} {LABEL} {BASEDIR}"
#end rule 

#==============================================================================
# Selection analyses
#==============================================================================

rule slac:
    input:
        input_msa = rules.combine.output.output,
        input_tree = rules.annotate.output.output_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.SLAC.json")
    shell:
        "mpirun -np {PPN} HYPHYMPI SLAC --alignment {input.input_msa} --samples 0 --tree {input.input_tree} --output {output.output}"
#end rule

rule bgm:
    input:
        input_msa = rules.combine.output.output,
        input_tree = rules.annotate.output.output_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.combined.fas.BGM.json")
    shell:
        "hyphy BGM --alignment {input.input_msa} --tree {input.input_tree} --output {output.output} --branches {LABEL}"
#end rule -- bgm

rule fel:
    input:
        input_msa = rules.combine.output.output,
        input_tree = rules.annotate.output.output_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.FEL.json")
    shell:
        "mpirun -np {PPN} HYPHYMPI FEL --alignment {input.input_msa} --tree {input.input_tree} --output {output.output} --branches {LABEL}"
#end rule

rule meme:
    input:
        input_msa = rules.combine.output.output,
        input_tree = rules.annotate.output.output_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.MEME.json")
    shell:
        "mpirun -np {PPN} HYPHYMPI MEME --alignment {input.input_msa} --tree {input.input_tree} --output {output.output} --branches {LABEL}"
#end rule

rule busted:
    input:
        input_msa = rules.combine.output.output,
        input_tree_clade = rules.annotate.output.output_clade_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.BUSTEDS.json")
    shell:
        "mpirun -np {PPN} HYPHYMPI BUSTED --alignment {input.input_msa} --tree {input.input_tree_clade} --output {output.output} --branches {LABEL} --starting-points 10"
#end rule 

rule relax:
    input:
        input_msa = rules.combine.output.output,
        input_tree_clade = rules.annotate.output.output_clade_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.RELAX.json")
    shell:
        "hyphy RELAX --alignment {input.input_msa} --models Minimal --tree {input.input_tree_clade} --output {output.output} --test {LABEL} --reference Reference --starting-points 10 --srv Yes"
#end rule

rule prime:
    input:
        input_msa = rules.combine.output.output,
        input_tree = rules.annotate.output.output_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.PRIME.json")
    shell:
        "hyphy PRIME --alignment {input.input_msa} --tree {input.input_tree} --output {output.output} --branches {LABEL}"
#end rule 

rule meme_full:
    input:
        input_msa = rules.combine.output.output,
        input_tree_full = rules.annotate.output.output_full_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.MEME-full.json")
    shell:
        "mpirun -np {PPN} HYPHYMPI MEME --alignment {input.input_msa} --tree {input.input_tree_full} --output {output.output} --branches {LABEL}"
#end rule 

rule fade:
    input:
        input_msa = rules.convert_to_protein.output.proteinput_fas,
        input_tree_clade = rules.annotate.output.output_clade_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.FADE.json")
    shell:
        "hyphy FADE --alignment {input.input_msa} --tree {input.input_tree_clade} --output {output.output} --branches {LABEL}"
#end rule 

rule cfel:
    input:
        input_msa = rules.combine.output.output,
        input_tree_clade = rules.annotate.output.output_clade_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.CFEL.json")
    shell:
        "mpirun -np {PPN} HYPHYMPI contrast-fel --alignment {input.input_msa} --tree {input.input_tree_clade} --output {output.output} --branch-set {LABEL} --branch-set Reference"
#end rule 

#==============================================================================
# MH Methods
#==============================================================================

rule fmm:
    input:
        input_msa = rules.combine.output.output,
        input_tree_clade = rules.annotate.output.output_clade_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.FMM.json")
    shell:
        "mpirun -np {PPN} HYPHYMPI FMM --alignment {input.input_msa} --tree {input.input_tree_clade} --output {output.output} --triple-islands Yes"
#end rule 

#==============================================================================
# Generate summary and annotation report
#==============================================================================
rule generate_report:
    input:
        expand(os.path.join(OUTDIR, "{GENE}.SLAC.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.FEL.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.MEME.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.MEME-full.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.PRIME.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.FADE.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.BUSTEDS.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.RELAX.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.CFEL.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.FMM.json"), GENE=genes)
    output:
        SUMMARY_JSON = os.path.join(OUTDIR, LABEL + "_summary.json"),
        ANNOTATION_JSON = os.path.join(OUTDIR, LABEL + "_annotation.json")
    shell:
         "bash scripts/process_json.sh {BASEDIR} {LABEL}"
#end rule

#==============================================================================
# End of file
#==============================================================================
