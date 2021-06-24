# Snakefile for SARS-CoV-2 clades analysis.
# @Author: Alexander Lucaci, Jordan Zehr, Stephen Shank

# HELPFUL
### !!! double `{` becomes escapes --> {{this}} 
### use log files to track status of jobs?
### snakemake --forceall --dag | dot -Tpdf > dag.pdf (build the dag)
### https://snakemake.readthedocs.io/en/stable/executing/cli.html

# Imports
import os
import sys
import json
import csv
from pathlib import Path

# Declares ------------------------------------------------------------
with open("snakemake_config.json", "r") as in_sc:
  config = json.load(in_sc)

with open("cluster.json", "r") as in_c:
  cluster = json.load(in_c)

# User defined settings -----------------------------------------------
BASEDIR = os.getcwd()
print(f'BASEDIR: {BASEDIR}')

# Which clades are you analyzing?
# It will look for this as a folder within data/, so make sure they match
LABEL = config["LABEL"]
GISAID_WG = config["GISAID_WG"] 

if not os.path.exists(os.path.join(BASEDIR, "data/" + LABEL)):
  os.mkdir(os.path.join(BASEDIR, "data/" + LABEL))
  print(f'MAKING DIR data/{LABEL}')
else:
  print(f'RUNNING ANALYSES FOR LABEL: {LABEL}')

# The script will also look for this specific whole genome fasta
# within data/{LABEL}/
### !!! need to check to see if this is a file, and alert if not there !!! ###
INPUT_WG = os.path.join(BASEDIR, "data/" + LABEL + "/" + GISAID_WG)
print(f"INPUT_WG: {INPUT_WG}")

# End -- User defined settings ----------------------------------------

genes = ["leader", "nsp2", "nsp3", "nsp4", "3C", "nsp6", "nsp7", "nsp8", "nsp9", "nsp10", "helicase", "exonuclease", "endornase", "S", "E", "M", "N", "ORF3a", "ORF6", "ORF7a", "ORF8" ,"RdRp", "methyltransferase"]
# for debugging or single gene analyses
#genes = ["S"]

# Reference sequence dirs
REF_SEQ_DIR = os.path.join(BASEDIR, "data/ReferenceSeq")
REF_ALN_DIR = os.path.join(BASEDIR, "data/ReferenceSetViPR")

# Set output directory
OUTDIR = os.path.join(BASEDIR, "results/" + LABEL)

# Create output dir.
Path(OUTDIR).mkdir(parents=True, exist_ok=True)
#OUTPUT_WG_FA = os.path.join(OUTDIR, GISAID_WG+".fa")

# Settings, these can be passed in or set in a config.json type file
PPN = cluster["__default__"]["ppn"] 

# Rule All ------------------------------------------------------------
## not working for now
#expand(os.path.join(OUTDIR, "{GENE}.BGM.json"), GENE=genes),
rule all:
    input:
        os.path.join(OUTDIR, GISAID_WG+".fa"),
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
        expand(os.path.join(OUTDIR, "{GENE}.ABSREL.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.BUSTEDS.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.RELAX.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.CFEL.json"), GENE=genes)
#end rule -- all

# Rules -- Main analysis ----------------------------------------------
## works
rule clean:
    input:
        in_wg = INPUT_WG
    output:
        out_wg = os.path.join(OUTDIR, GISAID_WG+".fa")
    shell:
       "bash scripts/cleaner.sh {input.in_wg} {output.out_wg}"
#end rule -- clean

# PROCESS QUERY SEQUENCES ---
## works
rule bealign_query:
    input:
        in_genome = rules.clean.output.out_wg,
        in_gene_ref_seq = os.path.join(REF_SEQ_DIR, "{GENE}.fas")
    output:
        output = os.path.join(OUTDIR, "{GENE}.query.bam")
    shell:
        "bealign -r {input.in_gene_ref_seq} -m HIV_BETWEEN_F {input.in_genome} {output.output}"
#end rule -- bealign

## works
rule bam2msa_query:
    input:
        in_bam = rules.bealign_query.output.output
    output:
        out_msa = os.path.join(OUTDIR, "{GENE}.query.msa.OG")
    shell:
        "bam2msa {input.in_bam} {output.out_msa}"        
#end rule -- bam2msa_query

## works
rule strike_ambigs_query:
   input:
       in_msa = rules.bam2msa_query.output.out_msa
   output:
       out_strike_ambigs = os.path.join(OUTDIR, "{GENE}.query.msa.SA")
   shell:
      "hyphy scripts/strike-ambigs.bf --alignment {input.in_msa} --output {output.out_strike_ambigs}"
#end rule -- strike_ambigs_query

## works
rule tn93_cluster_query:
    params:
        THRESHOLD_QUERY = config["threshold_query"],
        MAX_QUERY = config["max_query"]
    
    
    input:
        in_msa = rules.strike_ambigs_query.output.out_strike_ambigs
    output:
        out_fasta = os.path.join(OUTDIR, "{GENE}.query.compressed.fas"),
        out_json = os.path.join(OUTDIR, "{GENE}.query.json")
    shell:
        "python3 scripts/tn93_cluster.py --input {input.in_msa} --output_fasta {output.out_fasta} --output_json {output.out_json} --threshold {params.THRESHOLD_QUERY} --max_retain {params.MAX_QUERY}"
#end rule tn93_cluster_query

# Do the above for Reference sequences. --------------------------------------------
## works
rule bealign_ref:
    input:
        in_genome_ref = os.path.join(REF_ALN_DIR, "sequences.{GENE}_nuc.compressed.fas"),
        in_gene_ref_seq = os.path.join(REF_SEQ_DIR, "{GENE}.fas")
    output:
        output = os.path.join(OUTDIR, "{GENE}.reference.bam")
    shell:
        "bealign -r {input.in_gene_ref_seq} -m HIV_BETWEEN_F -K {input.in_genome_ref} {output.output}"
#end rule -- bealign_ref

## works
rule bam2msa_ref:
    input:
        in_bam = rules.bealign_ref.output.output
    output:
        out_msa = os.path.join(OUTDIR, "{GENE}.reference.msa.OG")
    shell:
        "bam2msa {input.in_bam} {output.out_msa}"
#end rule -- bam2msa_ref

## works
rule strike_ambigs_ref:
   input:
       in_msa = rules.bam2msa_ref.output.out_msa
   output:
       out_strike_ambigs = os.path.join(OUTDIR, "{GENE}.reference.msa.SA")
   shell:
      "hyphy scripts/strike-ambigs.bf --alignment {input.in_msa} --output {output.out_strike_ambigs}"
#end rule -- strike_ambigs_ref

## works 
rule tn93_cluster_ref:
    params:
        THRESHOLD_REF = config["threshold_ref"],
        MAX_REF = config["max_ref"]

    input:
        in_msa = rules.strike_ambigs_ref.output.out_strike_ambigs,
        in_gene_ref_seq = os.path.join(REF_SEQ_DIR, "{GENE}.fas")
    output:
        out_fasta = os.path.join(OUTDIR, "{GENE}.reference.compressed.fas"),
        out_json = os.path.join(OUTDIR, "{GENE}.reference.json")
    shell:
        "python3 scripts/tn93_cluster.py --input {input.in_msa} --output_fasta {output.out_fasta} --output_json {output.out_json} --threshold {params.THRESHOLD_REF} --max_retain {params.MAX_REF} --reference_seq {input.in_gene_ref_seq}"
#end rule tn93_cluster_ref

# Combine them, the alignment ----------------------------------------------------
## works 
rule combine:
    params:
      THRESHOLD_QUERY = config["threshold_query"]

    input:
        in_compressed_fas = rules.tn93_cluster_query.output.out_fasta,
        in_msa = rules.tn93_cluster_ref.output.out_fasta,
	      in_gene_ref_seq = os.path.join(REF_SEQ_DIR, "{GENE}.fas")
    output:
        output = os.path.join(OUTDIR, "{GENE}.combined.fas")
    shell:
        "python3 scripts/combine.py --input {input.in_compressed_fas} -o {output.output} --threshold {params.THRESHOLD_QUERY} --msa {input.in_msa} --reference_seq {input.in_gene_ref_seq}"
#end rule -- combine

# Convert to protein
## works
rule convert_to_protein:
    input:
        combined_fas = rules.combine.output.output
    output:
        protein_fas = os.path.join(OUTDIR, "{GENE}.AA.fas")
    shell:
        "hyphy conv Universal 'Keep Deletions' {input.combined_fas} {output.protein_fas}"
#end rule -- convert_to_protein

# Combined ML Tree
## works
rule raxml:
    params:
        threads = PPN
    input:
        combined_fas = rules.combine.output.output
    output:
        combined_tree = os.path.join(OUTDIR, "{GENE}.combined.fas.raxml.bestTree")
    shell:
        "raxml-ng --model GTR --msa {input.combined_fas} --threads {params.THREADS} --tree pars{{3}} --force"
#end rule -- raxml

## works ~~ a bit clunky
rule annotate:
    input:
       in_tree = rules.raxml.output.combined_tree,
       in_compressed_fas = rules.tn93_cluster_query.output.out_fasta
    output:
       out_int_tree = os.path.join(OUTDIR, "{GENE}.int.nwk"),
       out_clade_tree = os.path.join(OUTDIR, "{GENE}.clade.nwk"),
       out_full_tree = os.path.join(OUTDIR, "{GENE}.full.nwk")
    shell:
       "bash scripts/annotate.sh {input.in_tree} 'REFERENCE' {input.in_compressed_fas} {LABEL} {BASEDIR}"
#end rule annotate
    
######################################################################
#---------------------Selection analyses ----------------------------#
######################################################################

rule slac:
    input:
        in_msa = rules.combine.output.output,
        in_tree = rules.annotate.output.out_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.SLAC.json")
    shell:
        # also has --samples 0 but I exlude it here.
        "mpirun -np {PPN} HYPHYMPI  SLAC --alignment {input.in_msa} --tree {input.in_tree} --output {output.output}"
#end rule -- slac

rule bgm:
    input:
        in_msa = rules.combine.output.output,
        in_tree = rules.annotate.output.out_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.combined.fas.BGM.json")
    shell:
        "hyphy BGM --alignment {input.in_msa} --tree {input.in_tree} --output {output.output} --branches {LABEL}"
#end rule -- bgm

rule fel:
    input:
        in_msa = rules.combine.output.output,
        in_tree = rules.annotate.output.out_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.FEL.json")
    shell:
        "mpirun -np {PPN} HYPHYMPI FEL --alignment {input.in_msa} --tree {input.in_tree} --output {output.output} --branches {LABEL}"
#end rule -- fel

rule meme:
    input:
        in_msa = rules.combine.output.output,
        in_tree = rules.annotate.output.out_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.MEME.json")
    shell:
        "mpirun -np {PPN} HYPHYMPI MEME --alignment {input.in_msa} --tree {input.in_tree} --output {output.output} --branches {LABEL}"
#end rule -- MEME

# These are exlcuded from Minimal run (not implemented)
rule absrel:
    input:
        in_msa = rules.combine.output.output,
        in_tree = rules.annotate.output.out_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.ABSREL.json")
    shell:
        "hyphy ABSREL --alignment {input.in_msa} --tree {input.in_tree} --output {output.output} --branches {LABEL}"
#end rule -- absrel

## why not bustedS? ##
rule busted:
    input:
        in_msa = rules.combine.output.output,
        in_tree_clade = rules.annotate.output.out_clade_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.BUSTEDS.json")
    shell:
        "mpirun -np {PPN} HYPHYMPI BUSTED --alignment {input.in_msa} --tree {input.in_tree_clade} --output {output.output} --branches {LABEL} --starting-points 10"
#end rule -- busted

rule relax:
    input:
        in_msa = rules.combine.output.output,
        in_tree_clade = rules.annotate.output.out_clade_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.RELAX.json")
    shell:
        "hyphy RELAX --alignment {input.in_msa} --models Minimal --tree {input.in_tree_clade} --output {output.output} --test {LABEL} --reference Reference --starting-points 10 --srv Yes"
#end rule -- relax
# End exclusion --

rule prime:
    input:
        in_msa = rules.combine.output.output,
        in_tree = rules.annotate.output.out_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.PRIME.json")
    shell:
        "mpirun -np {PPN} HYPHYMPI PRIME --alignment {input.in_msa} --tree {input.in_tree} --output {output.output} --branches {LABEL}"
#end rule -- prime

rule meme_full:
    input:
        in_msa = rules.combine.output.output,
        in_tree_full = rules.annotate.output.out_full_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.MEME-full.json")
    shell:
        "mpirun -np {PPN} HYPHYMPI --alignment {input.in_msa} --tree {input.in_tree_full} --output {output.output} --branches {LABEL}"
#end rule -- meme_full

rule fade:
    input:
        in_msa = rules.convert_to_protein.output.protein_fas,
        in_tree_clade = rules.annotate.output.out_clade_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.FADE.json")
    shell:
        "hyphy FADE --alignment {input.in_msa} --tree {input.in_tree_clade} --output {output.output} --branches {LABEL}"
#end rule -- fade

# cFEL
rule cfel:
    input:
        in_msa = rules.combine.output.output,
        in_tree_clade = rules.annotate.output.out_clade_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.CFEL.json")
    shell:
        "hyphy contrast-fel --alignment {input.in_msa} --tree {input.in_tree_clade} --output {output.output} --branch-set {LABEL} --branch-set Reference"
#end rule -- cfel

# End of file
