# Snakefile for SARS-CoV-2 Clades analysis (RASCL).
# 2021
# @Author: Alexander Lucaci, Jordan Zehr, Stephen Shank

# Imports -------------------------------------------------------------
import os
import sys
import json
import csv
from pathlib import Path

# Declares ------------------------------------------------------------
with open("snakemake_config.json", "r") as in_sc:
  config = json.load(in_sc)
#end with

with open("cluster.json", "r") as in_c:
  cluster = json.load(in_c)
#end with

# User settings -------------------------------------------------------
BASEDIR = os.getcwd()

print(f'We are operating out of base directory: {BASEDIR}')

# Which clades are you analyzing?
LABEL = config["LABEL"] # add assert
WholeGenomeSeqs = config["WholeGenomeSeqs"] # add assert

# The script will also look for this specific whole genome fasta within data/{LABEL}/
#INPUT_WG = os.path.join(BASEDIR, "data", LABEL, WholeGenomeSeqs)
INPUT_WG = os.path.join(BASEDIR, "data", WholeGenomeSeqs)
print(f"Input whole genome fasta: {INPUT_WG}")

# End -- User defined settings ----------------------------------------

genes = ["leader", "nsp2", "nsp3", "nsp4", "3C", "nsp6", "nsp7", "nsp8", "nsp9", "nsp10", "helicase", "exonuclease", "endornase", "S", "E", "M", "N", "ORF3a", "ORF6", "ORF7a", "ORF8" ,"RdRp", "methyltransferase"]

# For debugging or single gene analyses
#genes = ["S"]

# Basename of INPUT_WG (used in rule clean)
INPUT_WG_basename = os.path.basename(INPUT_WG)

# Reference sequence dirs
REF_SEQ_DIR = os.path.join(BASEDIR, "data", "ReferenceSeq")
REF_ALN_DIR = os.path.join(BASEDIR, "data", "ReferenceSetViPR")

# Set output directory
OUTDIR = os.path.join(BASEDIR, "results", LABEL)

# Create output dir.
Path(os.path.join(BASEDIR,"results")).mkdir(parents=True, exist_ok=True)
Path(OUTDIR).mkdir(parents=True, exist_ok=True)

# Settings, these can be passed in or set in a config.json type file
PPN = cluster["__default__"]["ppn"] 

# Rule All ------------------------------------------------------------
rule all:
    input:
        os.path.join(OUTDIR, INPUT_WG_basename + ".fa"),
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
        os.path.join(OUTDIR, LABEL + "_summary.json"),
        os.path.join(OUTDIR, LABEL + "_annotation.json")
#end rule -- all

# Rules -- Main analysis ----------------------------------------------
rule clean:
    input:
        in_wg = INPUT_WG
    output:
        out_wg = os.path.join(OUTDIR, INPUT_WG_basename + ".fa")
    shell:
       "bash scripts/cleaner.sh {input.in_wg} {output.out_wg}"
#end rule -- clean

# PROCESS QUERY SEQUENCES
rule bealign_query:
    input:
        in_genome = rules.clean.output.out_wg,
        in_gene_ref_seq = os.path.join(REF_SEQ_DIR, "{GENE}.fas")
    output:
        output = os.path.join(OUTDIR, "{GENE}.query.bam")
    shell:
        "bealign -r {input.in_gene_ref_seq} -m HIV_BETWEEN_F {input.in_genome} {output.output}"
#end rule -- bealign

rule bam2msa_query:
    input:
        in_bam = rules.bealign_query.output.output
    output:
        out_msa = os.path.join(OUTDIR, "{GENE}.query.msa.OG")
    shell:
        "bam2msa {input.in_bam} {output.out_msa}"        
#end rule -- bam2msa_query

rule strike_ambigs_query:
   input:
       in_msa = rules.bam2msa_query.output.out_msa
   output:
       out_strike_ambigs = os.path.join(OUTDIR, "{GENE}.query.msa.SA")
   conda: 'environment.yml'
   shell:
      "hyphy scripts/strike-ambigs.bf --alignment {input.in_msa} --output {output.out_strike_ambigs}"
#end rule -- strike_ambigs_query

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
#end rule tn93_cluster_query`

# Do the above for background sequences. --------------------------------------------
rule bealign_ref:
    input:
        in_genome_ref = os.path.join(REF_ALN_DIR, "sequences.{GENE}_nuc.compressed.fas"),
        in_gene_ref_seq = os.path.join(REF_SEQ_DIR, "{GENE}.fas")
    output:
        output = os.path.join(OUTDIR, "{GENE}.reference.bam")
    shell:
        "bealign -r {input.in_gene_ref_seq} -m HIV_BETWEEN_F -K {input.in_genome_ref} {output.output}"
#end rule -- bealign_ref

rule bam2msa_ref:
    input:
        in_bam = rules.bealign_ref.output.output
    output:
        out_msa = os.path.join(OUTDIR, "{GENE}.reference.msa.OG")
    shell:
        "bam2msa {input.in_bam} {output.out_msa}"
#end rule -- bam2msa_ref

rule strike_ambigs_ref:
   input:
       in_msa = rules.bam2msa_ref.output.out_msa
   output:
       out_strike_ambigs = os.path.join(OUTDIR, "{GENE}.reference.msa.SA")
   conda: 'environment.yml'
   shell:
      "hyphy scripts/strike-ambigs.bf --alignment {input.in_msa} --output {output.out_strike_ambigs}"
#end rule -- strike_ambigs_ref

rule tn93_cluster_ref:
    params:
        THRESHOLD_REF = config["threshold_ref"],
        MAX_REF = config["max_ref"],
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
rule combine:
    params:
        THRESHOLD_QUERY = config["threshold_query"]
    input:
        in_compressed_fas = rules.tn93_cluster_query.output.out_fasta,
        in_msa = rules.tn93_cluster_ref.output.out_fasta,
	in_gene_ref_seq = os.path.join(REF_SEQ_DIR, "{GENE}.fas")
    output:
        output = os.path.join(OUTDIR, "{GENE}.combined.fas")
        #output_csv = os.path.join(OUTDIR, "{GENE}.combined.fas.csv")
    conda: 'environment.yml'
    shell:
        "python3 scripts/combine.py --input {input.in_compressed_fas} -o {output.output} --threshold {params.THRESHOLD_QUERY} --msa {input.in_msa} --reference_seq {input.in_gene_ref_seq}"
#end rule -- combine

# Convert to protein
rule convert_to_protein:
    input:
        combined_fas = rules.combine.output.output
    output:
        protein_fas = os.path.join(OUTDIR, "{GENE}.AA.fas")
    conda: 'environment.yml'
    shell:
        "hyphy conv Universal 'Keep Deletions' {input.combined_fas} {output.protein_fas}"
#end rule -- convert_to_protein

# Combined ML Tree
rule raxml:
    params:
        THREADS = PPN
    input:
        combined_fas = rules.combine.output.output
    output:
        combined_tree = os.path.join(OUTDIR, "{GENE}.combined.fas.raxml.bestTree")
    shell:
        "raxml-ng --model GTR --msa {input.combined_fas} --threads {params.THREADS} --tree pars{{3}} --force"
#end rule -- raxml

rule annotate:
    input:
       in_tree = rules.raxml.output.combined_tree,
       in_compressed_fas = rules.tn93_cluster_query.output.out_fasta
    output:
       out_int_tree = os.path.join(OUTDIR, "{GENE}.int.nwk"),
       out_clade_tree = os.path.join(OUTDIR, "{GENE}.clade.nwk"),
       out_full_tree = os.path.join(OUTDIR, "{GENE}.full.nwk")
    conda: 'environment.yml'
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
    conda: 'environment.yml'
    shell:
        "mpirun -np {PPN} HYPHYMPI  SLAC --alignment {input.in_msa} --samples 0 --tree {input.in_tree} --output {output.output}"
#end rule -- slac

rule bgm:
    input:
        in_msa = rules.combine.output.output,
        in_tree = rules.annotate.output.out_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.combined.fas.BGM.json")
    conda: 'environment.yml'
    shell:
        "hyphy BGM --alignment {input.in_msa} --tree {input.in_tree} --output {output.output} --branches {LABEL}"
#end rule -- bgm

rule fel:
    input:
        in_msa = rules.combine.output.output,
        in_tree = rules.annotate.output.out_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.FEL.json")
    conda: 'environment.yml'
    shell:
        "mpirun -np {PPN} HYPHYMPI FEL --alignment {input.in_msa} --tree {input.in_tree} --output {output.output} --branches {LABEL}"
#end rule -- fel

rule meme:
    input:
        in_msa = rules.combine.output.output,
        in_tree = rules.annotate.output.out_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.MEME.json")
    conda: 'environment.yml'
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
    conda: 'environment.yml'
    shell:
        "hyphy ABSREL --alignment {input.in_msa} --tree {input.in_tree} --output {output.output} --branches {LABEL}"
#end rule -- absrel

rule busted:
    input:
        in_msa = rules.combine.output.output,
        in_tree_clade = rules.annotate.output.out_clade_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.BUSTEDS.json")
    conda: 'environment.yml'
    shell:
        "mpirun -np {PPN} HYPHYMPI BUSTED --alignment {input.in_msa} --tree {input.in_tree_clade} --output {output.output} --branches {LABEL} --starting-points 10"
#end rule -- busted

rule relax:
    input:
        in_msa = rules.combine.output.output,
        in_tree_clade = rules.annotate.output.out_clade_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.RELAX.json")
    conda: 'environment.yml'
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
    conda: 'environment.yml'
    shell:
        "mpirun -np {PPN} HYPHYMPI PRIME --alignment {input.in_msa} --tree {input.in_tree} --output {output.output} --branches {LABEL}"
#end rule -- prime

rule meme_full:
    input:
        in_msa = rules.combine.output.output,
        in_tree_full = rules.annotate.output.out_full_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.MEME-full.json")
    conda: 'environment.yml'
    shell:
        "mpirun -np {PPN} HYPHYMPI MEME --alignment {input.in_msa} --tree {input.in_tree_full} --output {output.output} --branches {LABEL}"
#end rule -- meme_full

rule fade:
    input:
        in_msa = rules.convert_to_protein.output.protein_fas,
        in_tree_clade = rules.annotate.output.out_clade_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.FADE.json")
    conda: 'environment.yml'
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
    conda: 'environment.yml'
    shell:
        "hyphy contrast-fel --alignment {input.in_msa} --tree {input.in_tree_clade} --output {output.output} --branch-set {LABEL} --branch-set Reference"
#end rule -- cfel

rule generate_report:
    input:
        expand(os.path.join(OUTDIR, "{GENE}.CFEL.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.FADE.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.MEME-full.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.PRIME.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.ABSREL.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.BUSTEDS.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.RELAX.json"), GENE=genes)
    output:
        SUMMARY_JSON = os.path.join(OUTDIR, LABEL + "_summary.json"),
        ANNOTATION_JSON = os.path.join(OUTDIR, LABEL + "_annotation.json")
    conda: 'environment.yml'
    shell:
         "bash scripts/process_json.sh {BASEDIR} {LABEL}"
#end rule generate_report

