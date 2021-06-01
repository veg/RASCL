# Snakemake pipeline for SARS-CoV-2 clades analysis.
# @Author: Alexander Lucaci


# @Usage on hpc:
# (Current) snakemake -s Snakefile --cluster "qsub -V -l nodes=1:ppn=8 -q epyc2 -l walltime=999:00:00" --jobs 50 all --rerun-incomplete --keep-going

# Imports
import os
import sys
import json
import csv
from pathlib import Path

# Declares ------------------------------------------------------------

# User defined settings -----------------------------------------------
BASEDIR = "/home/aglucaci/SARS-CoV-2_Clades"
# Which clades are you analyzing?
# It will look for this as a folder within data/, so make sure they match
LABEL = "B-1-617-2"
# The script will also look for this specific whole genome fasta
# within data/{LABEL}/
GISAID_WG = "gisaid_hcov-19_2021_05_14_19.fasta"
# End -- User defined settings ----------------------------------------

INPUT_WG = os.path.join(BASEDIR, "data/" + LABEL + "/" + GISAID_WG)
print("# INPUT GISAID FILE:", INPUT_WG)

genes = ["leader", "nsp2", "nsp3", "nsp4", "3C", "nsp6", "nsp7", "nsp8", "nsp9", "nsp10", "helicase", "exonuclease", "endornase", "S", "E", "M", "N", "ORF3a", "ORF6", "ORF7a", "ORF8" ,"RdRp", "methyltransferase"]
# for debugging or single gene analyses
#genes = ["S"]

# Load config.json settings -------------------------------------------
REF_SEQ_DIR = os.path.join(BASEDIR, "data/ReferenceSeq")
REF_ALN_DIR = os.path.join(BASEDIR, "data/ReferenceSetViPR")

# Software dependencies
# This will be moved to a config.json 
# Or config.yaml type file.
SNAKEMAKE = "/home/aglucaci/anaconda3/bin/snakemake"
BEALIGN = "/home/aglucaci/anaconda3/bin/bealign"
BAM2MSA = "/usr/local/bin/bam2msa"
HYPHYMP = "/home/aglucaci/hyphy-develop/hyphy"

# Modify to parallelize
# NP setting needs to match snakemake evokation 
HYPHYMPI = "mpirun -np 8 /home/aglucaci/hyphy-develop/HYPHYMPI"

RES = "/home/aglucaci/hyphy-develop/res"
PYTHON = "/home/aglucaci/anaconda3/bin/python"
CLEANER = os.path.join(BASEDIR, "scripts/cleaner.sh")
ANNOTATE = os.path.join(BASEDIR, "scripts/annotate.sh")
TN93_CLUSTER = os.path.join(BASEDIR, "scripts/tn93_cluster.py")
COMBINE = os.path.join(BASEDIR, "scripts/combine.py")
RAXML = "/usr/local/bin/raxml-ng"

# Batch files for custom HyPhy procedures
STRIKE_AMBIGS_BF = os.path.join(BASEDIR, "scripts/strike-ambigs.bf")

# Annotate the tree
ANNOTATOR = os.path.join(BASEDIR, "scripts/annotator.bf")

# Set output directory
OUTDIR = os.path.join(BASEDIR, "results/" + LABEL)

# Create output dir.
Path(OUTDIR).mkdir(parents=True, exist_ok=True)
OUTPUT_WG_FA = os.path.join(OUTDIR, GISAID_WG+".fa")

# Settings, these can be passed in or set in a config.json type file
MAX_REF = 200
MAX_QUERY = 500
THRESHOLD_QUERY = 0.0005
THRESHOLD_REF = 0.0005
THREADS = 8

# Rule All ------------------------------------------------------------
rule all:
    input:
        OUTPUT_WG_FA,
        expand(os.path.join(OUTDIR, "{GENE}.query.bam"), GENE=genes),                     
        expand(os.path.join(OUTDIR, "{GENE}.query.msa.OG"), GENE=genes),                     
        expand(os.path.join(OUTDIR, "{GENE}.query.msa.SA"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.query.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.query.compressed.fas"), GENE=genes),
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
        expand(os.path.join(OUTDIR, "{GENE}.SLAC.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.BGM.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.FEL.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.MEME.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.MEME-full.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.PRIME.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.FADE.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.ABSREL.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.BUSTED.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.RELAX.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.CFEL.json"), GENE=genes)
#end rule -- all

# Rules -- Main analysis ----------------------------------------------
rule clean:
    input:
        INPUT = INPUT_WG
    output:
        output = OUTPUT_WG_FA
    shell:
       "bash " + CLEANER + " {input.INPUT} {output.output}"
#end rule -- clean

# PROCESS QUERY SEQUENCES ---
rule bealign_query:
    input:
        in_genome = rules.clean.output.output,
        in_gene_ref_seq = os.path.join(REF_SEQ_DIR, "{GENE}.fas")
    output:
        output = os.path.join(OUTDIR, "{GENE}.query.bam")
    shell:
        "{BEALIGN} -r {input.in_gene_ref_seq} -m HIV_BETWEEN_F {input.in_genome} {output.output}"
#end rule -- bealign

rule bam2msa_query:
    input:
        in_bam = rules.bealign_query.output.output
    output:
        out_msa = os.path.join(OUTDIR, "{GENE}.query.msa.OG")
    shell:
        "{BAM2MSA} {input.in_bam} {output.out_msa}"        
#end rule -- bam2msa_query

rule strike_ambigs_query:
   input:
       in_msa = rules.bam2msa_query.output.out_msa
   output:
       out_strike_ambigs = os.path.join(OUTDIR, "{GENE}.query.msa.SA")
   shell:
      "{HYPHYMP} LIBPATH={RES} {STRIKE_AMBIGS_BF} --alignment {input.in_msa} --output {output.out_strike_ambigs}"
#end rule -- strike_ambigs_query

rule tn93_cluster_query:
    input:
        in_msa = rules.strike_ambigs_query.output.out_strike_ambigs
    output:
        out_fasta = os.path.join(OUTDIR, "{GENE}.query.compressed.fas"),
        out_json = os.path.join(OUTDIR, "{GENE}.query.json")
    shell:
        "{PYTHON} {TN93_CLUSTER} --input {input.in_msa} --output_fasta {output.out_fasta} --output_json {output.out_json} --threshold {THRESHOLD_QUERY} --max_retain {MAX_QUERY}"
#end rule tn93_cluster_query

# Do the above for Reference sequences. --------------------------------------------
rule bealign_ref:
    input:
        in_genome_ref = os.path.join(REF_ALN_DIR, "sequences.{GENE}_nuc.compressed.fas"),
        in_gene_ref_seq = os.path.join(REF_SEQ_DIR, "{GENE}.fas")
    output:
        output = os.path.join(OUTDIR, "{GENE}.reference.bam")
    shell:
        "{BEALIGN} -r {input.in_gene_ref_seq} -m HIV_BETWEEN_F -K {input.in_genome_ref} {output.output}"
#end rule -- bealign_ref

rule bam2msa_ref:
    input:
        in_bam = rules.bealign_ref.output.output
    output:
        out_msa = os.path.join(OUTDIR, "{GENE}.reference.msa.OG")
    shell:
        "{BAM2MSA} {input.in_bam} {output.out_msa}"
#end rule -- bam2msa_ref

rule strike_ambigs_ref:
   input:
       in_msa = rules.bam2msa_ref.output.out_msa
   output:
       out_strike_ambigs = os.path.join(OUTDIR, "{GENE}.reference.msa.SA")
   shell:
      "{HYPHYMP} LIBPATH={RES} {STRIKE_AMBIGS_BF} --alignment {input.in_msa} --output {output.out_strike_ambigs}"
#end rule -- strike_ambigs_ref

rule tn93_cluster_ref:
    input:
        in_msa = rules.strike_ambigs_ref.output.out_strike_ambigs,
        in_gene_ref_seq = os.path.join(REF_SEQ_DIR, "{GENE}.fas")
    output:
        out_fasta = os.path.join(OUTDIR, "{GENE}.reference.compressed.fas"),
        out_json = os.path.join(OUTDIR, "{GENE}.reference.json")
    shell:
        "{PYTHON} {TN93_CLUSTER} --input {input.in_msa} --output_fasta {output.out_fasta} --output_json {output.out_json} --threshold {THRESHOLD_REF} --max_retain {MAX_REF} --reference_seq {input.in_gene_ref_seq}"
#end rule tn93_cluster_ref

# Combine them, the alignment ----------------------------------------------------
rule combine:
    input:
        in_compressed_fas = rules.tn93_cluster_query.output.out_fasta,
        in_msa = rules.tn93_cluster_ref.output.out_fasta,
	in_gene_ref_seq = os.path.join(REF_SEQ_DIR, "{GENE}.fas")
    output:
        output = os.path.join(OUTDIR, "{GENE}.combined.fas")
    shell:
        "{PYTHON} {COMBINE} --input {input.in_compressed_fas} -o {output.output} --threshold {THRESHOLD_QUERY} --msa {input.in_msa} --reference_seq {input.in_gene_ref_seq}"
#end rule -- combine

# Convert to protein
rule convert_to_protein:
    input:
        combined_fas = rules.combine.output.output
    output:
        protein_fas = os.path.join(OUTDIR, "{GENE}.AA.fas")
    shell:
        "{HYPHYMP} LIBPATH={RES} conv Universal 'Keep Deletions' {input.combined_fas} {output.protein_fas}"
#end rule -- convert_to_protein

# Combined ML Tree
rule raxml:
    input:
        combined_fas = rules.combine.output.output
    output:
        combined_tree = os.path.join(OUTDIR, "{GENE}.combined.fas.raxml.bestTree")
    shell:
        # --tree pars{3} fucks up here due to the curly brackets
        "{RAXML} --model GTR --msa {input.combined_fas} --threads " + str(THREADS) + " --force"
#end rule -- raxml

#rule annotate:
#    input:
#        in_tree = rules.raxml.output.combined_tree,
#        in_compressed_fas = rules.tn93_cluster_query.output.out_fasta
#    output:
#        out_int_tree = os.path.join(OUTDIR, "{GENE}.int.nwk"),
#        out_clade_tree = os.path.join(OUTDIR, "{GENE}.clade.nwk"),
#        out_full_tree = os.path.join(OUTDIR, "{GENE}.full.nwk"),
#	 out_prefix = os.path.join(OUTDIR, "{GENE}.")
#    shell:
#        "{HYPHYMP} LIBPATH={RES} {ANNOTATOR} {input.in_tree} 'REFERENCE' {input.in_compressed_fas} {LABEL} {output.out_prefix}"
#end rule -- annotate

rule annotate:
    input:
       in_tree = rules.raxml.output.combined_tree,
       in_compressed_fas = rules.tn93_cluster_query.output.out_fasta
    output:
       out_int_tree = os.path.join(OUTDIR, "{GENE}.int.nwk"),
       out_clade_tree = os.path.join(OUTDIR, "{GENE}.clade.nwk"),
       out_full_tree = os.path.join(OUTDIR, "{GENE}.full.nwk")
    shell:
       "bash {ANNOTATE} {input.in_tree} 'REFERENCE' {input.in_compressed_fas} {LABEL}"
#end rule annotate
    

# #####################################################################
# Selection analyses -------------------------------------------------#
# #####################################################################
rule slac:
    input:
        in_msa = rules.combine.output.output,
        in_tree = rules.annotate.output.out_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.SLAC.json")
    shell:
        # also has --samples 0 but I exlude it here.
        "{HYPHYMP} LIBPATH={RES} SLAC --alignment {input.in_msa} --tree {input.in_tree} --output {output.output}"
#end rule -- slac

rule bgm:
    input:
        in_msa = rules.combine.output.output,
        in_tree = rules.annotate.output.out_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.BGM.json")
    shell:
        "{HYPHYMP} LIBPATH={RES} BGM --alignment {input.in_msa} --tree {input.in_tree} --output {output.output} --branches {LABEL}"
#end rule -- bgm

rule fel:
    input:
        in_msa = rules.combine.output.output,
        in_tree = rules.annotate.output.out_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.FEL.json")
    shell:
        "{HYPHYMPI} LIBPATH={RES} FEL --alignment {input.in_msa} --tree {input.in_tree} --output {output.output} --branches {LABEL}"
#end rule -- fel

rule meme:
    input:
        in_msa = rules.combine.output.output,
        in_tree = rules.annotate.output.out_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.MEME.json")
    shell:
        "{HYPHYMPI} LIBPATH={RES} MEME --alignment {input.in_msa} --tree {input.in_tree} --output {output.output} --branches {LABEL}"
#end rule -- MEME

# These are exlcuded from Minimal run (not implemented)
rule absrel:
    input:
        in_msa = rules.combine.output.output,
        in_tree = rules.annotate.output.out_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.ABSREL.json")
    shell:
        "{HYPHYMP} LIBPATH={RES} ABSREL --alignment {input.in_msa} --tree {input.in_tree} --output {output.output} --branches {LABEL}"
#end rule -- absrel

rule busted:
    input:
        in_msa = rules.combine.output.output,
        in_tree_clade = rules.annotate.output.out_clade_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.BUSTED.json")
    shell:
        "{HYPHYMP} LIBPATH={RES} BUSTED --alignment {input.in_msa} --tree {input.in_tree_clade} --output {output.output} --branches {LABEL} --starting-points 10"
#end rule -- busted

rule relax:
    input:
        in_msa = rules.combine.output.output,
        in_tree_clade = rules.annotate.output.out_clade_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.RELAX.json")
    shell:
        "{HYPHYMP} LIBPATH={RES} RELAX --alignment {input.in_msa} --models Minimal --tree {input.in_tree_clade} --output {output.output} --test {LABEL} --reference Reference --starting-points 10 --srv Yes"
#end rule -- relax
# End exclusion --

rule prime:
    input:
        in_msa = rules.combine.output.output,
        in_tree = rules.annotate.output.out_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.PRIME.json")
    shell:
        "{HYPHYMPI} LIBPATH={RES} PRIME --alignment {input.in_msa} --tree {input.in_tree} --output {output.output} --branches {LABEL}"
#end rule -- prime

rule meme_full:
    input:
        in_msa = rules.combine.output.output,
        in_tree_full = rules.annotate.output.out_full_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.MEME-full.json")
    shell:
        "{HYPHYMPI} LIBPATH={RES} MEME --alignment {input.in_msa} --tree {input.in_tree_full} --output {output.output} --branches {LABEL}"
#end rule -- meme_full

rule fade:
    input:
        in_msa = rules.convert_to_protein.output.protein_fas,
        in_tree_clade = rules.annotate.output.out_clade_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.FADE.json")
    shell:
        "{HYPHYMP} LIBPATH={RES} FADE --alignment {input.in_msa} --tree {input.in_tree_clade} --output {output.output} --branches {LABEL}"
#end rule -- fade

# cFEL
rule cfel:
    input:
        in_msa = rules.combine.output.output,
        in_tree_clade = rules.annotate.output.out_clade_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.CFEL.json")
    shell:
        "{HYPHYMPI} LIBPATH={RES} contrast-fel --alignment {input.in_msa} --tree {input.in_tree_clade} --output {output.output} --branch-set {LABEL} --branch-set Reference"
#end rule -- cfel

# End of file
