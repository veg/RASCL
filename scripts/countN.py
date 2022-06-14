# Imports
import sys, csv, json, datetime
from collections import Counter
from Bio import SeqIO

# Declares
variants = {}
byPosition = []
N          = 0
dates      = Counter()

# Input parameters
#input_msa = snakemake.params.input_msa
#input_all = snakemake.params.input_all
#input_uniq = snakemake.params.in_uniq
#output1 = snakemake.params.output1
#output2 = snakemake.params.output2

# Input parameters
input_msa  = sys.argv[1]
input_all  = sys.argv[2]
input_uniq = sys.argv[3]
output1    = sys.argv[4]
output2    = sys.argv[5]

"""
        input_msa = rules.strike_ambigs_query.output.out_strike_ambigs,
        input_all = rules.cluster_processor_t0.output.output,
        in_uniq =  rules.cluster_processor_consensus.output.output,

"""

# Main ---
#with open(sys.argv[1]) as handle:
with open(input_msa) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        N += 1
        try:
            D = datetime.datetime.strptime (record.id.split ("|")[-1], "%Y-%m-%d").strftime ("%Y%m%d")
        except Exception as e:
            D = None
        #end try

        dates [D] += 1

        S = str (record.seq)
        for i in range (0,len(record.seq),3):
            if N == 1:
                byPosition.append (Counter())
            codon = S[i:i+3]
            byPosition[i//3][codon] += 1
        #end inner for
    #end outer for
#end with

variants ['N']      = N            
variants ['counts'] = byPosition
variants ['dates'] = dates
hbp        = []
H          = 0

mapping = {
    'A' : 0,
    'C' : 1,
    'G' : 2,
    'T' : 3,
    'N' : 4,
    '-' : 5
}

dotplot = []

#with open(sys.argv[2]) as handle:
with open(input_uniq) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        H += 1
        S = str (record.seq)
        
        dotplot.append ([mapping[k] if k in mapping else 6 for k in S])
        
        for i in range (0,len(record.seq),3):
            if H == 1:
                hbp.append (Counter())
            codon = S[i:i+3]
            hbp[i//3][codon] += 1
        #end for
    #end for
#end with          

variants ['H']      = H            
variants ['haplotypes'] = hbp
hbpa        = []
HA          = 0

#with open(sys.argv[3]) as handle:
with open(input_all) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        HA += 1
        S = str (record.seq)
               
        for i in range (0,len(record.seq),3):
            if HA == 1:
                hbpa.append (Counter())
            codon = S[i:i+3]
            hbpa[i//3][codon] += 1
        #end for
    #end for
#end with


variants ['HA']             = HA            
variants ['all-haplotypes'] = hbpa
    
# Output files

#json.dump (variants, sys.stdout)
#json.dump (dotplot, sys.stderr)

with open(output1, 'w') as f:
    json.dump (variants, f)


with open(output2, 'w') as fh:
    json.dump (dotplot, fh)



# End of file



        
                
