import sys, csv, json
from collections import Counter
from Bio import SeqIO


# Arguments
#ranges = [int(k) for k in sys.argv[2].split (',')]
#start = sys.argv[2]
#end = sys.argv[3]

input_file = snakemake.params.in_wg
#start = snakemake.params.start_site
#end = snakemake.params.end_site
gene = snakemake.params.gene

#output_file = sys.argv[1]
output_file = snakemake.params.output


# Gene map
# Based on coordinates of SC2.
# For example in https://github.com/veg/SARS-CoV-2/blob/master/scripts/extract_genes.sh

#genes = ["leader", "nsp2", "nsp3", "nsp4", "3C", "nsp6", "nsp7", "nsp8", "nsp9", "nsp10", "helicase", "exonuclease", "endornase", "S", "E", "M", "N", "ORF3a", "ORF6", "ORF7a", "ORF8" ,"RdRp", "methyltransferase"]


SC2_GENE_MAP = {"S":      {"start": 20000, "end": 26000},
                "M":      {"start": 25000, "end": 30000},
                "N":      {"start": 26000, "end": 35000},
                "E":      {"start": 25500, "end": 27000},
                "ORF3a":  {"start": 25500, "end": 27000},
                "ORF6":   {"start": 26000, "end": 30000},
                "ORF7a":  {"start": 26000, "end": 35000},
                "ORF8":   {"start": 26000, "end": 35000},
                "leader": {"start": 1,     "end": 1000},
                "nsp2":   {"start": 500,   "end": 3000},
                "nsp3":   {"start": 2000,  "end": 10000},
                "nsp4":   {"start": 8000,  "end": 11000},
                "3C":     {"start": 9000,  "end": 12000},
                "nsp6":   {"start": 10500,  "end": 12500},
                "nsp7":   {"start": 11500,  "end": 12500},
                "nsp8":   {"start": 11500,  "end": 13000},
                "nsp9":   {"start": 12000,  "end": 13500},
                "nsp10":  {"start": 12500,  "end": 14000},
                "RdRp":   {"start": 13000,  "end": 17000},
                "helicase":   {"start": 15500,  "end": 18500},
                "exonuclease":  {"start": 17500,  "end": 20000},
                "endornase":  {"start": 19000,  "end": 21000},
                "methyltransferase":  {"start": 20000,  "end": 22000}
               }

start, end = 0, 0

start = SC2_GENE_MAP[gene]["start"]
end = SC2_GENE_MAP[gene]["end"]

#with open(sys.argv[1]) as handle:
with open(input_file) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        s = str (record.seq)
        if len (s) > 28000:
            #print ('>%s\n%s\n' % (record.id, s[ranges[0]:ranges[1]]))
            
            #print ('>%s\n%s\n' % (record.id, s[start:end]))
            with open(output_file, "a") as fh:
                print ('>%s\n%s\n' % (record.id, s[start:end]), file=fh)
            #end with
            fh.close()
        #end if
    #end for
#end with

handle.close()


        
                
