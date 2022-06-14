import sys, csv, json
from   operator import itemgetter, attrgetter

from Bio import SeqIO

sequences = {}

#input_file = snakemake.params.input
#output_file = snakemake.params.output

input_file  = sys.argv[1]
output_file = sys.argv[2]



#with open(sys.argv[1]) as handle:))
with open(input_file) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        S = str (record.seq)
        if not S in sequences:
            sequences[S] = set()
        #end if    
        sequences[S].add (record.id)
    #end for
#end with

output = []

for seq, ids in sequences.items():
    m = [k for k in ids]
    output.append ({
        'size' : len (ids),
        'members'  : m,
        'centroid' : ">%s||%d\n%s" % (m[0], len (ids), seq)
    })       
#end for

# Change
#json.dump (output, sys.stdout)

of = open(output_file, "w")
json.dump (output, of)
of.close()


# End of file
