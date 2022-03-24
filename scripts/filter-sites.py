import sys, csv, json
from collections import Counter
from Bio import SeqIO


# Arguments
#ranges = [int(k) for k in sys.argv[2].split (',')]

#start = sys.argv[2]
#end = sys.argv[3]

input_file = snakemake.params.in_wg

start = snakemake.params.start_site
end = snakemake.params.end_site


#output_file = sys.argv[1]
output_file = snakemake.params.output


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


        
                
