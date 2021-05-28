# TN93 Cluster script


# Imports -------------------------------------------------------------
import os
import sys
import argparse
import json
import shutil

# Declares
# Argparse here
arguments = argparse.ArgumentParser(description='Cluster an MSA with genetic distance (TN93)')

arguments.add_argument('-i', '--input',            help = 'MSA file to process',                                  required = True, type = str )
arguments.add_argument('-o', '--output_fasta',           help = 'Output json file',                                     required = True, type = str)
arguments.add_argument('-j', '--output_json',           help = 'Output json file',                                     required = True, type = str)
arguments.add_argument('--threshold',              help = 'Distance threshold for clustering query sequences',    required = True, type = float)
#arguments.add_argument('--step',                   help = 'Distance threshold for clustering query sequences',    required = True, type = float)
arguments.add_argument('-m', '--max_retain',    help = 'The maximum number of sequences to retain',               required = True, type = int)
arguments.add_argument('-r', '--reference_seq',    help = 'The maximum number of sequences to retain',               required = False, type = str)

settings = arguments.parse_args()

# Output is actuall the compressed.fas, thats what we want.
# Input is the MSA.
# also pass in {GENE}.{query/reference}.json filename
_ref_seq_name = ""
if settings.reference_seq:
    with open (settings.reference_seq) as fh:
        for l in fh:
            if l[0] == '>':
                _ref_seq_name = l[1:].split (' ')[0].strip()
                break
            #end if
        #end for
    #end with
#end if
            
print ("Reference seq_name %s" % _ref_seq_name)


# files
cluster_json = settings.output_json          #output file # This is actually a json file
compressed_fasta = settings.output_fasta           # .compressed.fas
msa_strike_ambigs = settings.input

# values
threshold = settings.threshold
#step = settings.step
max_toRetain = settings.max_retain

# Need to load this from a config.json file.
task_runners = {}
task_runners['tn93-cluster'] = "/usr/local/bin/tn93-cluster"

# Helper functions -----------------------------------------------------

def run_command (exec, arguments, filename, tag):
    result = os.system (" ".join ([exec] + arguments))
    return os.path.getmtime(filename)
#end method

def cluster_to_fasta (in_file, out_file, ref_seq = None):
    with open (in_file, "r") as fh:
        cluster_json = json.load (fh)
        #print (colored('Running ... converting representative clusters to .FASTA\n', 'cyan'))
        check_uniq = set ()
        print("# Saving to fasta:", out_file)
        with open (out_file, "w") as fh2:
            for c in cluster_json:
                cc = c['centroid'].split ('\n')
                if ref_seq:
                    if ref_seq in c['members']:
                        cc[0] = ">" + ref_seq
                        print ("\n".join (cc), file = fh2)
                        continue
                    #end if
                #end if
                seq_id = cc[0]
                while seq_id in check_uniq:
                        seq_id = cc[0] + '_' + ''.join(random.choices ('0123456789abcdef', k = 10))
                #end while
                check_uniq.add (seq_id)
                print (seq_id,"\n",cc[1], file = fh2)
            #end for
         #end with
    #end with
    return (os.path.getmtime(out_file), len(cluster_json))
#end method

# Main subroutine -----------------------------------------------------

while True:
    input_stamp = run_command (task_runners['tn93-cluster'], ['-f', '-o', cluster_json, '-t', "%g" % threshold, msa_strike_ambigs], cluster_json, "extract representative clusters at threshold %g" % threshold)
    
    
    if _ref_seq_name != "":
        input_stamp, cluster_count = cluster_to_fasta (cluster_json, compressed_fasta, _ref_seq_name)
    else:
        input_stamp, cluster_count = cluster_to_fasta (cluster_json, compressed_fasta) # changes the json to fasta, also returns the count len().
    #end if
    
    #print("# Current number of sequences", cluster_count)
    if cluster_count <= max_toRetain:
        #shutil.copy (msa_SA, msa)
        break
    else:
       step = threshold * 0.25 
       threshold += step
    #end if
#end while

sys.exit(0)
# End of file

# /usr/local/bin/tn93-cluster -f -o /home/aglucaci/SARS-CoV-2/clades/B-1-1-519/3C.query.json -t 0.0001 /home/aglucaci/SARS-CoV-2/clades/B-1-1-519/3C.query.msa.bk
