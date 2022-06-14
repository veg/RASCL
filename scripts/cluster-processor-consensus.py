
# Imports
import sys, json
from Bio import SeqIO
from   operator import itemgetter, attrgetter
from collections import Counter

# Declares
master     = {}
attributes = {}
sep_char = '/'

SL = 0
nucs = set (["A","C","G","T"])

#cluster_min = snakemake.params.cluster_min
#input_msa = snakemake.params.input_msa
#input_uniq_all = snakemake.params.input_uniq_all
#input_t1_clusters = snakemake.params.input_t1_clusters
#input_t0_clusters = snakemake.params.input_t0_clusters
#input_u_clusters = snakemake.params.input_u_clusters
#output = snakemake.params.output



input_msa         = sys.argv[1]
input_uniq_all    = sys.argv[2]
input_t1_clusters = sys.argv[3]
input_t0_clusters = sys.argv[4]
input_u_clusters  = sys.argv[5]
#cluster_min       = sys.argv[6]
output            = sys.argv[6]



# Helper functions
def seq_id (id):
    return id.split (sep_char)[2]
#end method
    
def cluster_label (core, ids):
    md  = None
    mxd = None
    countries = Counter ()
    for sid in ids:
        if len(sid):
            sd = attributes[sid]
            if sd[1] and len (sd[1]) > 7:
                if md:
                    md = min (md, sd[1])
                else:
                    md = sd[1]
                if mxd:
                    mxd = max (md, sd[1])
                else:
                    mxd = sd[1]
            countries[sd[0]] += 1
        
    countries = sorted (countries.items(), key = itemgetter (1), reverse = True )    
    
    id = "cluster_%s/%s-%s/%s||%d" % (core, 
                                     str (md), 
                                     str (mxd), 
                                     "%d_%s" % (len (countries), "_".join ([k [0] for k in countries[:3]])), 
                                     len(ids))
    return id
#end method
  
def do_sequence (seq, by_position):
    for k in range (SL):
        by_position[k][seq[k*3:k*3+3]] += 1
    #end for
#end method
    
def compute_consensus (ids):
    print ("Generating the consensus of %d sequences " % len (ids), file = sys.stderr)
    by_position = [Counter () for k in range (SL)]
    for sid in ids:
        if len (sid) > 0:
            do_sequence (master[sid], by_position)
            #print (by_position)    
    consensus = []
    for k in by_position:
        all = [[n,c] for n,c in k.items()]
        resolved = [n for n in all if len ([nn for nn in n[0] if nn not in nucs]) == 0]
        if len (resolved) >= 1:
            max_count = max (resolved, key = itemgetter (1))
            if max_count [1] > 1:
                consensus.append (max (resolved, key = itemgetter (1))[0])
                continue
        consensus.append (max (all, key = itemgetter (1))[0])
            
    return "".join (consensus)    
#end method

# Main subroutine -----------------------------------------------------
with open(input_msa, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        S = str (record.seq)
        SL = len (S) // 3
        try:
            tags_pipe = record.id.split ('|')
            tags = record.id.split ('/')
            #hCoV-19/USA/WA-CDC-UW22013065303/2022|2022-01-30|2022-02-12||228615
            epi_id = seq_id(record.id)
            attributes[epi_id] = [tags[1],tags_pipe[2]]
            master[epi_id] = S
        except Exception as e:
            #print (record.id, file = sys.stderr)
            pass
        #end try
    #end for
#end with

for_analysis = {}
counts = {}

with open(input_uniq_all, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        try:
            for_analysis [seq_id (record.id)] = []
            counts [seq_id (record.id)] = record.id.split ('|')[-1]
        except Exception as e:
            #print (tags_slash)
            pass
        #end try
    #end for
#end with
            
#print (counts, file = sys.stderr)
                        
#for f in sys.argv[4:]: 
for f in [input_t1_clusters, input_t0_clusters, input_u_clusters]:              
    with open(f) as handle:
        cluster_info = json.load (handle)
        ci = {}
        N = 0
        for i, c in enumerate (cluster_info):
            for m in c["members"]:
                try:
                    N += int (m.split ('||')[1])
                except:
                    N += 1
                try:
                    ci[seq_id(m)] = i
                except Exception as e:
                    pass
        #end inner for

        for seq, members in for_analysis.items():
            check_list = set (members)
            check_list.add (seq)
            others = set ()
            for s in check_list:
                if s in ci:
                    for m in cluster_info[ci[s]]['members']:
                        try:
                            others.add (seq_id(m))
                        except:
                            print ("Malformed %s" % m, file = sys.stderr)
        
            for_analysis[seq] = list (check_list.union (others))
        #end inner for
    #end with            
#end for 
#S = 0 

#minL = int (sys.argv[1])
#minL = int(cluster_min)
minL = int(3)
# group by consensus because some of the sequences may compress to the same consensus
by_consensus = {}


for seq, members in for_analysis.items():
   if len (members) >= minL:
        #print (seq, len (members))
        #S += len (members) 
        #if (len (members) < 10):
        consensus = compute_consensus (members)
        if consensus in by_consensus:
            by_consensus[consensus]["members"].extend (members)
        else:
            by_consensus[consensus] = {'seq' : seq, 'members' : members}
        #end inner if
        #print (">%s\n%s\n" % (cluster_label (seq, members), compute_consensus (members)))
    #end if
#end for  
#print (S)       
 
for consensus, members in by_consensus.items():
    with open(output, "a") as fh:
        print (">%s\n%s\n" % (cluster_label (members['seq'], members["members"]), consensus), file=fh)
    #end with
    fh.close()
#end for       

 
            