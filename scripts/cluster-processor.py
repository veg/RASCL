import sys, json

blacklist = set ()

input_file = snakemake.params.input

input_file2 = ""
input_file2 = snakemake.params.input2
output_file = snakemake.params.output
output_file2 = snakemake.params.output2


#if len(sys.argv) > 2:
if input_file2 != "":
    with open (input_file2, "r") as fh:
        cluster_json = json.load (fh)
        largest = max ([len(k['members']) for k in cluster_json])
        for k in cluster_json:
            if k["size"] != largest:
                for n in k["members"]:
                    blacklist.add (n)
                #end for
            #end if
        #end for
    #end with   
#end if
   

#with open (sys.argv[1], "r") as fh:
with open (input_file, "r") as fh:
    cluster_json = json.load (fh)
    for k in cluster_json:
        N = 0
        for m in k['members']:
            s = m.split ("||")
            if (len (s) > 1):
                N += int (s[1])
            else:
                N += 1
            #end if
        #end for
        seq = k['centroid'].split ("\n")
        if seq[0].split ('>')[1] in blacklist:
            print ("%s" % "\n".join (seq), file = sys.stderr)
        else:
            id = seq[0].split ("||")[0] + "||" + str (N)
            #print (N, file = sys.stderr)
            with open(output_file, "a") as fh:
                print (id, file=fh)
                print (seq[1].strip(), file=fh)
            #end with
            fh.close()
        #end if
        #NN += N
    #end for
#end with
    
