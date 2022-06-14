"""
Combine analysis results for Omicron selection

Authors:
    Sergei L Kosakovsky Pond (spond@temple.edu)
    Alexander G Lucaci (alexander.lucaci@temple.edu)

Version:
    v0.0.1 (2021-01-17)
    v0.0.2 (2022-06-11)

"""
# Imports -------------------------------------------------------------
import argparse
import csv
import random
import os
import json
import sys
import datetime, math
from collections import Counter
from operator import itemgetter

# CL Parse ------------------------------------------------------------
arguments = argparse.ArgumentParser(description='Combine analysis results for Omicron selection')
arguments.add_argument('-i', '--input',  help = 'Directory to scan', required = True, type = str)
arguments.add_argument('-o', '--output',  help = 'Save JSON here', required = True, type = str)
settings = arguments.parse_args()

# Lambda functions ----------------------------------------------------
msa_name = lambda date: date.strftime ("%m%d%Y")
file_name = lambda b,f,e: os.path.join (b,f + e)

# Declares ------------------------------------------------------------
required_extensions = {
    "variants"  : ".S.VariantCounts.json",
    "fel"       : ".S.query.msa.SA.uniq2.fas.FEL.json",
    "meme"      : ".S.query.msa.SA.uniq2.fas.MEME.json",
    "busted"    : ".S.query.msa.SA.uniq2.fas.BUSTEDS.json",
    "slac"      : ".S.query.msa.SA.uniq2.fas.SLAC.json",
    "bgm"       : ".S.query.msa.SA.uniq2.fas.BGM.json",
    "Cluster 1" : ".S.cluster1.json",
    "Cluster 2" : ".S.cluster2.json",
    "Cluster 3" : ".S.cluster3.json"
}

max_date = datetime.datetime (2000,1,1)
combined_data                 = {}
max_p                         = 0.05
nucs = set (['A','C','G','T'])
ever_selected = set ()
site_info     = []
dates         = []
substitutions = None
clusters   = {}

# Helper functions ----------------------------------------------------
def convert_d (d, key):
    result = []
    for i in range (len (d)):
        result.append ([d[str(i)][key],d[str(i)]['proportion']])
    #end for
    return result
#end method

# Main subroutine -----------------------------------------------------
print("# Exploring input directory:", settings.input)
print("# Saving results to directory:", settings.output)

#for root, dirs, files in os.walk(settings.input):
#    print(root, dirs, files)
#sys.exit(0)

for root, dirs, files in os.walk(settings.input):
    dir_name = os.path.basename(root) 
    #print("# Analysis date:", [dir_name.split("-")[1]])
    #dir_date  = datetime.datetime.strptime (str(dir_name.split("-")[1]), "%m%d%Y")
    #print("# Directory date:", dir_date)

    try:
        #dir_date  = datetime.datetime.strptime (dir_name, "%Y-%m-%d")
        dir_date  = datetime.datetime.strptime (dir_name.split("-")[1], "%m%d%Y")
        print("# Directory date:", dir_date)

        base_name = msa_name(dir_date)
        print("# basename:", base_name)

        required_names = set ([ settings.input + "/lineage_BA_2.fasta" + re for re in required_extensions.values()])
        print("# Required names:", required_names)   

        #matched_names = [k for k in files if k in required_names]
        matched_names = [k for k in required_names]
        print()
        print("# Matched files:", matched_names)

        if len (matched_names) < len (required_names):
            continue
        #end if

        is_max = False
        print()
        print("# Checking max date:", max_date) 
 
        if dir_date > max_date:
            print("# This is the latest analysis:", dir_date)
            max_date = dir_date
            is_max = True
        #end if
            
        directory_record = {}
        N_rich = set ()
        problematic = set ()
        
        #with open (file_name (root, base_name, required_extensions["variants"]), "r") as fh:
        variants_file = os.path.join(settings.input, "lineage_BA_2.fasta" + required_extensions["variants"])
        print()
        print("# Loading variants_file:", variants_file)

        with open (variants_file ,"r") as fh:
            #print (file_name (root, base_name, required_extensions["variants"]), file = sys.stderr)
            print(variants_file, file = sys.stderr)
            variant_json = json.load (fh)
            directory_record["N"] = variant_json["N"]
            directory_record["H"] = variant_json["H"]
            directory_record["HA"] = variant_json ["HA"]
            
            ## describe haplotype data; identify variants that are 
            # gappy or have too many Ns
             
            allN = 0
            for site, variants in enumerate (variant_json ["haplotypes"]):
                
                Ns   = 0
                gaps = 0
                fs   = 0
                amb  = 0
                
                for codon, ccount in variants.items():
                    if codon == 'NNN': 
                        Ns += ccount
                    elif codon == '---':
                        gaps += ccount
                    else:
                        codon_letters = set ([k for k in codon])
                        d = codon_letters.difference (nucs)
                        if len (d) > 0:
                            if '-' in d:
                                fs += ccount
                            else:
                                amb += ccount
                            #end if
                    #end if
                #end for

                if Ns * 10 > variant_json["H"] or amb * 20 > variant_json["H"]:
                    N_rich.add (site + 1)
                #end if

                if gaps * 20 > variant_json["H"] and gaps * 3 < variant_json["H"]:
                    problematic.add (site + 1)
                #end if

                if gaps * 2 > variant_json["H"]:
                    problematic.add (site + 1)
                #end if

                if fs * 100 > variant_json["H"]:
                    problematic.add (site + 1)
                #end if

                allN += Ns
            
            directory_record ["NNN"] = allN / variant_json["H"] 
        #end with        
       
        # Process SLAC  ---------------------------------------------------------------------------
        slac_file = os.path.join(settings.input, "lineage_BA_2.fasta" + required_extensions["slac"])
        print()
        print("# Processing SLAC file:", slac_file)

        #with open (file_name (root, base_name, required_extensions["slac"]), "r") as fh:
        with open (slac_file, "r") as fh:
            print (file_name (root, base_name, required_extensions["slac"]), file = sys.stderr)
            slac_info = json.load (fh)
            L    = len (slac_info["MLE"]["content"]["0"]["by-site"]["RESOLVED"])
            S    = sum ([k[2] for k in slac_info["MLE"]["content"]["0"]["by-site"]["RESOLVED"]]) / L
            S2   = sum ([k[2]*k[2] for k in slac_info["MLE"]["content"]["0"]["by-site"]["RESOLVED"]])
            NS   = sum ([k[3] for k in slac_info["MLE"]["content"]["0"]["by-site"]["RESOLVED"]]) / L
            NS2   = sum ([k[3]*k[3] for k in slac_info["MLE"]["content"]["0"]["by-site"]["RESOLVED"]])
            
            cutoff_S = S + math.sqrt (S2/L - S*S)*2
            cutoff_NS = NS + math.sqrt (NS2/L - NS*NS)*2
            
            outliers = set ([i+1 for i, k in enumerate (slac_info["MLE"]["content"]["0"]["by-site"]["RESOLVED"]) if k[2]>=cutoff_S and k[3] >= cutoff_NS])
            for site, slac_data in enumerate (slac_info["MLE"]["content"]["0"]["by-site"]["RESOLVED"]):
                if len (site_info) <= site:
                    site_info.append ({'fel' : [], 'meme' : [], 'slac' : []})
                
                site_info[site]['slac'].append (
                    [0,0]
                )
                
            for branch, slac_data in slac_info["branch attributes"]["0"].items():
                if branch[0:4] == "Node":
                    for site,k in enumerate (slac_data["nonsynonymous substitution count"][0]):
                        site_info[site]['slac'][-1][1] += k
                        site_info[site]['slac'][-1][0] += slac_data["synonymous substitution count"][0][site]
                    #end for
                #end if
            #end for
                
            problematic.update (outliers)
        
        directory_record ['issues'] = {'problematic' : sorted (list(problematic)), 'N rich' : sorted (list (N_rich))}

        # Process MEME ----------------------------------------------------------------------------
        meme_file = os.path.join(settings.input, "lineage_BA_2.fasta" + required_extensions["meme"])
        print()
        print("# Processing MEME file:", meme_file)
        
        with open (meme_file, "r") as fh:
        #with open (file_name (root, base_name, required_extensions["meme"]), "r") as fh:
            #print (file_name (root, base_name, required_extensions["meme"]), file = sys.stderr)
            print(meme_file, file=sys.stderr)
            meme_info = json.load (fh)
            directory_record ['omega'] = {
                    "leaves" : meme_info["fits"]["Global MG94xREV"]["Rate Distributions"]["non-synonymous/synonymous rate ratio for *background*"][0][0],
                    "internal" : meme_info["fits"]["Global MG94xREV"]["Rate Distributions"]["non-synonymous/synonymous rate ratio for *test*"][0][0]
            }
            for site, meme_data in enumerate (meme_info["MLE"]["content"]["0"]):
                    
                site_info[site]['meme'].append (
                    [meme_data[6],meme_data[7],meme_data[3], meme_data[4]]
                    #pv, branches, omega+, p+
                )
                
                if meme_data[6] <= max_p:
                    ever_selected.add (site)
                #end if
            #end for
        #end with

        # Process FEL -----------------------------------------------------------------------------
        fel_file = os.path.join(settings.input, "lineage_BA_2.fasta" + required_extensions["fel"])
        print()
        print("# Processing FEL file:", fel_file)

        with open (fel_file, "r") as fh:
        #with open (file_name (root, base_name, required_extensions["fel"]), "r") as fh:
            print (fel_file, file = sys.stderr)
            #print (file_name (root, base_name, required_extensions["fel"]), file = sys.stderr)
            fel_info = json.load (fh)
            for site, fel_data in enumerate (fel_info["MLE"]["content"]["0"]):
                     
                site_info[site]['fel'].append (
                    [fel_data[4],fel_data[0],fel_data[1],fel_data[5],fel_data[6],fel_data[7],fel_data[8]]
                )
                #pv, alpha, beta, branch length, omega, LB, UB
                
                if fel_data[4] <= max_p:
                    ever_selected.add (site)
                #end if
            #end for
        #end with

        
        if is_max:
            print()
            # Process slac subs file --------------------------------------------------------------
            try:
                # lineage_BA_2.fasta.S.query.msa.SA.uniq2.fas.subs.json               
                subs_file = os.path.join(settings.input, "lineage_BA_2.fasta.S.query.msa.SA.uniq2.fas.subs.json")
                print("# Openings SLAC Mapper subs:", subs_file)
                with open (subs_file, "r") as fh:
                #with open (file_name (root, base_name, ".subs.json"), "r") as fh:
                    print (subs_file, file = sys.stderr)
                    print("# Parsing subs")
                    substitutions = json.load (fh)
                    print("# Parsing complete")
                #end with
            except Exception as e:
                print("# Error in subs file")
                print (e, file = sys.stderr)
            #end try

            print("# Moving on to cluster files")
 
            # Process clusters -------------------------------------------------------------------- 
            # lineage_BA_2.fasta.S.cluster1.json 
            clusters = {}
            print("# Parsing cluster files")
            for e in ["Cluster 1","Cluster 2","Cluster 3"]:
                print("# Creating cluster file variable")
                    
                cluster_file = os.path.join(settings.input, "lineage_BA_2.fasta" + required_extensions[e])


                print("# Opening cluster file:", cluster_file)
                with open (cluster_file) as fh:
                #with open (file_name (root, base_name, required_extensions[e])) as fh:
                   #print (file_name (root, base_name, required_extensions[e]), file = sys.stderr)
                   print (cluster_file, file = sys.stderr)
                   clusters[e] = json.load (fh)
                #end with
            #end for
            # Process BGM file --------------------------------------------------------------------       
 
            bgm_file = os.path.join(settings.input, "lineage_BA_2.fasta" + required_extensions["bgm"])
            pairs = []
            if os.path.exists(bgm_file) and os.path.getsize(bgm_file) > 0:
                print("# Opening BGM json:", bgm_file) 
                with open (bgm_file, "r") as fh:
                #with open (file_name (root, base_name, required_extensions["bgm"]), "r") as fh:
                    print (bgm_file, file = sys.stderr)
                    #print (file_name (root, base_name, required_extensions["bgm"]), file = sys.stderr)
                    bgm_info = json.load (fh)
                    #pairs    = []
                    if "MLE" in bgm_info:
                        for pair_info in bgm_info["MLE"]["content"]:
                            if pair_info [4] >= 0.8:
                                pairs.append (pair_info)
                            #end if
                        #end foor
                    #end if
                #end with
            #end if

            # Process BUSTEDS file ----------------------------------------------------------------
            busted_file = os.path.join(settings.input, "lineage_BA_2.fasta" + required_extensions["busted"])

            with open (busted_file, "r") as fh:
            #with open (file_name (root, base_name, required_extensions["busted"]), "r") as fh:
                #print (file_name (root, base_name, required_extensions["busted"]), file = sys.stderr)
                print (busted_file, file = sys.stderr)
                busted_info = json.load (fh)
                #print (busted_info["fits"]["Unconstrained model"]["Rate Distributions"])
                busted = {
                    'p' : busted_info ["test results"]["p-value"],
                    'leaves' : convert_d (busted_info["fits"]["Unconstrained model"]["Rate Distributions"]["Background"], "omega"),
                    'internal' : convert_d (busted_info["fits"]["Unconstrained model"]["Rate Distributions"]["Test"], "omega"),
                    'srv' : convert_d (busted_info["fits"]["Unconstrained model"]["Rate Distributions"]["Synonymous site-to-site rates"], "rate"),
                }

            #end with
            counts = variant_json ["counts"]
            haplos = variant_json ["haplotypes"]
            haplos_all = variant_json ["all-haplotypes"]
            seq_dates = variant_json ["dates"]
        #end if 

        dates.append (base_name)
        combined_data [base_name] = directory_record
                 
    except Exception as e:
        print (dir_name, e, file = sys.stderr)
        pass
    #end try
#end for

             
sorted_dates = sorted ([[k,i] for i,k in enumerate (dates)], key = itemgetter (0))

combined_data ['analysis dates'] = [k[0] for k in sorted_dates]

#print("# BUSTEDS:", busted)
#combined_data ['busted'] = busted

site_selection_info = {}
            
for k in sorted (list (ever_selected)):
   site_selection_info [k+1] = {
        'fel' :  [site_info[k]['fel'][d[1]] for d in sorted_dates],
        'meme' : [site_info[k]['meme'][d[1]] for d in sorted_dates],
        'slac' : [site_info[k]['slac'][d[1]] for d in sorted_dates],
   }
                
combined_data ['p'] = max_p
combined_data ['sites'] = site_selection_info

combined_data ['pairs'] = pairs
combined_data ['counts'] = counts
combined_data ['haplos'] = haplos
combined_data ['haplos_all'] = haplos_all
combined_data ['dates'] = seq_dates
combined_data ['subs'] = substitutions
combined_data ['clusters'] = clusters
  
with open (settings.output, 'w') as fh:
    json.dump (combined_data, fh)

# End of file

