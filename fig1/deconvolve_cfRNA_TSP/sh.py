"""
Sevahn Vorperian
18 July 2020
Iteratively write the 12 thuy samples for deconvolution with sevDecon
"""
import pandas as pd
import numpy as np

########## ENTER ARGUMENTS HERE ###############
memory = "10000"
time = "300:00" # make 48:00:00 for 2 days if jackknifing
scriptName = "deconvolve"
importName = "sev"
pre = "decon"
hyper = "-10Xscaled-allNCI"
env = "snakemake" # name of the conda environment
jackknife = "False"
deconMethod = "nuSVR" # "nuSVR" or "nnls" or "QP"
cpmThresh = 1

########## NO NEED TO EDIT PAST THIS LINE ###########
#thuy = "/home/groups/quake/sevahn/deconvolution/sevDecon/NPC_ThuyAnalysis.csv"
#thuy = "/oak/stanford/groups/quake/sevahn/molecstetho/all_molec_stetho_dat/analysis/molecstetho_htseq_merged.csv" 
#samps = pd.read_csv(thuy, index_col = (0, 1), sep = ",").columns
#thuy = "/home/groups/quake/sevahn/deconvolution/sevDecon/NPC_ThuyAnalysis.csv"
#thuy = "/oak/stanford/groups/quake/sevahn/molecstetho/all_molec_stetho_dat/analysis/molecstetho_htseq_merged.csv" 
thuy = "/oak/stanford/groups/quake/sevahn/alzheimers/ad_cpmOnly_postQC_unstranded_FINAL.csv"

samps = list(pd.read_csv(thuy, index_col = 0).columns)
biologRepSRR = {}
for i in samps:
    biologRepSRR[i] = [i]

#for srr in buffy + plasma: 
for bioRep in biologRepSRR: 
    srr = biologRepSRR[bioRep]
    print(srr)
    shname = bioRep + hyper + ".sh"
    file = open(shname, "w+")
    file.write("#!/bin/bash")
    file.write("\n")
    file.write("#SBATCH --job-name=" + pre + bioRep + hyper + "\n")
    file.write("#SBATCH --output=" + pre + bioRep + hyper + ".out" + "\n")
    file.write("#SBATCH --error=" + pre + bioRep + hyper + ".err" + "\n")
    file.write("#SBATCH --qos=normal\n")
    file.write("#SBATCH --ntasks=1\n")
    file.write("#SBATCH --partition=quake,normal\n")
    file.write("#SBATCH --cpus-per-task=1\n")
    file.write("#SBATCH --mem=" + str(memory) + "\n")
    file.write("#SBATCH --time=" + time + "\n")
    #file.write("#SBATCH --mail-type=ALL\n")
    #file.write("#SBATCH --mail-user=sevahn@stanford.edu\n")
    file.write("\n\n\n")
    file.write("source activate " + env + "\n")
    
    file.write("\n")
    file.write('python3 -c \"import ' + scriptName + " as " + importName + "; " + importName + ".main(" + str(cpmThresh) + ", \'" + thuy + "\', " + str(srr) + " , \'" + "NNLS" + "\', \'" + bioRep + "\',  jackknife = " + jackknife + ')\"')
    file.write("\n")
    file.write('python3 -c \"import ' + scriptName + ' as ' + importName + '; ' + importName + ".main(" + str(cpmThresh) + ", \'" + thuy + '\', ' + str(srr) + ' , \'' + 'QP' + '\', \'' + bioRep + '\',  jackknife = ' + jackknife + ')\"')
    file.write('\n')
    file.write('python3 -c \"import ' + scriptName + ' as ' + importName + '; ' + importName + ".main(" + str(cpmThresh) + ",  \'" + thuy + '\', ' + str(srr) + ' , \'' + 'nuSVR' + '\', \'' + bioRep + '\',  jackknife = ' + jackknife + ')\"')
    file.close()
