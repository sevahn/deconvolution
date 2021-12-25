"""
Generate individual slurm jobs for each sample to be deconvolved!
"""


import pandas as pd
import numpy as np

########## ENTER ARGUMENTS HERE ###############

# path to your CPM-normalized counts
# IMPORTANT: first row should be gene names, column names should be samples
sampPath = "samples.csv"

env = "snakemake" # name of the conda environment 

partition = "owners,normal" # partition on supercomputer for jobs to be run 

### other params to tweak (need be) 
memory = "10000" # job memory
time = "300:00" # job run time (for jobs with large nu values they take longer) 
scriptName = "deconvolve" # script name is deconvolve.py, change if different for you
cpmThresh = 1 # minimum gene expression level in sample (can increase/decrease need be)
deconMethod = "nuSVR" # specify to "nuSVR" or "nnls" or "QP" (whichever kind of deconvolution you'd like to perform)

########## DO NOT EDIT PAST THIS LINE ###########
importName = "deconv"
pre = "decon"
hyper = "cfRNA_deconv"
jackknife = "False" 

# read in list of samples, where the columns are the sample names
samps = list(pd.read_csv(sampPath, index_col = 0).columns)

biologRepSRR = {}
for i in samps:
    biologRepSRR[i] = [i]

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
    file.write("#SBATCH --partition=" + partition + "\n")
    file.write("#SBATCH --cpus-per-task=1\n")
    file.write("#SBATCH --mem=" + str(memory) + "\n")
    file.write("#SBATCH --time=" + time + "\n")
    file.write("\n\n\n")
    file.write("source activate " + env + "\n")
    
    file.write("\n")

    # line to call nu-SVR deconvolution 
    file.write('python3 -c \"import ' + scriptName + ' as ' + importName + '; ' + importName + ".main(" + str(cpmThresh) + ",  \'" + sampPath + '\', ' + str(srr) + ' , \'' + 'nuSVR' + '\', \'' + bioRep + '\',  jackknife = ' + jackknife + ')\"')
    file.close()
