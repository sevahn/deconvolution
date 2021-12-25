import pandas as pd
import numpy as np
import scanpy as sc

import os

# merge together all the csv files, one for support vectors and the other for coefficients
sampData = "/oak/stanford/groups/quake/sevahn/alzheimers/ad_cpmOnly_postQC_unstranded_FINAL.csv"
#sampData = "/oak/stanford/groups/quake/sevahn/molecstetho/all_molec_stetho_dat/analysis/molecstetho_htseq_merged.csv" 
samps = pd.read_csv(sampData, sep = ",", index_col = [0,1])
sampNames = samps.columns
basePath = os.getcwd() + "/"
print(basePath)


supVecBase = "_coarsegrain___supportVectors.csv"
coefBase = "_deconvolutionCoefs.csv"

allSuppVec = pd.DataFrame()
allCoefs = pd.DataFrame()

for bioRep in sampNames:
    print(bioRep)
    if bioRep in ['2520', '2499','2503']: continue 
    print(basePath + bioRep + "/" + bioRep +  supVecBase)
    print(basePath + bioRep + "/" + bioRep + coefBase)
    print("")
    sampSV = pd.read_csv(basePath + bioRep + "/" + bioRep +  supVecBase, sep = ",", index_col = 0)
    sampCoef = pd.read_csv(basePath + bioRep + "/" + bioRep + coefBase, sep = ",", index_col = 0)
    allCoefs = pd.concat([sampCoef, allCoefs], ignore_index = False)
    allSuppVec = pd.concat([sampSV, allSuppVec], axis = "columns", ignore_index = False)

#allCoefs.index = sampNames

fprefix = "allAD_unstranded_truseq3_sigmatv18.0_"
fend = "_12142021.csv"

allSuppVec.to_csv(fprefix + "PLASMA_allSV" + fend, sep = ",", header = True, index = False)
allCoefs.to_csv(fprefix + "PLASMA_allCoefs" + fend, sep = ",", header = True, index = True)

