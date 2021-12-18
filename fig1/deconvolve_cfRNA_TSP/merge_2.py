# merge together all the csv files, one for support vectors and the other for coefficients
# this script goes into each of the deconvolved sample directories
# and merges the coefficient files and the support vector files
# the best coefficient will then be chosen locally 

import pandas as pd
import numpy as np
import scanpy as sc

import os
# dataset name
datasetName = "NCI_cfRNA" # name of dataset
fend = "_12142021.csv" # name of file ending...
# path to your data
sampData = "/oak/stanford/groups/quake/sevahn/alzheimers/ad_cpmOnly_postQC_unstranded_FINAL.csv"
sampNames = samps.columns
basePath = os.getcwd() + "/"


#### do not edit below this line

# read in samples
samps = pd.read_csv(sampData, sep = ",", index_col = [0,1])

supVecBase = "_coarsegrain___supportVectors.csv"
coefBase = "_deconvolutionCoefs.csv"

allSuppVec = pd.DataFrame()
allCoefs = pd.DataFrame()

for bioRep in sampNames:
    sampSV = pd.read_csv(basePath + bioRep + "/" + bioRep +  supVecBase, sep = ",", index_col = 0)
    sampCoef = pd.read_csv(basePath + bioRep + "/" + bioRep + coefBase, sep = ",", index_col = 0)
    allCoefs = pd.concat([sampCoef, allCoefs], ignore_index = False)
    allSuppVec = pd.concat([sampSV, allSuppVec], axis = "columns", ignore_index = False)



allSuppVec.to_csv(datasetName  + "_PLASMA_allSV" + fend + ".csv", sep = ",", header = True, index = False)
allCoefs.to_csv(datasetName + "_PLASMA_allCoefs" + fend + ".csv", sep = ",", header = True, index = True)

