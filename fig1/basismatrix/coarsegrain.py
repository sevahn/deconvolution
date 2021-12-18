"""
Sevahn Vorperian, Quake Lab
27 July 2021
sample up to 30 of each cell type for basis matrix generation with CIBERSORTx
reaassign cell type labels based on cut tree
"""

import numpy as np
import scanpy as sc
import pandas as pd
from random import sample

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

import os

basePath = "/oak/stanford/groups/quake/sevahn/tsp/trainvaltest/" 
adataName = "TSP1_TSP15_scvi_donor-method_normalized-log1p-scaled_annotated_clean.h5ad"

cutTreeDF = "joinedLabels_cutPerCompartment_cleanTSP_sigMat18_08012021.csv"

goodAnnot = "cell_ontology_class"
ss2or10x = "10X" # ['10X', 'smartseq2']
fend = "_v18_08012021.txt"

# do not edit below this line

############# READ IN NECESSARY FILES ###############
cgDir = os.getcwd() + "/"
master = sc.read_h5ad(basePath + adataName) 

# set the counts to be the decontX counts -- no CPM normalization, the basis matrix feature in CIBERSORTx handles that
master.X = master.layers['decontXcounts']

# get the disassociation genes from Supp Table 5 and exclude from subsequent analysis
dissassociationGenes = pd.read_excel(basePath + "SuppTab5_GenesAffectedByDissociation.xlsx", dtype = str)
dissassociationGenes = dissassociationGenes["Gene_Symbol"].values.tolist()
keepGenes = np.setdiff1d(master.var_names, dissassociationGenes) 

# read in the coarsegrained labels, this is the output of cuttree.ipynb
cutTreeDF = pd.read_csv(cutTreeDF, index_col = 0)

#############################################################################
#      DOWNSAMPLE EACH CELL TYPE PRIOR TO 'COARSEGRAINING' ANNOTATION       # 
#############################################################################
adata = master[:, keepGenes]
adata = adata[adata.obs["method"] == ss2or10x]
adata = adata[adata.obs["tissue"] != "Eye"]

# eliminate the extraenous observations, like the immune cell types and the too broad epithelial and endothelial cell type annotations
keepCells = cutTreeDF.index.tolist()
adata = adata[adata.obs[goodAnnot].isin(keepCells)]
print("uniq cell ", np.unique(adata.obs["cell_ontology_class"]))

######## SAMPLE 30 OF EACH CELL TYPE ########
downsampledCells = pd.DataFrame(index = adata.var_names)
for cell in keepCells:
    obsThisCell = adata[adata.obs[goodAnnot] == cell]
    numObs = obsThisCell.shape[0]
    maxObs = 30 # get maximum 50 cells per type
    if numObs >= maxObs:
        indexSample = sample(range(0, numObs), maxObs)
    else:
        indexSample = range(0, numObs)

    # stuff to append
    data = obsThisCell.X[indexSample, :].T
    cellLabs = obsThisCell.obs.iloc[indexSample,:][goodAnnot].to_list()

    # make a dataframe
    thisCellDF = pd.DataFrame(columns = cellLabs, index = adata.var_names, data = data.todense())
    print(cell, " ", thisCellDF.shape)
    print('all = ', downsampledCells.shape) 
    if cell == "endothelial cell":
        thisCellDF.to_csv("endothelial_fuck.csv",sep = ",", index = True, header = True)
        #continue

    downsampledCells = downsampledCells.join(thisCellDF)
    
# drop the all zero rows (e.g. zero genes)
#downsampledCells = downsampledCells.loc[~(downsampledCells == 0).all(axis = 1)]

downsampledCells.to_csv("cleanTSP_30each_noLabelChange" + fend, index = True, header = True)

#############################################################################
#                     NOW PERFORM THE COARSEGRAINING                        #
#############################################################################

######## SIMPLIFY THE ANNOTATIONS ########

# convert to a dictionary, called lumpDict
print("cutTreeDF.head() ", cutTreeDF.head())
lumpDict = pd.DataFrame.to_dict(cutTreeDF, orient = "dict")["joinedLabel"]

downsampledCells.rename(columns = lumpDict, inplace = True)

# save output
downsampledCells.to_csv("cg-cleanTSP-" + ss2or10x + fend, sep = "\t")
