# Sevahn Vorperian
# Quake Lab, Stanford University
# Code to generate per-organ compartment dendrograms on Tabula Sapiens
# July 30 2021

import numpy as np
import scanpy as sc
import pandas as pd
from random import sample
from scipy.cluster.hierarchy import dendrogram, linkage

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

import os

datPath = "/oak/stanford/groups/quake/sevahn/tsp/trainvaltest/" 
adataName = "TSP1_TSP15_scvi_donor-method_normalized-log1p-scaled_annotated_clean.h5ad"
cellAnnotCol = "cell_ontology_class" # col of adata.obs with the cell annotations
basePath = "/oak/stanford/groups/quake/sevahn/tsp/trainvaltest/"
fend = "_cleanTSP_07312021"

# get the disassociation genes from Supp Table 5 and exclude from subsequent analysis
dissassociationGenes = pd.read_excel(datPath + "SuppTab5_GenesAffectedByDissociation.xlsx", dtype = str)
dissassociationGenes = dissassociationGenes["Gene_Symbol"].values.tolist()

##### DO NOT EDIT BELOW THIS LINE #####
adatPath = datPath + adataName
adata = sc.read_h5ad(adatPath)
os.system("echo got adata")

adata.X = adata.layers["decontXcounts"]
sc.pp.normalize_total(adata, target_sum = 1e6)
sc.pp.log1p(adata)
os.system("echo normalized and logged :\) ")

# drop the disassociation genes
keepGenes = np.setdiff1d(adata.var_names, dissassociationGenes)
adata = adata[:, keepGenes]

# only do the 10X data
adata = adata[adata.obs["method"] == "10X"]

# drop the eyes
adata = adata[adata.obs["tissue"] != "Eye"]

# collapse the cells
tspCellTypes = list(np.unique(adata.obs[cellAnnotCol])) 

# handle extra muscle annots
allMuscle = [i for i in tspCellTypes if "muscle" in i] + ['myometrial cell']
goodMuscle = ["cardiac muscle cell", "cell of skeletal muscle", "smooth muscle cell"]
dropMuscle = np.setdiff1d(allMuscle, goodMuscle).tolist()

# handle fibroblasts
allFibro = [i for i in tspCellTypes if 'fibro' in i]
dropFibro = np.setdiff1d(allFibro, ['fibroblast']).tolist()

# 'corneal keratocyte' is a type of keratocyte 
# 'connective tissue cell' is too broad
# 'lacrimal gland functional unit cell' is several cell types
# 'ciliary body' is also a collection of multiple cell types 
annotExclude = ["epithelial cell", "ocular surface cell", 'radial glial cell', 'lacrimal gland functional unit cell', 'connective tissue cell', 'corneal keratocyte', 'ciliary body'] 
annotExclude += dropMuscle
annotExclude += dropFibro
goodAnnot = np.setdiff1d(adata.obs[cellAnnotCol].values.tolist(), annotExclude)

# subset the adata to the good annotations
adata = adata[adata.obs[cellAnnotCol].isin(goodAnnot)]
os.system("echo " + str(adata.obs.shape))

# RUN THE PCA ON THE adata object with the correct transformed gene counts,
# excluded diassociation genes, and the correct annotations
sc.pp.pca(adata)
os.system("echo got PCA")

# total number of unique cell types
tspCellTypes, uniqCell = np.unique(adata.obs[cellAnnotCol], return_counts = True)

##### BROADEN LABELS FOR CERTAIN CELL TYPES ##### 
cgColName = "coarsegrain_cell"

# compute the dendrogram
zKey = "dendrogram_" + cgColName 

compartments = np.unique(adata.obs["compartment"].dropna())
print("COMPARTMENTS = ", compartments)

for compartment in compartments:
    if compartment == "germ line": continue # "sperm" is the only annotation here
     
    adataSubset = adata[adata.obs["compartment"] == compartment]
    cellsInComp = np.unique(adataSubset.obs["cell_ontology_class"]).tolist()

    if compartment == "epithelial":
        cellsInComp += ['hillock-club cell of prostate epithelium']
        cellsInComp.remove("mesothelial cell") # this is stromal only

    # do off the joined cell types, not the compartment, since some cell types have mis-annotated compartments
    adataSubset = adata[adata.obs["cell_ontology_class"].isin(cellsInComp)]
    
    saveEnd = compartment + "-" + fend 
    sc.tl.dendrogram(adataSubset, n_pcs = 50, linkage_method = "complete", key_added = zKey, optimal_ordering = False, groupby = cellAnnotCol)
    plt.rcParams['figure.figsize'] = [20,20]


    # save the linkage df
    Z = adataSubset.uns[zKey]["linkage"]
    zDF = pd.DataFrame(Z)
    zDF.to_csv("tsp_linkageDF_" + saveEnd + ".csv")

    ivlList =  adataSubset.uns[zKey]["dendrogram_info"]['ivl']
    ivlDF = pd.DataFrame(index = ivlList)
    ivlDF.to_csv("tsp_ivl_" + saveEnd + ".csv")

    # save the heights of the nodes
    dcoordMat = adataSubset.uns[zKey]["dendrogram_info"]["dcoord"]
    dcoordDF = pd.DataFrame(data = dcoordMat)
    dcoordDF.to_csv("tsp_dcoord_" + saveEnd + ".csv")

    leafList = adataSubset.uns[zKey]["dendrogram_info"]['leaves']
    leafDF = pd.DataFrame(index = leafList)
    leafDF.to_csv("tsp_leafNames_" + saveEnd + ".csv")

    sc.pl.dendrogram(adataSubset, dendrogram_key = zKey, groupby = cellAnnotCol, save = "tsp1p2_PCA_" + saveEnd + ".png")



