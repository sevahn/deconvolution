# Sevahn Vorperian
# Quake Lab @ Stanford University
# Deconvolve cfRNA into relative fractional contributions of cell type specific RNA
# December 2021

# the essentials
import numpy as np
import pandas as pd
import scanpy as sc

# scale mixture/matrix to zero mean/unit var
from sklearn import preprocessing

# for computing RMSE 
from sklearn.metrics import mean_squared_error

# for pearson correlation  
from scipy import stats

# deconvolution functions
from sklearn.svm import NuSVR # for nu-SVR
from scipy.optimize import nnls # for NNLS
from cvxopt import matrix # for QP 
from cvxopt import solvers # for QP

# for file writing and timing deconvolution
import os
import time

# PATH FOR BASIS MATRIX - adjust for wherever this is for you 
BASISMATRIX_PATH = os.getcwd() + "/" + "tsp_v1_basisMatrix.txt"

######### no need to edit below this line #########

def processMixture(mixturePath):
    """
    Read in mixtures, eliminate genes that are all zero counts in sample, 
    normalize to CPM.

    Parameters
    ----------
    mixturePath: str
        file path for mixture file
    
    Returns
    -------
    thuyCPM: pd.DataFrame (genes x mixtures)
        dataframe with CPM counts as values in matrix
    """
    thuyDat = pd.read_csv(mixturePath, sep = ",", index_col = 0)
    
    # no need to CPM normalize, data passed in is already CPM-normalized
    thuyCPM = thuyDat

    # drop the all-zero rows from thuyCPM
    thuyCPM = thuyCPM[(thuyCPM.T != 0).any()]
    return(thuyCPM)

def nuSVR(sigMat, mixtures, nuVal, cVal, droppedGene):
    """
    Perform nu-SVR on mixtures given reference signature matrix and mixtures for deconvolution.
    Code inspired from AutoGeneS.

    Parameters
    ----------
    sigMat : pd.DataFrame (genes x cellTypes)
        Reference signature matrix for deconvolution (genes x cell types) in CPM space.
        Rows are genes and cols are the different cell/tissue types
    mixtures : pd.DataFrame (genes x mixtures)
        cfRNA mixtures for deconvolution in CPM space.
        Rows are genes and cols are the different mixtures.
    nuVal: float
        nu-SVR hyperparameter, is the upper bound on margin errors and lower bound on fraction of support vectors 
    cVal: float
        nu-SVR hyperparameter, is the regularization strength 
    droppedGene: str
        name of gene to drop from basis matrix if jackknifing, commonly "" 
    Returns
    -------
    coefs : numpy array of (Mixtures x cellTypes)
        nu-SVR hyperplane coefficients >= 0 (else normalized to zero)
    svDict : dictionary
        keys = mixtures and values are lists of genes used to define the hyper place
    """
    clf = NuSVR(C = cVal, nu = nuVal, kernel = "linear")
    
    # coefs = (mixture x cell types matrix)
    # possessing the hyperplane coefficients of the sigMatrix columns 
    coefs = np.zeros((mixtures.shape[1], sigMat.shape[1]))
    svDict = {}
    fileprefix = ""
    for m in range(mixtures.shape[1]):
        fileprefix += mixtures.columns[m]
        thisMix = mixtures.iloc[:, m]
        t0 = time.time()  
        clf.fit(sigMat, thisMix)
        os.system("echo done fitting -- total fit time was " + str(time.time() - t0))
        coefficients = clf.coef_
        supportVecs = clf.support_ # mask of the indices of the support vectors
        
        # get the genes that were used in the deconvolution
        goodGenes = sigMat.index[supportVecs].to_list()
        
        sampName = mixtures.columns[m]
        os.system("echo sampName in nuSVR " + sampName)
        svDict[sampName] = goodGenes
        coefs[m] = coefficients[0]

        # make predictions on mixture
        preds = clf.predict(sigMat)
        intercept = clf.intercept_[0] * np.ones(sigMat.shape[0])
 
    fileprefix += "_coarsegrain_" + droppedGene + "_" 
    return(coefs, svDict, fileprefix, preds, intercept) 

def NNLS(sigMat, mixture, sumToOne = False):
    # call as is
    mixture = np.squeeze(mixture)
    coefs, resid = nnls(sigMat, mixture, maxiter = 10 ** 3)
    if sumToOne == True:
        # add sum to one constraint
        print("sum to one in NNLS not currently supported")
    
    # save the coefs/resid as a dataframe
    coefs = coefs.reshape((1, len(coefs)))
    return(coefs, resid) 

def scale(sigMat, mixture):
    """
    Scale basis matrix and mixture to zero mean and unit variance to improve runtime performance

    Parameters
    ---------
    sigMat: pd.DataFrame 
    Basis matrix that will be used to deconvolve samples

    mixture: pd.DataFrame 
    Mixture that will be deconvolved
 
    Returns
    --------
    scaledMix: pd.DataFrame
    Mixture scaled to zero mean and unit variance
    
    scaledMat: pd.DataFrame
    Basis matrix scaled to zero mean and unit variance

    sigMat: pd.DataFrame
    Basis matrix that's unscaled
    
    mixture: pd.DataFrame
    Unscaled mixture 
    """
    
    # scale along the columns (e.g. a given sample or cell type)
    scaledMat = preprocessing.scale(sigMat.values)
    scaledMat = pd.DataFrame(data = scaledMat, index = sigMat.index, columns = sigMat.columns) 
    
    # also scale the mixture separately to zero mean and unit variance
    scaledMix = preprocessing.scale(mixture)
    scaledMix = pd.DataFrame(data = scaledMix, index = mixture.index, columns = mixture.columns)
    return(scaledMix, scaledMat, sigMat, mixture)

def QP(sigMat, mixture):
    """
    Perform quadratic programming. Referenced guide here for implementation:
    https://courses.csail.mit.edu/6.867/wiki/images/a/a7/Qp-cvxopt.pdf

    Parameters:
    -----------
    @param sigMat = pd.DataFrame of the scaled basis matrix
    @param mixture = pd.DataFrame of the scaled mixture

    Returns:
    -----------
    coefs, a numpy array of the learned coefficient weights
    """
    # set up matrix product
    P = matrix(sigMat.values.T.dot(sigMat.values))
    q = matrix(-1 * sigMat.values.T.dot(mixture.values))
    
    # constraint of coefficients summing to one
    dimMix = mixture.shape[0]
    numCells = sigMat.shape[1] 
    A = matrix(np.ones(numCells, dtype = "float").reshape((1, numCells)))
    b = matrix(np.ones(1, dtype = "float").reshape((1,1)))
    
    # constraint of coefficients each greater than or equal to zero (nonnegativity)
    # eg product with negative identity matrix is >= 0 
    G = matrix(-1.0 * np.identity(numCells))
    h = matrix(0.0, (numCells, 1))
    
    # solve the problem
    soln = solvers.qp(P, q, G, h, A, b) 
    coefs = soln["x"]
    coefs = np.asarray(coefs)
    coefs = coefs.reshape((1, len(coefs)))
    return(coefs)

def deconvolve(scaledMixtures, scaledSigMat, deconMethod, droppedGene, biologRepName):
    """
    Perform deconvolution and return performance metrics (pearson R and RMSE)

    Parameters
    ----------
    adataPath: str
        File path for adata
    
    cpmThresh: int 
        Minimum CPM threshold for counts in mixture to be considered; values less than this will be excluded for deconvolution.
 
    deconMethod: str (either, 'qp', 'nnls', or 'nusvr')
        Deconvolution method to perform on samples 
    
    droppedGene: str
        Name of gene to be excluded from basis matrix if performing jackknifing. 
        Currently set to "" (e.g. no jackknifing) 
    
    biologRepName: str
        Name of biological replicate that will be deconvolved for labeling in the output support vector/coefficient csv files
 
    Returns
    -------
    nothing ("foo"): outputs relevant files along the way. 
    (e.g. support vectors per nu/C combo + corresponding learned coefs)
    """
    ######### DECONVOLUTION TIME!!! ##########
    srrName = biologRepName
 
    thisDir = os.getcwd()                       
    savePath = os.path.join(thisDir, srrName)
    if not os.path.exists(savePath):
        os.mkdir(savePath)
    savePath += "/"

    if droppedGene != "":
        coefFile = srrName + "_" + deconMethod.upper()  + "_Coefs_JackKnifingGenes.csv"
        srrName += "-" + droppedGene 
    else:
        coefFile = srrName + "_deconvolutionCoefs.csv" 
        srrName += "-" + deconMethod.upper()
 
    ### INITIALIZE FILE THAT WILL BE SAVED ##
    if os.path.exists(savePath + coefFile):
        aggCoefDF = pd.read_csv(savePath + coefFile, sep = ",", index_col = 0)
    else:
        os.system("creating coefficient file")
        aggCoefDF = pd.DataFrame(columns = scaledSigMat.columns.tolist() + ["r", "rmse"]) 
   
    if deconMethod.upper() == "QP":
        os.system("echo Performing Quadratic Programming")
        coefs = QP(scaledSigMat, scaledMixtures)
        coefDF = pd.DataFrame(data = coefs, index = [srrName], columns = scaledSigMat.columns)        
        groundTruth = np.asarray(scaledMixtures.values, dtype = "float").squeeze()
       
        # non-negativity & sum-to-one is implicit in QP function, use coefs as-is 
        preds = np.asarray(scaledSigMat.values.dot(coefs.T), dtype = "float").squeeze()
        pearsonR_QP, pearsonP_QP = stats.pearsonr(groundTruth, preds)
        rmse_QP = np.sqrt(mean_squared_error(preds, groundTruth))

        # concat to coefficient dataframe
        coefDF["r"] = [pearsonR_QP]
        coefDF["rmse"] = [rmse_QP]
        aggCoefDF = pd.concat([aggCoefDF, coefDF])

    if deconMethod.upper() == "NNLS":
        os.system("echo performing NNLS")
        coefs, resid = NNLS(scaledSigMat, scaledMixtures, sumToOne = False)
        if droppedGene != "": dfSampName = srrName + "-" + droppedGene
        else: dfSampName = srrName
        coefDF = pd.DataFrame(data = coefs, index = [dfSampName], columns = scaledSigMat.columns)
        
        groundTruth = np.asarray(scaledMixtures.values, dtype = "float").squeeze()
        
        normCoefs = coefDF # make a copy of coefs
        normCoefs[normCoefs < 0] = 0 # set < 0 -> 0, implicit in NNLS
        
        # enforce sum to one
        normCoefs = normCoefs / normCoefs.sum(axis = "columns").values[0]
        
        # get dot product on predictions
        preds_norm = np.asarray(scaledSigMat.values.dot(normCoefs.T), dtype = "float").squeeze()
        pearsonR_NNLS, pearsonP_NNLS = stats.pearsonr(groundTruth, preds_norm)
        coefDF["rmse"] = [np.sqrt(mean_squared_error(preds_norm, groundTruth))]
        coefDF['r'] = [pearsonR_NNLS] 
         
        aggCoefDF = pd.concat([aggCoefDF, coefDF])
        
    if deconMethod.upper() == "NUSVR":
        os.system("echo Performing nu-SVR")
        nuVals = [0.05, 0.1, 0.15, 0.25, 0.5, 0.75]
        cVals = [0.1, 0.5, 0.75, 1, 10]
        coefDF, svDictDF, allCorrDF = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
 
        # do a grid search for C/nu combinations
        scaledMixtureList = np.ndarray.flatten(np.asarray(scaledMixtures.values.tolist())).tolist()
        for cVal in cVals:
            for nuVal in nuVals:
                os.system("echo nu = " +  str(nuVal))
                os.system("echo C = " +  str(cVal))
                coefs, svDict, fileprefix, clfPred, intercept = nuSVR(scaledSigMat, scaledMixtures, nuVal, cVal, droppedGene)
                nuSampName = [srrName + "-" + droppedGene + "-nu=" + str(nuVal) + "-C=" + str(cVal)] 
                
                # iteratively write to the coefficient file
                coefUpdateDF  = pd.DataFrame(data = coefs, columns = scaledSigMat.columns, 
                                index = nuSampName)
                coefDF = pd.concat([coefDF, coefUpdateDF])

                # account for non-uniform length of values in svDict
                # update the key in the dictionary to account for the nu value and save it
                
                # update the keys of the dictionary to the full sample name
                svDict[nuSampName[0]] = svDict.pop(srrName.split("-")[0]) 
                svUpdate = pd.DataFrame(dict([ (k, pd.Series(v)) for k, v in svDict.items()]))
                svDictDF = pd.concat([svDictDF, svUpdate], axis = "columns")

                # get the correlations
                corrDF = pd.DataFrame(index = nuSampName, columns = ["r", "rmse"])
                
                # get the normalized coefs and compute RMSE on predictions
                normCoef = coefs # make copy of coefficients
                normCoef[normCoef < 0] = 0 # set the coefs < 0 -> 0
                normCoef = normCoef / normCoef.sum() # normalize to get the fractional contribs
                predExpr_normCoefs  = np.asarray(np.asmatrix(scaledSigMat).dot(normCoef.T)).squeeze() # use the whole basis matrix
                rmse_nuSVR = np.sqrt(mean_squared_error(scaledMixtureList, predExpr_normCoefs))
                pearsonR_nuSVR, pearsonP_nuSVR = stats.pearsonr(scaledMixtureList, predExpr_normCoefs)
                
                corrDF['rmse'] = [rmse_nuSVR]
                corrDF['r'] = [pearsonR_nuSVR]

                scaledMixtureList = np.ndarray.flatten(np.asarray(scaledMixtures.values.tolist())).tolist()
                clfPred = clfPred.tolist()
                
                allCorrDF = pd.concat([corrDF, allCorrDF])
        
        coefDF = pd.concat([coefDF, allCorrDF], axis = "columns") # append along columns
        aggCoefDF = pd.concat([aggCoefDF, coefDF])
        
        # put everything in its own folder
        svDictDF.to_csv(savePath + fileprefix + "_supportVectors.csv", sep = ",",
                 header = True)                                                                            
    aggCoefDF.to_csv(savePath + coefFile, sep = ",", header = True, index = True) 
    return("foo")


def main(cpmThresh, mixturePath, techrepsrrlist, deconMethod, biologRepName, jackknife = False):
    """
    Main deconvolution function that calls necessary helper functions to perform deconvolution
    Parameters
    ---------
    cpmThresh: int
        CPM values less than this are set to zero

    mixturePath: str
        Filepath to counts matrix of samples for deconvolution 

    techrepsrrlist: list of strings
        List of technical replicates of sample names in mixturePath that will be averaged

    deconMethod: str (either: 'nusvr', 'qp', 'nnls'; case insensitive)
        Deconvolution method to apply on samples

    biologRepName: str
        Name of sample that will be deconvolved 
    
    jackknife: bool
        If True, drop one gene at a time from basis matrix for performing deconvolution

    Returns
    --------
    Nothing, the functions that it calls return/save files as needed
    """
    # read in the bulk RNA counts table for deconvolution 
    thuyCPM = processMixture(mixturePath)
    
    # get the mean of technical replicates and apply cpmThreshold
    thuyCPM = thuyCPM[techrepsrrlist].mean(axis = 1).to_frame() # get the average expression of the gene across the rows
    thuyCPM[thuyCPM >= cpmThresh].dropna()
    
    thuyCPM = thuyCPM[(thuyCPM.T != 0).any()]
    thuyCPM.columns = [biologRepName] 
   
    # read in the basis matrix 
    basisMatrix = pd.read_csv(BASISMATRIX_PATH, sep = "\t", index_col = 0) 
   
    # deconvolution can only be performed on the gene intersection of what's in the mixture and the
    # row space of the basis matrix (e.g. genes) -- reindex both to the gene intersection only 
    intersection = np.intersect1d(thuyCPM.index, basisMatrix.index)
    relThuy = thuyCPM.loc[intersection, :]
    relBM = basisMatrix.loc[intersection, :]
    
    # drop gene duplicates if present 
    val, counts = np.unique(relThuy.index, return_counts = True)
    dups = val[counts > 1].tolist()
    relThuy = relThuy.drop(dups)
    relBM = relBM.drop(dups)
    
    # perform preprocessing about the index to improve runtime performance
    sigMat = relBM
    mixture = relThuy
 
    scaledMixture, scaledSigMat, unscaledSigMat, unscaledMixture  = scale(sigMat, mixture)
    droppedGene = ""

    # it's deconvolution time!
    deconvolve(scaledMixture, scaledSigMat, deconMethod, droppedGene, biologRepName)
    os.system("echo Deconvolution is Complete.")
