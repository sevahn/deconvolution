# Generation of Tabula Sapiens v1 Basis Matrix for Systems Level Deconvolution 

1. run  `dendo.py` on Tabula Sapiens object
2. read in outputted .csv files to `treecutter.ipynb`
3. feed in grouped annotation .csv file from `treecutter.ipynb` to coarsegrain.py
4. pass output of `coarsegrain.py` CIBERSORTx 'Signature Matrix' feature; output is basis matrix (`deconvolve_cfRNA_TSP/tsp_v1_basisMatrix.txt`)
5. viz basis matrix/acquire L2 condition number from Step 4 using the Fig 1E script in the previous directory
6. ready to go for deconvolution - see `deconvolve_cfRNA_TSP` directory for implementation details/general comments.

 
