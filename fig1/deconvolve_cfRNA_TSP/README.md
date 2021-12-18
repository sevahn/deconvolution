# nuSVR Deconvolution of Cell Free RNA using Tabula Sapiens v 1.0

To deconvolve samples:
Step 0. Unzip the basis matrix

Step 1. Create the requisite environment using 'cfrna_deconv.yml'

Step 2. The file path to your samples should correspond to a file where the first column is a list of genes and each subsequent column corresponds to a single sample

Step 3. open `sh_1.py` and indicate this file path in the header. Fill in the requested fields, close (`:wq`), and run.
	- if deconvolving a lot of samples, just run `sh_1.py` as its own job, otherwise you can run in the terminal (e.g. `python3 sh_1.py`)
	**the output will be a set of `.sh` files, each corresponding to a sample**

Step 4. If you have your own custom basis matrix, update the header of `deconvolve.py` with its path. Again, the first column should be the list of genes and each subsequent column should correspond to a single sample. The units must match between the two (e.g. both CPM or both TPM, etc) and there must be no log-transformation.
 
Step 5. Launch all sample deconvolution jobs simultaneously: `for i in *.sh ; do sbatch $i ; done`

Step 6. When your jobs complete, run `merge_2.py`, this will write out two files, one of all the combined coefficients and one of all the support vectors corresponding to a given nu/C combination for a  given sample.
 
Couple deconvolution notes:
- Do *not* pass in samples in log-transformed space. For a full proof for why this is inappropriate, please check out this reference: (Zhong, Y. & Liu, Z., Nature Methods 2012)
- The fractions that come out of this program denote the relative fractional contributions of cell type specific RNA for the 23 tissues from which these cell types originate (if using the Tabula Sapiens v1.0 basis matrix) 
- If you are looking for signal from a specific tissue-specific cell type in cfRNA, it is advised to perform signature scoring in conjunction with systems-level deconvolution.
