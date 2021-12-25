#!/bin/bash
#SBATCH --job-name=decon3295-10Xscaled-allNCI
#SBATCH --output=decon3295-10Xscaled-allNCI.out
#SBATCH --error=decon3295-10Xscaled-allNCI.err
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --partition=quake,normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000
#SBATCH --time=300:00



source activate snakemake

python3 -c "import deconvolve as sev; sev.main(1, '/oak/stanford/groups/quake/sevahn/alzheimers/ad_cpmOnly_postQC_unstranded_FINAL.csv', ['3295'] , 'NNLS', '3295',  jackknife = False)"
python3 -c "import deconvolve as sev; sev.main(1, '/oak/stanford/groups/quake/sevahn/alzheimers/ad_cpmOnly_postQC_unstranded_FINAL.csv', ['3295'] , 'QP', '3295',  jackknife = False)"
python3 -c "import deconvolve as sev; sev.main(1,  '/oak/stanford/groups/quake/sevahn/alzheimers/ad_cpmOnly_postQC_unstranded_FINAL.csv', ['3295'] , 'nuSVR', '3295',  jackknife = False)"