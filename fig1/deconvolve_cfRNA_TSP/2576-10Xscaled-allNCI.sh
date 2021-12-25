#!/bin/bash
#SBATCH --job-name=decon2576-10Xscaled-allNCI
#SBATCH --output=decon2576-10Xscaled-allNCI.out
#SBATCH --error=decon2576-10Xscaled-allNCI.err
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --partition=quake,normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000
#SBATCH --time=300:00



source activate snakemake

python3 -c "import deconvolve as sev; sev.main(1, '/oak/stanford/groups/quake/sevahn/alzheimers/ad_cpmOnly_postQC_unstranded_FINAL.csv', ['2576'] , 'NNLS', '2576',  jackknife = False)"
python3 -c "import deconvolve as sev; sev.main(1, '/oak/stanford/groups/quake/sevahn/alzheimers/ad_cpmOnly_postQC_unstranded_FINAL.csv', ['2576'] , 'QP', '2576',  jackknife = False)"
python3 -c "import deconvolve as sev; sev.main(1,  '/oak/stanford/groups/quake/sevahn/alzheimers/ad_cpmOnly_postQC_unstranded_FINAL.csv', ['2576'] , 'nuSVR', '2576',  jackknife = False)"