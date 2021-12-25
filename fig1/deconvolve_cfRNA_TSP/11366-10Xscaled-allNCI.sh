#!/bin/bash
#SBATCH --job-name=decon11366-10Xscaled-allNCI
#SBATCH --output=decon11366-10Xscaled-allNCI.out
#SBATCH --error=decon11366-10Xscaled-allNCI.err
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --partition=quake,normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000
#SBATCH --time=300:00



source activate snakemake

python3 -c "import deconvolve as sev; sev.main(1, '/oak/stanford/groups/quake/sevahn/alzheimers/ad_cpmOnly_postQC_unstranded_FINAL.csv', ['11366'] , 'NNLS', '11366',  jackknife = False)"
python3 -c "import deconvolve as sev; sev.main(1, '/oak/stanford/groups/quake/sevahn/alzheimers/ad_cpmOnly_postQC_unstranded_FINAL.csv', ['11366'] , 'QP', '11366',  jackknife = False)"
python3 -c "import deconvolve as sev; sev.main(1,  '/oak/stanford/groups/quake/sevahn/alzheimers/ad_cpmOnly_postQC_unstranded_FINAL.csv', ['11366'] , 'nuSVR', '11366',  jackknife = False)"