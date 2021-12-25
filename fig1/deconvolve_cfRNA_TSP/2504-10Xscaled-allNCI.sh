#!/bin/bash
#SBATCH --job-name=decon2504-10Xscaled-allNCI
#SBATCH --output=decon2504-10Xscaled-allNCI.out
#SBATCH --error=decon2504-10Xscaled-allNCI.err
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --partition=quake,normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000
#SBATCH --time=300:00



source activate snakemake

python3 -c "import deconvolve as sev; sev.main(1, '/oak/stanford/groups/quake/sevahn/alzheimers/ad_cpmOnly_postQC_unstranded_FINAL.csv', ['2504'] , 'NNLS', '2504',  jackknife = False)"
python3 -c "import deconvolve as sev; sev.main(1, '/oak/stanford/groups/quake/sevahn/alzheimers/ad_cpmOnly_postQC_unstranded_FINAL.csv', ['2504'] , 'QP', '2504',  jackknife = False)"
python3 -c "import deconvolve as sev; sev.main(1,  '/oak/stanford/groups/quake/sevahn/alzheimers/ad_cpmOnly_postQC_unstranded_FINAL.csv', ['2504'] , 'nuSVR', '2504',  jackknife = False)"