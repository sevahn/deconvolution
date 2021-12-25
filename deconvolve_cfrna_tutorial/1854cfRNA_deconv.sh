#!/bin/bash
#SBATCH --job-name=decon1854cfRNA_deconv
#SBATCH --output=decon1854cfRNA_deconv.out
#SBATCH --error=decon1854cfRNA_deconv.err
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --partition=owners,normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000
#SBATCH --time=300:00



source activate snakemake

python3 -c "import deconvolve as deconv; deconv.main(1,  'samples.csv', ['1854'] , 'nuSVR', '1854',  jackknife = False)"