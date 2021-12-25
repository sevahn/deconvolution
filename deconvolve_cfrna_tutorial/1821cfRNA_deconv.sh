#!/bin/bash
#SBATCH --job-name=decon1821cfRNA_deconv
#SBATCH --output=decon1821cfRNA_deconv.out
#SBATCH --error=decon1821cfRNA_deconv.err
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --partition=owners,normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000
#SBATCH --time=300:00



source activate snakemake

python3 -c "import deconvolve as deconv; deconv.main(1,  'samples.csv', ['1821'] , 'nuSVR', '1821',  jackknife = False)"
