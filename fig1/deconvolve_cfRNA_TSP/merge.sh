#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --output=merge.out
#SBATCH --error=merge.err
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --partition=quake,normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000
#SBATCH --time=300:00



source activate snakemake


#python3 sh.py
python3 merge.py
#python3 predictions.py 
#python3 pearson.py 
