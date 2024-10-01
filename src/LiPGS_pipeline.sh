#!/bin/bash
#SBATCH --job-name=LiPGS_pipeline         # Job name
#SBATCH --output=job_%j.out        # Output file (%j expands to job ID)
#SBATCH --error=job_%j.err         # Error file (%j expands to job ID)
#SBATCH -p batch                    # partition 
#SBATCH -N 1                        # number of nodes
#SBATCH -n 1                        
#SBATCH --mem=64G                   # Memory required per node (adjust based on your data)
#SBATCH --time=5:00:00             # Time limit (hours:minutes:seconds)

module load Anaconda3/2024.06-1
source activate PGS

# Run Snakemake, assuming Snakefile is in the same directory
snakemake create_root_dir extract_phenotype split_data

# rm -rf .snakemake/ extract_phenotype split_data --configfile config.yaml 
