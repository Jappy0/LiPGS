#!/bin/bash
#SBATCH --job-name=PRScs         # Job name
#SBATCH --output=job_%j.out        # Output file (%j expands to job ID)
#SBATCH --error=job_%j.err         # Error file (%j expands to job ID)
#SBATCH -p batch                    # partition 
#SBATCH -N 1                        # number of nodes
#SBATCH -n 1                        
#SBATCH --cpus-per-task=12          # Number of CPU cores per task
#SBATCH --mem=64G                   # Memory required per node (adjust based on your data)
#SBATCH --time=5:00:00             # Time limit (hours:minutes:seconds)

module load Anaconda3/2024.06-1
source activate PGS

python split_data.py --phenotype_file /hpcfs/users/a1236780/Repos/LiPGS/output/ConLiGen_european_phenotypes_binary.txt --plink_prefix /hpcfs/users/a1236780/Repos/LiPGS/ConLiGen/Conligen_filtered_withrsid --output_dir /hpcfs/users/a1236780/Repos/LiPGS/output/data/ --num_groups 5