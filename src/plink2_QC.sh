#!/bin/bash
#SBATCH --job-name=Plink2_QC        # Job name
#SBATCH --output=job_%j.out        # Output file (%j expands to job ID)
#SBATCH --error=job_%j.err         # Error file (%j expands to job ID)
#SBATCH -p batch                    # partition 
#SBATCH -N 1                        # number of nodes
#SBATCH -n 1                    
#SBATCH --mem=128G                   # Memory required per node (adjust based on your data)
#SBATCH --time=50:00:00             # Time limit (hours:minutes:seconds)

module load Anaconda3/2024.06-1
source activate PGS

python plink2_QC.py

echo "PLINK2 QC execution completed."
